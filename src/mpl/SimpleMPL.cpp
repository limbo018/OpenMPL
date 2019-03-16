/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin, Qi Sun
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

/*******************************************
in conflict_num() : uncomment color assertmessage
*******************************************/
#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"

#include <stack>
#include <time.h>
#include <boost/graph/graphviz.hpp>
#include <boost/timer/timer.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>

// only valid when gurobi is available 
#if GUROBI == 1
#include <limbo/algorithms/coloring/ILPColoring.h>
#include <limbo/algorithms/coloring/LPColoring.h>
#include <limbo/algorithms/coloring/MISColoring.h>
#endif
// only valid when lemon cbc api is available 
#if LEMONCBC == 1
#include <limbo/algorithms/coloring/ILPColoringLemonCbc.h>
#endif
#if CSDP == 1
#include <limbo/algorithms/coloring/SDPColoringCsdp.h>
#endif
#include <limbo/algorithms/coloring/BacktrackColoring.h>

#ifdef DEBUG_NONINTEGERS
std::vector<unsigned int> vLP1NonInteger; 
std::vector<unsigned int> vLP1HalfInteger; 
std::vector<unsigned int> vLP2NonInteger; 
std::vector<unsigned int> vLP2HalfInteger; 
std::vector<unsigned int> vLPEndNonInteger; 
std::vector<unsigned int> vLPEndHalfInteger; 
std::vector<unsigned int> vLPNumIter; 
#endif

SIMPLEMPL_BEGIN_NAMESPACE

SimpleMPL::SimpleMPL()
{
    this->reset(true);
}
SimpleMPL::~SimpleMPL()
{
    if (m_db) delete m_db;
}
void SimpleMPL::run(int argc, char** argv)
{
    this->reset(false);
	this->read_cmd(argc, argv);
	this->read_gds();
	this->solve();
#ifdef _DGOUT
	return;
#endif
	this->report();
	this->write_gds();
}
void SimpleMPL::reset(bool init)
{
    // release memory and set to initial value 
    if (!init)
    {
        if (m_db) delete m_db;
        std::vector<uint32_t>().swap(m_vVertexOrder);
        std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
        std::vector<uint32_t>().swap(m_vCompId);
        std::vector<uint32_t>().swap(m_vColorDensity);
        std::vector<std::pair<uint32_t, uint32_t> >().swap(m_vConflict);
    }
    m_db = NULL;
    m_comp_cnt = 0;
}
void SimpleMPL::read_cmd(int argc, char** argv)
{
    // a little bit tricky here 
    // in order to support run-time switch of layoutdb_type
    // we need to construct a dummy ControlParameter for CmdParser
    // then construct actual layoutdb_type according to the option of CmdParser

    // construct a dummy ControlParameter
    ControlParameter tmpParms;
	// read command 
	CmdParser cmd (tmpParms);
	// check options 
    mplAssertMsg(cmd(argc, argv), "failed to parse command");

    // construct actual layout database according to the options 
    if (tmpParms.shape_mode == ShapeModeEnum::RECTANGLE)
        m_db = new LayoutDBRect;
    else 
        m_db = new LayoutDBPolygon;
    // get options from dummy ControlParameter
    tmpParms.swap(m_db->parms);
}

void SimpleMPL::read_gds()
{
    char buf[256];
    mplSPrint(kINFO, buf, "reading input files takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	mplPrint(kINFO, "Reading input file %s\n", m_db->input_gds().c_str());
	// read input gds file 
	GdsReader reader (*m_db);
    mplAssertMsg(reader(m_db->input_gds()), "failed to read %s", m_db->input_gds().c_str());
	// must call initialize after reading 
	m_db->initialize_data();
	// report data 
	m_db->report_data();
}

void SimpleMPL::write_gds()
{
    char buf[256];
    mplSPrint(kINFO, buf, "writing output file takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	if (m_db->output_gds().empty()) 
	{
        mplPrint(kWARN, "Output file not specified, no file generated\n");
		return;
	}
	// write output gds file 
	if (m_db->use_stitch())
	{
		GdsWriter writer;
		mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, StitchRelation, m_mAdjVertex);
	}
	else
	{
		GdsWriter writer;
		mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
}

void SimpleMPL::solve()
{
    // skip if no uncolored layer 
    if (m_db->parms.sUncolorLayer.empty())
        return;

    char buf[256];
    mplSPrint(kINFO, buf, "coloring takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	if (m_db->vPatternBbox.empty())
	{
        mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}

	this->construct_graph();
    if (m_db->simplify_level() > 0) // only perform connected component when enabled 
        this->connected_component();
    else 
    {
        uint32_t vertex_num = m_vVertexOrder.size();
        uint32_t order_id = 0;
        for (uint32_t v = 0; v != vertex_num; ++v)
        {
            m_vCompId[v] = 0;
            m_vVertexOrder[v] = order_id++;
        }
        m_comp_cnt = 1;
    }

	// create bookmark to index the starting position of each component 
	// TODO: if the vBookmark generation is time-consuming, some optimization algorithm may be proposed (set is better than vector)
	std::vector<uint32_t>().swap(vBookmark);
	vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}
#if DEBUG_LIWEI
	std::cout<<"BookMark is: ";
	for (std::vector<char>::const_iterator i = vBookmark.begin(); i != vBookmark.end(); ++i)
    	std::cout << *i << ' ';
	std::cout<<std::endl;
#endif
	if (m_db->use_stitch()) //returns whether use stitch
	{
		this->cal_boundaries();
		this->setVddGnd(); //perhapes we should also consider pshape->getPointNum()==4
		
		clock_t begin = clock();
		
		this->lgSimplification();
		GdsWriter writer;
		writer.write_Simplification(m_db->output_gds() + "_lgSimplification.gds", *m_db, m_vCompId, m_mAdjVertex, in_DG, isVDDGND, true);
		this->projection();		///< vBookmark has already been updated in projection()
		
		clock_t end = clock();
		mplPrint(kINFO, "Projection takes  %f.\n", (double)(end - begin) / CLOCKS_PER_SEC);

		std::vector<uint32_t>().swap(dgCompId);
		dgCompId.assign(m_db->vPatternBbox.size(), 0);
		globalCompId = 1;
		std::cout << "======= Nodes in DG : " << std::count(in_DG.begin(), in_DG.end(), true) << std::endl;
		// update VddGnd information
		this->setVddGnd();
		vdd_multi_comp.resize(m_db->vPatternBbox.size());
		this->dgSimplColoring();
		if (m_db->use_stitch())
		{
			GdsWriter writer;
			writer.write_Simplification(m_db->output_gds() + "_dgSimplification.gds", *m_db, dgCompId, StitchRelation, std::vector<bool>(), isVDDGND, false);
		}
		std::cout << "This graph finally been decomposed into : " << globalCompId << std::endl;

		for (uint32_t i = 0; i < vdd_multi_comp.size(); i++)
		{
			if (isVDDGND[i])
			{
				std::cout << "VDD " << i << " in : ";
				for (std::vector<uint32_t>::iterator it = vdd_multi_comp[i].begin(); it != vdd_multi_comp[i].end(); it++)
					std::cout << *it << " ";
				std::cout << std::endl;
			}
		}
#ifdef _DGOUT
		return;
#endif

	}
	
	mplPrint(kINFO, "Solving %u independent components...\n", m_comp_cnt);
	
	// thread number controled by user option 
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
    for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
    {
        // construct a component 
        std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+vBookmark[comp_id];
        std::vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+vBookmark[comp_id+1] : m_vVertexOrder.end();
        // solve component 
        // pass iterators to save memory 
        this->solve_component(itBgn, itEnd, comp_id);
	}
	
#ifdef DEBUG_NONINTEGERS
    mplPrint(kNONE, "vLP1NonInteger vLP1HalfInteger vLP2NonInteger vLP2HalfInteger vLPEndNonInteger vLPEndHalfInteger vLPNumIter\n"); 
    try 
    {
        // I make it simple, so it may go out of range 
        for (uint32_t i = 0; i < vLPNumIter.size(); ++i)
        {
            mplPrint(kNONE, "%u %u %u %u %u %u %u\n", 
                    vLP1NonInteger.at(i), vLP1HalfInteger.at(i), 
                    vLP2NonInteger.at(i), vLP2HalfInteger.at(i), 
                    vLPEndNonInteger.at(i), vLPEndHalfInteger.at(i), 
                    vLPNumIter.at(i)); 
        }
    }
    catch (std::exception const& e)
    {
        mplPrint(kERROR, "%s\n", e.what()); 
    }
    mplPrint(kNONE, "sum of %lu: %u %u %u %u %u %u %u\n", 
            vLPNumIter.size(), 
            limbo::sum(vLP1NonInteger.begin(), vLP1NonInteger.end()), limbo::sum(vLP1HalfInteger.begin(), vLP1HalfInteger.end()), 
            limbo::sum(vLP2NonInteger.begin(), vLP2NonInteger.end()), limbo::sum(vLP2HalfInteger.begin(), vLP2HalfInteger.end()), 
            limbo::sum(vLPEndNonInteger.begin(), vLPEndNonInteger.end()), limbo::sum(vLPEndHalfInteger.begin(), vLPEndHalfInteger.end()), 
            limbo::sum(vLPNumIter.begin(), vLPNumIter.end()));
#endif
}

void SimpleMPL::lgSimplification()
{
	std::vector<bool>().swap(in_DG);
	std::vector<bool>().swap(articulation_vec);
	in_DG.assign(m_db->vPatternBbox.size(), false);
	articulation_vec.assign(m_db->vPatternBbox.size(), false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		// projection on each component
		this->lgSimplification(itBgn, itEnd, comp_id);
	}
}

void SimpleMPL::lgSimplification(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
#ifdef QDEBUG
	std::cout << "Now in component " << comp_id << " projection.\n";
#endif
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
#if LIWEI_BEFORE
	simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;
#else
	std::cout<<"LG simplify strategy should be HIDE_SMALL_DEGREE"<<std::endl;
#endif
	

	graph_simplification_type gs(dg, m_db->color_num());
	
	gs.set_isVDDGND(vdd_set);
	gs.simplify(simplify_strategy);
	std::vector<bool> projLocal(mGlobal2Local.size(), false);
	std::cout << "mGlobal2Local : " << mGlobal2Local.size() << std::endl;
	std::cout << "hidden_vertices() : " << gs.hidden_vertices().size() << std::endl;

	std::vector<vertex_descriptor> all_articulations;
	gs.get_articulations(all_articulations);
	std::set<vertex_descriptor> temp_set(all_articulations.begin(), all_articulations.end());

	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<vertex_descriptor> vSimpl2Orig;
		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		for (uint32_t i = 0; i < vSimpl2Orig.size(); i++)
			projLocal[vSimpl2Orig[i]] = true;
	}

	for (std::map<uint32_t, uint32_t>::iterator it = mGlobal2Local.begin(); it != mGlobal2Local.end(); it++)
	{
		if (projLocal[it->second])
			in_DG[it->first] = true;
		if (temp_set.find(it->second) != temp_set.end())
			articulation_vec[it->first] = true;
	}
}

void SimpleMPL::dgSimplColoring()
{
	std::vector<std::vector<uint32_t> > comp2subcomp(m_comp_cnt);
	std::vector<std::map<vertex_descriptor, std::set<uint32_t> > > m_mArtiPoints(m_comp_cnt);
	std::vector<std::vector<vertex_descriptor> > m_mCompVertex;
	std::vector<std::vector<vertex_descriptor> > HiddenVertices;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		std::vector<std::vector<vertex_descriptor> > comp_vertex;
		std::map<vertex_descriptor, std::set<uint32_t> >& arti_point = m_mArtiPoints[comp_id];
		// Simplification on each component
		this->dgSimplColoring(itBgn, itEnd, comp_id, arti_point, comp_vertex);
	}
#ifdef _DGOUT
	return;
#endif 
}


/*
for (std::vector<vertex_descriptor>::iterator it = sub_vSimpl2Orig.begin(); it != sub_vSimpl2Orig.end(); it++)
{
	// std::cout << "*it : " << *it << std::endl;
	// std::cout << "vSimpl2Orig.at(*it) : " << vSimpl2Orig.at(*it) << std::endl;
	vertex_descriptor curpid = mLocal2Global.at(vSimpl2Orig.at(*it));
	// std::cout << curpid << std::endl;
				
	dgCompId[curpid] = globalCompId;
	if (isVDDGND[curpid])
	{
		vdd_multi_comp[curpid].push_back(globalCompId);
	}
}
globalCompId++;
*/

void SimpleMPL::dgSimplColoring(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::map<vertex_descriptor, std::set<uint32_t> >& m_ArtiPoint, std::vector<std::vector<vertex_descriptor> >& m_CompVertex)
{
	std::cout << "\n\n\n**************************************************************\n";
	std::cout << "Debug| dgSimplification comp " << comp_id << std::endl;
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
	simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

	graph_simplification_type gs(dg, m_db->color_num());
	gs.precolor(vColor.begin(), vColor.end());

	std::vector<uint32_t> mLocal2Global(mGlobal2Local.size());
	for (std::map<uint32_t, uint32_t>::iterator it = mGlobal2Local.begin(); it != mGlobal2Local.end(); it++)
		mLocal2Global[it->second] = it->first;

	gs.set_isVDDGND(vdd_set);
	gs.simplify(simplify_strategy);

	///< big component : input big component (the only several big components in original input gds file).
	///< small component : small component divided from big component after conducting IVR and biconnected decomposition, may contains one or two VDDGNDs, is corresponding one row in big component
	///< sub small component : sub small component divided from small component after only conducting biconnected decomposition, divide one row into sevel small sub components
	
	std::vector<std::vector<int8_t> > big_comp_color(gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSmall2Big (gs.num_component());
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

	double total_acc_obj_value = 0;

	for (uint32_t small_comp_id = 0; small_comp_id < gs.num_component(); ++small_comp_id)
	{
		std::cout << "\n==================================================\n";
		std::cout << "Debug| dgSimplification sub_comp " << small_comp_id << std::endl;
		graph_type small_g;
		
		std::vector<int8_t>& small_comp_color = big_comp_color[small_comp_id];	///< used to store coloring results of current small component
		
		std::vector<vertex_descriptor>& vSimplSmall2OriBig = mSmall2Big[small_comp_id];		///< used to store the simplificaiton mapping information from current small component to original big component

		gs.simplified_graph_component(small_comp_id, small_g, vSimplSmall2OriBig);

		uint32_t temp = 0;
		std::set<vertex_descriptor> small_vdd_set;
		std::vector<int8_t> small_vColor(num_vertices(small_g), -1);
		///< construct small vColor
		///< also get vdd_set for current sub component
		for (std::vector<vertex_descriptor>::iterator it = vSimplSmall2OriBig.begin(); it != vSimplSmall2OriBig.end(); it++, temp++)
		{
			if (gs.whether_VDDGND(*it))
				small_vdd_set.insert(temp);
			small_vColor[temp] = vColor[*it];
		}

		///< ***************************************************///
		///< Now simplify current small component
		///< **************************************************///

		graph_simplification_type simplified_small_g(small_g, m_db->color_num());	///< construct small graph simplification object
		///< merge VDD nodes in sg, and conduct biconnected devision
		simplified_small_g.precolor(small_vColor.begin(), small_vColor.end());
		simplified_small_g.mergeVDD(small_vdd_set);
		simplified_small_g.simplify(graph_simplification_type::BICONNECTED_COMPONENT);
		std::vector<std::vector<int8_t> > small_comp_color_2 (simplified_small_g.num_component());	///< used to store the coloring results of sub small components, and recover small_comp_color;
		std::vector<std::vector<vertex_descriptor> > sub_small2small(simplified_small_g.num_component());	///< store the mapping information from sub small components to current small component
		
		double acc_obj_value = 0;

		for (uint32_t ss_comp_id = 0; ss_comp_id < simplified_small_g.num_component(); ++ss_comp_id)
		{
			std::cout << "\n===============================\n";
			std::cout << "Debug| dgSimplification ss_comp " << ss_comp_id << std::endl;
			std::cout << "Debug| dgSimplification with nodes : ";
			graph_type sub_small_simpl_sg;
			
			std::vector<vertex_descriptor>& sub_vSimpl2OrigSmall = sub_small2small[ss_comp_id];
			simplified_small_g.simplified_graph_component(ss_comp_id, sub_small_simpl_sg, sub_vSimpl2OrigSmall);
			
			std::vector<int8_t>& sub_small_comp_color = small_comp_color_2[ss_comp_id];	///< used to store 
			sub_small_comp_color.assign(num_vertices(sub_small_simpl_sg), -1);

			typedef lac::Coloring<graph_type> coloring_solver_type;
			coloring_solver_type* pcs = create_coloring_solver(sub_small_simpl_sg);

			boost::graph_traits<graph_type>::vertex_iterator vi, vie;
			for (boost::tie(vi, vie) = vertices(sub_small_simpl_sg); vi != vie; vi++)
			{
				vertex_descriptor v = *vi;
				int8_t color = small_vColor[sub_vSimpl2OrigSmall[v]];
				if (color >= 0 && color < m_db->color_num())
				{
					pcs->precolor(v, color);
					sub_small_comp_color[v] = color; ///< necessary for merge_K4 coloring
				}
			}
			double obj_value1 = (*pcs)();
			double obj_value2 = std::numeric_limits<double>::max();
#ifndef DEBUG_NONINTEGERS
			if (obj_value1 >= 1 && boost::num_vertices(sub_small_simpl_sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::SDP_CSDP) && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // merge K4 is not performed
				obj_value2 = merge_K4_coloring(sub_small_simpl_sg, sub_small_comp_color);
#endif
			if (obj_value1 < obj_value2)
			{
				acc_obj_value += obj_value1;
				///< collect  coloring results from simplified graph
				for (boost::tie(vi, vie) = vertices(sub_small_simpl_sg); vi != vie; vi++)
				{
					vertex_descriptor v = *vi;
					int8_t color = pcs->color(v);
					mplAssert(color >= 0 && color < m_db->color_num());
					sub_small_comp_color[v] = color;
				}
			}
			else // no need to update sub_small_comp_color, since it has already been updated by merge_K4_coloring
				acc_obj_value += obj_value2;
		}
	}
}

uint32_t SimpleMPL::merge_K4_coloring(SimpleMPL::graph_type const& dg, std::vector<int8_t>& vColor) const
{
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs(dg, m_db->color_num());
	gs.precolor(vColor.begin(), vColor.end());
	if (m_db->color_num() == 3)
		gs.max_merge_level(3);
	else if (m_db->color_num() == 4)
		gs.max_merge_level(2);
	gs.simplify(graph_simplification_type::MERGE_SUBK4);
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

	std::vector<std::vector<int8_t> > mSubColor(gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig(gs.num_component());
	double acc_obj_value = 0;

	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];
		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		vSubColor.assign(num_vertices(sg), -1);

		typedef lac::Coloring<graph_type> coloring_solver_type;
		coloring_solver_type* pcs = create_coloring_solver(sg);

		boost::graph_traits<graph_type>::vertex_iterator vi, vie;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = vColor[vSimpl2Orig[v]];
			if (color >= 0 && color < m_db->color_num())
			{
				pcs->precolor(v, color);
				vSubColor[v] = color;
			}
		}
		double obj_value1 = (*pcs)();

		acc_obj_value += obj_value1;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; vi++)
		{
			vertex_descriptor v = *vi;
			int8_t color = pcs->color(v);
			vSubColor[v] = color;
		}

		delete pcs;
	}

	gs.recover(vColor, mSubColor, mSimpl2Orig);

	return acc_obj_value;
}

void SimpleMPL::setVddGnd()
{
	std::vector<bool>().swap(isVDDGND);
	int threshold = 6000;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	isVDDGND.assign(vertex_num, false);
	for (uint32_t i = 0; i < vertex_num; i++)
	{
#if LIWEI_BEFORE
		if (boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) >= threshold)
		{
			isVDDGND[i] = true;
			// m_db->vPatternBbox[i]->color(m_db->color_num() - 1);
		}
#else
		if (boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) >= 0.64*(boundaries[1]-boundaries[0])||\
		boost::polygon::delta(*m_db->vPatternBbox[i], gtl::VERTICAL) >= 0.64*(boundaries[3]-boundaries[2]))
		{
			isVDDGND[i] = true;
			// m_db->vPatternBbox[i]->color(m_db->color_num() - 1);
		}
		if(boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) == 6000 ){
			isVDDGND[i] = false;
			std::cout<<"dirty coding for oracle cases"<<std::endl;
		}
#endif
#if 0
		else
		{
			coordinate_difference width = boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL);
			coordinate_difference height = boost::polygon::delta(*m_db->vPatternBbox[i], gtl::VERTICAL);
		}
#endif
	}
	return;
}

void SimpleMPL::cal_boundaries(){
	std::vector<uint32_t>().swap(boundaries);
	assert(boundaries.empty());
	uint32_t vertex_num = m_db->vPatternBbox.size();
	uint32_t tmp_bound[4] = {0};
	for (uint32_t i = 0; i < vertex_num; i++){
		tmp_bound[0] = std::max(gtl::xl(*vPatternBbox[i]),tmp_bound[0]);
		tmp_bound[1] = std::max(gtl::xh(*vPatternBbox[i]),tmp_bound[1]);
		tmp_bound[2] = std::max(gtl::yl(*vPatternBbox[i]),tmp_bound[2]);
		tmp_bound[3] = std::max(gtl::yh(*vPatternBbox[i]),tmp_bound[3]);
	}
	boundaries.push_back(tmp_bound[0]);
	boundaries.push_back(tmp_bound[1]);
	boundaries.push_back(tmp_bound[2]);
	boundaries.push_back(tmp_bound[3]);
#if DEBUG_LIWEI
	std::cout<<"boundries of layout is calculated and the value is"<<boundaries[0]<<','\
																	 boundaries[1]<<','\
																	 boundaries[2]<<','\
																	 boundaries[3]<<','\
	<<std::endl;
#endif
	return;
}

void SimpleMPL::projection()
{
	uint32_t vertex_num = m_db->vPatternBbox.size();
	std::vector<rectangle_pointer_type> rect_vec = m_db->polyrect_patterns();	///< original rectangle list
	std::vector<uint32_t> Poly_Rect_begin = m_db->PolyRectBgnLoc();				///< original polygons mapping to rectangles
	std::vector<uint32_t> Poly_Rect_end;
	Poly_Rect_end.resize(vertex_num);

	for (uint32_t i = 0; i < vertex_num - 1; i++)
		Poly_Rect_end[i] = Poly_Rect_begin[i + 1] - 1;
	Poly_Rect_end[vertex_num - 1] = rect_vec.size() - 1;

	std::vector<uint32_t> new_vCompId_vec;
	std::vector<uint32_t> new_vertex_order;
	std::vector<rectangle_pointer_type> new_rect_vec;				///< store the new rectangles
	std::vector<uint32_t> new_Rect2ParentPoly;						///< map from new rectangles to new parent polygons
	std::vector<std::vector<uint32_t> >().swap(ori2new_polygon);	///< map from original polygons to new polygons
	std::vector<uint32_t>().swap(new2ori_polygon);					///< map from new polygons to old polygons
	ori2new_polygon.resize(vertex_num);
	
	uint32_t new_polygon_id = -1;						///< new polygon id
	uint32_t new_rectangle_id = -1;						///< new rectangle id

	for (uint32_t i = 0; i < vertex_num; i++)
	{
		uint32_t v = m_vVertexOrder[i];
		uint32_t comp_id = m_vCompId[v];

		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[v];
		uint32_t pid = pPattern->pattern_id();

		uint32_t start_idx = Poly_Rect_begin[pid];
		uint32_t end_idx = Poly_Rect_end[pid];
		
		std::vector<uint32_t>& nei_vec = m_mAdjVertex[pid];
#if DEBUG_LIWEI
		std::cout<<"pid is"<<pid<<", v is"<<v<<".They should be equal"<<std::endl;
#endif
		if (in_DG[pid] && articulation_vec[pid] == false && isVDDGND[pid] == false)
		{
			std::vector<rectangle_pointer_type> poss_nei_vec;
			for (std::vector<uint32_t>::iterator it = nei_vec.begin(); it != nei_vec.end(); it++)
			{
				if (in_DG[*it] == false) continue;
				uint32_t s_idx = Poly_Rect_begin[*it];
				uint32_t e_idx = Poly_Rect_end[*it];
				poss_nei_vec.insert(poss_nei_vec.end(), rect_vec.begin() + s_idx, rect_vec.begin() + e_idx + 1);
			} ///< for nei_vec

			std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_split;	///< the generated rectangles and their parent original rectangle id
			///< traverse all rectangles in current polygon
			for (uint32_t j = start_idx; j <= end_idx; j++)
			{
				rectangle_type rect(*rect_vec[j]);
				std::vector<rectangle_pointer_type> split;
				splitRectangle(rect, split, poss_nei_vec);
#if DEBUG_LIWEI
				std::cout<<"rect pattern id is, this one is to check whether the rect id in same pattren is same (it should be same)"<<rect.pattern_id()<<std::endl;
#endif
				for (uint32_t i = 0; i < split.size(); i++)
					rect_split.push_back(std::make_pair(split[i], rect.pattern_id()));
			} ///< for

			uint32_t pivot = new_polygon_id;
			std::vector<uint32_t> new_polygon_id_list;
			// reconstruct polygons, also generate stitch relationships
			reconstruct_polygon(new_polygon_id, new_polygon_id_list, rect_split);

			assert(new_polygon_id_list.size() == rect_split.size());

			for (uint32_t i = 0; i < rect_split.size(); i++)
			{
				rect_split[i].first->pattern_id(++new_rectangle_id);
				new_rect_vec.push_back(rect_split[i].first);
				new_Rect2ParentPoly.push_back(new_polygon_id_list[i]);
				
				// insert new polygon
				if (pivot != new_polygon_id_list[i])
				{
					pivot = new_polygon_id_list[i];
					new2ori_polygon.push_back(pid);
					ori2new_polygon[pid].push_back(pivot);
					new_vertex_order.push_back(pivot);
					new_vCompId_vec.push_back(comp_id);
				}
			} ///< for rect_list
		} ///< if, stitch generation core
		else
		{
			new_polygon_id += 1;
			for (uint32_t j = start_idx; j <= end_idx; j++)
			{
				rectangle_pointer_type rect = rect_vec[j];
				rect->pattern_id(++new_rectangle_id);
				if (m_db->gen_stitch())
					rect->color(new_polygon_id % 7);
				new_rect_vec.push_back(rect);
				new_Rect2ParentPoly.push_back(new_polygon_id);
			}
			new2ori_polygon.push_back(pid);
			new_vertex_order.push_back(new_polygon_id);
			new_vCompId_vec.push_back(comp_id);
			ori2new_polygon[pid].push_back(new_polygon_id);
			StitchRelation.push_back(std::vector<uint32_t>());
		} ///< else (in_DG)
	} ///< for all vertices

	m_db->refresh(new_rect_vec, new_Rect2ParentPoly);

	updateConflictRelation();


	std::vector<uint32_t>().swap(m_vCompId);
	m_vCompId.swap(new_vCompId_vec);

	std::vector<uint32_t>().swap(m_vVertexOrder);
	m_vVertexOrder.swap(new_vertex_order);

	std::vector<uint32_t>().swap(vBookmark);
	vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

	return;
}

void SimpleMPL::updateConflictRelation()
{
	// now update adjacency list, we still need original adjacency list.
	// here we use set to remove duplicated neighboring polygons
	std::vector<std::set<uint32_t> >new_mAdjVertex;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	new_mAdjVertex.resize(m_db->vPatternBbox.size());
	uint32_t edge_num = 0;
	/*
	// generate new index information
	std::vector<uint32_t> new_poly_rect_begin = m_db->PolyRectBgnLoc();
	std::vector<uint32_t> new_poly_rect_end;
	
	assert(vertex_num == m_db->vPatternBbox.size());

	new_poly_rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		new_poly_rect_end[i] = new_poly_rect_begin[i + 1] - 1;
	new_poly_rect_end[vertex_num - 1] = new_rect_vec.size() - 1;
	*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif
	for (uint32_t v = 0; v < vertex_num; ++v)
	{
		rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
		std::set<uint32_t>& vAdjVertex = new_mAdjVertex[v];

		// find patterns connected with pPattern 
		// query tPatternBbox in m_db
		rectangle_type rect(*pPattern);
		// bloat pPattern with minimum coloring distance 
		gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
		gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
		for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
			itq != m_db->tPatternBbox.qend(); ++itq)
		{
			rectangle_pointer_type const& pAdjPattern = *itq;
			// Stitch relationships have already generated. If these two polygons have the same original polygon, ignore them.
			if (new2ori_polygon[pAdjPattern->pattern_id()] == new2ori_polygon[pPattern->pattern_id()])
				continue;
			if (pAdjPattern != pPattern) // skip pPattern itself 
			{
				mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
				// we consider euclidean distance
				// use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
				coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
				if (distance < m_db->coloring_distance)
					vAdjVertex.insert(pAdjPattern->pattern_id());
			}
		}
		vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
		edge_num += vAdjVertex.size();
	}
	std::cout << "After stitch insertion, it has " << edge_num << " edges.\n";

	///< update m_mAdjVertex
	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex.resize(new_mAdjVertex.size());
	for (uint32_t i = 0; i < new_mAdjVertex.size(); i++)
		m_mAdjVertex[i].insert(m_mAdjVertex[i].end(), new_mAdjVertex[i].begin(), new_mAdjVertex[i].end());

	/*
	// traverse all the original polygons to get the original neighbor list
	for (uint32_t i = 0, ie = m_mAdjVertex.size(); i < ie; i++)
	{
		// list all the possible polygon neighbors
		std::vector<uint32_t> poss_nei;
		bool split_flag = false;
		for (uint32_t j = 0, je = m_mAdjVertex[i].size(); j < je; j++)
		{
			uint32_t nei_id = m_mAdjVertex[i][j];
			for (uint32_t it = 0; it < ori2new_polygon[nei_id].size(); it++)
				poss_nei.push_back(ori2new_polygon[nei_id][it]);
		}
		// traverse all the newly-generated nodes in polygon i
		for (std::vector<uint32_t>::iterator it = ori2new_polygon[i].begin(); it != ori2new_polygon[i].end(); it++)
		{
			uint32_t start_idx = new_poly_rect_begin[*it];
			uint32_t end_idx = new_poly_rect_end[*it];
			for (std::vector<uint32_t>::iterator nei_poly = poss_nei.begin(); nei_poly != poss_nei.end(); nei_poly++)
			{
				if (*nei_poly > *it)
					continue;
				for (uint32_t now_rect = start_idx; now_rect <= end_idx; now_rect++)
				{
					uint32_t nei_start_idx = new_poly_rect_begin[*nei_poly];
					uint32_t nei_end_idx = new_poly_rect_end[*nei_poly];
					// traverse all rectangles in current newly-generated neighbor polygon
					for (uint32_t nei_rect = nei_start_idx; nei_rect <= nei_end_idx; nei_rect++)
					{
						uint32_t nei_rect_pid = new_rect_vec[nei_rect]->pattern_id();
						coordinate_difference distance = boost::geometry::distance(*new_rect_vec[now_rect], *new_rect_vec[nei_rect_pid]);
						if (distance < m_db->coloring_distance)
						{
							new_mAdjVertex[*it].insert(*nei_poly);
							new_mAdjVertex[*nei_poly].insert(*it);
							break;
						}// if
					} // for
				} // for
			} // for
		} // for
		
	} // for
*/
}

void SimpleMPL::reconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t>& new_polygon_id_list, std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list)
{
	std::vector<std::vector<uint32_t> > stitch_list;
	uint32_t start = polygon_id + 1;
	std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_temp;
	std::vector<uint32_t> new_polygon_id_temp;

	uint32_t rect_list_size = rect_list.size();
	std::vector<bool> visited(rect_list_size, false);

	std::vector<uint32_t>().swap(new_polygon_id_list);
	new_polygon_id_list.assign(rect_list_size, std::numeric_limits<uint32_t>::max());

	// traverse all the rectangles to generate new IDs
	for (uint32_t i = 0; i< rect_list_size; i++)
	{
		if (new_polygon_id_list[i] == std::numeric_limits<uint32_t>::max())
		{
			polygon_id += 1; //TODO: here maybe a bug, += 1 should be behind
			new_polygon_id_list[i] = polygon_id;
		}
		visited[i] = true;
		for (uint32_t j = i + 1; j < rect_list_size; j++)
		{
			if (visited[j] == false && rect_list[i].second != rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				visited[j] = true;
				new_polygon_id_list[j] = new_polygon_id_list[i];
			}
		}
	}

	// for rectangles in the same original rectangle, we need to introduce a stitch between them, and record the stitch relationships.
	stitch_list.resize(polygon_id - start + 1);
	for (uint32_t i = 0; i < rect_list.size(); i++)
	{
		for (uint32_t j = i + 1; j < rect_list.size(); j++)
		{
			if (rect_list[i].second == rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				stitch_list[new_polygon_id_list[i] - start].push_back(new_polygon_id_list[j]);
				stitch_list[new_polygon_id_list[j] - start].push_back(new_polygon_id_list[i]);
			}
		}
	}
	// rectangles in the same polygon should be abutting
	for (uint32_t i = start, ie = polygon_id; i <= ie; ++i)
	{
		for (uint32_t j = 0; j < new_polygon_id_list.size(); ++j)
		{
			if (new_polygon_id_list[j] == i)
			{
				new_polygon_id_temp.push_back(i);
				rect_temp.push_back(rect_list[j]);
			}
		}
	}
	rect_list.swap(rect_temp);
	new_polygon_id_list.swap(new_polygon_id_temp);
	StitchRelation.insert(StitchRelation.end(), stitch_list.begin(), stitch_list.end());
}

void SimpleMPL::splitRectangle(rectangle_type & pRect, std::vector<rectangle_pointer_type>& split, std::vector<rectangle_pointer_type> nei_Vec)
{
	bool hor = gtl::delta(pRect, gtl::HORIZONTAL) >= gtl::delta(pRect, gtl::VERTICAL);
	std::vector<coordinate_type> vstitches;
	coordinate_difference width = gtl::xh(pRect) - gtl::xl(pRect);

	coordinate_type lower_boundary;
	coordinate_type upper_boundary;

	std::set<coordinate_type> vset;
	std::vector<coordinate_type> vPossibleStitches;
	if (hor)
	{
		lower_boundary = gtl::xl(pRect);
		upper_boundary = gtl::xh(pRect);
	}
	else {
		lower_boundary = gtl::yl(pRect);
		upper_boundary = gtl::yh(pRect);
	}
	vset.insert(lower_boundary);
	vset.insert(upper_boundary);

	std::vector<rectangle_type> vInterSect;
	// generate intersections
	for (std::vector<rectangle_pointer_type>::iterator it = nei_Vec.begin(); it != nei_Vec.end(); it++)
	{
		coordinate_difference distance = boost::geometry::distance(pRect, *(*it));
		if (distance >= m_db->coloring_distance) continue;
		rectangle_type temp(*(*it));
		gtl::bloat(temp, gtl::HORIZONTAL, m_db->coloring_distance);
		gtl::bloat(temp, gtl::VERTICAL, m_db->coloring_distance);
		
		boost::polygon::intersect(temp, pRect, true);
		vInterSect.push_back(temp);
		if (hor)
		{
			vset.insert(gtl::xl(temp));
			vset.insert(gtl::xh(temp));
		}
		else
		{
			vset.insert(gtl::yl(temp));
			vset.insert(gtl::yh(temp));
		}
	}

	for (std::set<coordinate_type>::iterator it = vset.begin(); it != vset.end(); it++)
		vPossibleStitches.push_back(*it);
	std::sort(vPossibleStitches.begin(), vPossibleStitches.end());

	GenerateStitchPosition(pRect, vInterSect, vPossibleStitches, vstitches);

#ifdef STITCH_FILTER
	// check the stitch positions' legalities
	// if the position is very colse to the rectangle's boundary, it's illegal.
	coordinate_type threshold = 20;
	std::vector<coordinate_type> temp;
	for (std::vector<coordinate_type>::iterator it = vstitches.begin(); it != vstitches.end(); it++)
	{
		coordinate_type dis_low = std::abs(*it - lower_boundary);
		coordinate_type dis_up = std::abs(*it - upper_boundary);
		if (dis_low >= threshold && dis_up >= threshold)
			temp.push_back(*it);
	}
	std::vector<coordinate_type>().swap(vstitches);
	vstitches.swap(temp);
#endif

#ifdef STITCH_LEGAL
	if (hor)
	{
		for(int )
	}
#endif
	// split rectangles according to the stitch positions.
	if (vstitches.size() <= 0)
	{
		rectangle_pointer_type new_Pattern = new rectangle_type(pRect);
		new_Pattern->color(pRect.color());
		split.push_back(new_Pattern);
	}
	else
	{
		if (hor)
		{
			vstitches.insert(vstitches.begin(), gtl::xl(pRect));
			vstitches.push_back(gtl::xh(pRect));
			for (uint32_t j = 0; j < vstitches.size() - 1; j++)
			{
				rectangle_pointer_type new_Pattern = new rectangle_type(
					vstitches[j], gtl::yl(pRect), vstitches[j + 1], gtl::yh(pRect));
				new_Pattern->color(pRect.color());
				split.push_back(new_Pattern);
			}
		}
		else
		{
			vstitches.insert(vstitches.begin(), gtl::yl(pRect));
			vstitches.push_back(gtl::yh(pRect));
			for (uint32_t j = 0; j < vstitches.size() - 1; j++)
			{
				rectangle_pointer_type new_Pattern = new rectangle_type(
					gtl::xl(pRect), vstitches[j], gtl::xh(pRect), vstitches[j + 1]);
				new_Pattern->color(pRect.color());
				split.push_back(new_Pattern);
			}
		}
	}
	mplAssert(split.size() > 0);

	return;
}

void SimpleMPL::GenerateStitchPosition(const rectangle_type pRect, std::vector<rectangle_type> vInterSect,
	std::vector<coordinate_type>& vPossibleStitches, std::vector<coordinate_type>& vstitches)
{
	bool ishor = gtl::delta(pRect, gtl::HORIZONTAL) >= gtl::delta(pRect, gtl::VERTICAL);
	coordinate_type lower, upper;
	if (ishor) {
		lower = gtl::xl(pRect);
		upper = gtl::xh(pRect);
	}
	else {
		lower = gtl::yl(pRect);
		upper = gtl::yh(pRect);
	}
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, int > > vStages;
	for (uint32_t i = 1; i < vPossibleStitches.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(vPossibleStitches[i - 1], vPossibleStitches[i]), 0));
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		coordinate_type left = vStages[i].first.first;
		coordinate_type right = vStages[i].first.second;
		for (uint32_t j = 0; j < vInterSect.size(); j++)
		{
			if (ishor)
			{
				if (left < gtl::xl(vInterSect[j])) continue;
				if (right > gtl::xh(vInterSect[j])) continue;
				vStages[i].second++;
			}
			else
			{
				if (left < gtl::yl(vInterSect[j])) continue;
				if (right > gtl::yh(vInterSect[j])) continue;
				vStages[i].second++;
			}
		}
	}
	if (vStages.size() <= 0) return;
	if (vStages[0].second != 0)
		vStages.insert(vStages.begin(), std::make_pair(std::make_pair(lower, lower), 0));
	if (vStages[vStages.size() - 1].second != 0)
		vStages.push_back(std::make_pair(std::make_pair(upper, upper), 0));
#ifdef _DEBUG_QI
	std::cout << "DEBUG_QI| vStages = ";
	for (uint32_t i = 0; i < vStages.size(); i++)
		std::cout << vStages[i].second << " ";
	std::cout << std::endl;
#endif
	std::vector<int32_t> vzeroids;
	for (int32_t i = 0; i < vStages.size(); i++)
	{
		if (vStages[i].second > 0) continue;
		vzeroids.push_back(i);
	}

	vstitches.clear();
	for (int32_t i = 0; i < vzeroids.size() - 1; i++)
	{
		int32_t pos1 = vzeroids[i];
		int32_t pos2 = vzeroids[i + 1];
		if (0 == i) assert(0 == pos1);
		else if (2 == pos1) // remove the useless stitch
		{
			bool bFind = true;
			if (vStages.size() < 5) bFind = false;
			else if (i != 1) bFind = false;
			else if (1 != vStages[1].second) bFind = false;
			else if (1 != vStages[3].second) bFind = false;
			else if (0 != vStages[4].second) bFind = false;
			if (bFind == true) continue;
			coordinate_type pos = (vStages[2].first.first + vStages[2].first.second) / 2;
			vstitches.push_back(pos);
		}
		else if (pos1 == vStages.size() - 3)
		{
			assert(i == vzeroids.size() - 2);
			bool bFind = true;
			int zsize = vzeroids.size();
			if (vStages.size() < 5) bFind = false;
			else if (i != zsize - 2) bFind = false;
			else if (1 != vStages[pos1 + 1].second) bFind = false;
			else if (1 != vStages[pos1 - 1].second) bFind = false;
			else if (0 != vStages[pos1 - 2].second) bFind = false;
			if (bFind == true) continue;
			coordinate_type pos = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}
		else
		{
			coordinate_type pos = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}

		double maxValue = 0.9;
		int32_t posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		for (int32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			int mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
			int diff = std::abs(vStages[i + 1].second - vStages[i - 1].second);
			double value = (double)mind + (double)diff*0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			coordinate_type pos = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
			vstitches.push_back(pos);
		}
	}
	std::sort(vstitches.begin(), vstitches.end());

#ifdef _DEBUG_PROJECTION
	if (lower == -2275 && upper == -2090)
	{
		std::cout << "\nDEBUG_PROJECTION| output the vzeroids : ";
		for (int i = 0; i < vzeroids.size(); i++)
			std::cout << vzeroids[i] << " ";
		std::cout << std::endl;
	}
	std::cout << "DEBUG_PROJECTION| output the vStages : ";
	for (uint32_t i = 0; i < vStages.size(); i++)
		std::cout << vStages[i].second << "[ " << vStages[i].first.first << "," << vStages[i].first.second << "] ;    ";
	std::cout << std::endl;

	std::cout << "DEBUG_PROJECTION| output the vstitches : ";
	for (uint32_t i = 0; i < vstitches.size(); i++)
		std::cout << vstitches[i] << " ";
	std::cout << std::endl;
#endif
}


void SimpleMPL::report() const 
{
    mplPrint(kINFO, "Conflict number = %u\n", conflict_num());
	for (int32_t i = 0, ie = m_db->color_num(); i != ie; ++i)
        mplPrint(kINFO, "Color %d density = %u\n", i, m_vColorDensity[i]);
	int count = 0;
	for (uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
	{
		if (m_db->vPatternBbox[i]->color() >= m_db->color_num())
			count++;
	}
	mplPrint(kINFO, "Invalid color number : %d\n", count);
}

void SimpleMPL::construct_graph()
{
	mplPrint(kINFO, "Constructing graph for %lu patterns...\n", m_db->vPatternBbox.size());
	// construct vertices 
	// assume vertices start from 0 and end with vertex_num-1
	uint32_t vertex_num = m_db->vPatternBbox.size();
	uint32_t edge_num = 0;
	m_vVertexOrder.resize(vertex_num, std::numeric_limits<uint32_t>::max()); 

	// construct edges 
	m_mAdjVertex.resize(vertex_num);
	if (m_db->hPath.empty()) // construct from distance 
        edge_num = construct_graph_from_distance(vertex_num);
	else // construct from conflict edges in hPath
        edge_num = construct_graph_from_paths(vertex_num);
	m_vCompId.resize(vertex_num, std::numeric_limits<uint32_t>::max());
	m_vColorDensity.assign(m_db->color_num(), 0);
	m_vConflict.clear();
    edge_num = edge_num>>1;

	// report statistics 
	mplPrint(kINFO, "%u vertices, %u edges\n", vertex_num, edge_num);
}

uint32_t SimpleMPL::construct_graph_from_distance(uint32_t vertex_num)
{
    uint32_t edge_num = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif
    for (uint32_t v = 0; v < vertex_num; ++v)
    {
        rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];

        // find patterns connected with pPattern 
        // query tPatternBbox in m_db
        rectangle_type rect (*pPattern);
        // bloat pPattern with minimum coloring distance 
        gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
        gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
        for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
                itq != m_db->tPatternBbox.qend(); ++itq)
        {
            rectangle_pointer_type const& pAdjPattern = *itq;
            if (pAdjPattern != pPattern) // skip pPattern itself 
            {
                mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
                // we consider euclidean distance
                // use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
                coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
                if (distance < m_db->coloring_distance)
                    vAdjVertex.push_back(pAdjPattern->pattern_id());
            }
        }
        vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
        // parallel with omp reduction here 
        edge_num += vAdjVertex.size();
    }
    return edge_num;
}

uint32_t SimpleMPL::construct_graph_from_paths(uint32_t vertex_num)
{
    // at the same time, estimate a coloring distance from conflict edges 
    m_db->coloring_distance = 0;
    for (set<int32_t>::const_iterator itLayer = m_db->parms.sPathLayer.begin(), itLayerE = m_db->parms.sPathLayer.end(); itLayer != itLayerE; ++itLayer)
    {
        if (!m_db->hPath.count(*itLayer)) continue;

        std::vector<path_type> const& vPath = m_db->hPath[*itLayer];

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
        for (uint32_t i = 0; i < vPath.size(); ++i)
        {
            path_type const& path = vPath[i];

            point_type const& p1 = path.low();
            point_type const& p2 = path.high();
            std::list<rectangle_pointer_type> vPattern1;
            std::list<rectangle_pointer_type> vPattern2;
            m_db->tPatternBbox.query(bgi::contains(p1), std::back_inserter(vPattern1));
            m_db->tPatternBbox.query(bgi::contains(p2), std::back_inserter(vPattern2));
#ifdef DEBUG
            if (vPattern1.size() > 1)
            {
                for (std::list<rectangle_pointer_type>::const_iterator it = vPattern1.begin(), ite = vPattern1.end(); it != ite; ++it)
                    mplPrint(kDEBUG, "multiple patterns found for layer %d, (%d, %d, %d, %d)\n", 
                            (*it)->layer(), gtl::xl(**it), gtl::yl(**it), gtl::xh(**it), gtl::yh(**it));
            }
            if (vPattern2.size() > 1)
            {
                for (std::list<rectangle_pointer_type>::const_iterator it = vPattern2.begin(), ite = vPattern2.end(); it != ite; ++it)
                    mplPrint(kDEBUG, "multiple patterns found for layer %d, (%d, %d, %d, %d)\n", 
                            (*it)->layer(), gtl::xl(**it), gtl::yl(**it), gtl::xh(**it), gtl::yh(**it));
            }
#endif
            mplAssert(vPattern1.size() == 1);
            mplAssert(vPattern2.size() == 1);
            rectangle_pointer_type const& pPattern1 = vPattern1.front();
            rectangle_pointer_type const& pPattern2 = vPattern2.front();
            coordinate_difference distance = gtl::length(path);

#ifdef _OPENMP
#pragma omp critical(dataupdate)
#endif
            {
                // if the input gds file contains duplicate edges 
                // we will have duplicate edges here 
                m_mAdjVertex[pPattern1->pattern_id()].push_back(pPattern2->pattern_id());
                m_mAdjVertex[pPattern2->pattern_id()].push_back(pPattern1->pattern_id());
                // estimate coloring distance 
                m_db->coloring_distance = std::max(distance, m_db->coloring_distance);
            }
        }
    }
    // eliminate all duplicates in adjacency list 
    uint32_t edge_num = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif 
    for (uint32_t v = 0; v < vertex_num; ++v)
    {
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];
        set<uint32_t> sAdjVertex (vAdjVertex.begin(), vAdjVertex.end());
        vAdjVertex.assign(sAdjVertex.begin(), sAdjVertex.end()); 
        vAdjVertex.swap(vAdjVertex); // shrink to fit 
        // parallel with omp reduction
        edge_num += vAdjVertex.size();
    }
    m_db->parms.coloring_distance_nm = m_db->coloring_distance*(m_db->unit*1e+9);
    mplPrint(kINFO, "Estimated coloring distance from conflict edges = %lld (%g nm)\n", m_db->coloring_distance, m_db->coloring_distance_nm());

    return edge_num;
}

void SimpleMPL::connected_component()
{
	mplPrint(kINFO, "Computing connected components...\n");
	uint32_t comp_id = 0;
	uint32_t order_id = 0; // position in an ordered pattern array with grouped components 

	uint32_t vertex_num = m_vVertexOrder.size();
	// here we employ the assumption that the graph data structure always has vertices labeled with 0~vertex_num-1
	// m_vVertexOrder only saves an order of it 
    for (uint32_t v = 0; v != vertex_num; ++v)
    {
        if (m_vCompId[v] == std::numeric_limits<uint32_t>::max()) // not visited
        {
            depth_first_search(v, comp_id, order_id);
            comp_id += 1;
        }
    }
	// record maximum number of connected components
    m_comp_cnt = comp_id;

#ifdef DEBUG 
	// check visited 
	for (uint32_t v = 0; v != vertex_num; ++v)
		mplAssert(m_vCompId[v] != std::numeric_limits<uint32_t>::max()); 
#endif

    // reorder with order_id 
    // maybe there is a way to save memory 
    std::vector<uint32_t> vTmpOrder (m_vVertexOrder.size());
    for (uint32_t v = 0; v != vertex_num; ++v)
        vTmpOrder[m_vVertexOrder[v]] = v;
    std::swap(m_vVertexOrder, vTmpOrder);

#ifdef DEBUG
	// check ordered 
	for (uint32_t v = 1; v != vertex_num; ++v)
		mplAssert(m_vCompId[m_vVertexOrder[v-1]] <= m_vCompId[m_vVertexOrder[v]]);
#endif
}

void SimpleMPL::depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& order_id)
{
    std::stack<uint32_t> vStack; 
	vStack.push(source);

	while (!vStack.empty())
	{
		uint32_t current = vStack.top();
		vStack.pop();
		if (m_vCompId[current] == std::numeric_limits<uint32_t>::max()) // not visited 
		{
			m_vCompId[current] = comp_id; // set visited 
			m_vVertexOrder[current] = order_id++; // update position 

			for (std::vector<uint32_t>::const_iterator it = m_mAdjVertex[current].begin(); it != m_mAdjVertex[current].end(); ++it)
			{
				uint32_t const& child = *it;
				mplAssert(current != child);
				vStack.push(child);
			}
		}
	}
}

///// helper functions for SimpleMPL::solve_component
///// since we may want to call solve_graph_coloring() multiple times 
///// it is better to wrap it as an independent function 

/// create coloring solver pointer according to algorithm type
lac::Coloring<SimpleMPL::graph_type>* SimpleMPL::create_coloring_solver(SimpleMPL::graph_type const& sg) const
{
    typedef lac::Coloring<graph_type> coloring_solver_type;
    coloring_solver_type* pcs = NULL;
    switch (m_db->algo().get())
    {
#if GUROBI == 1
        case AlgorithmTypeEnum::ILP_GURBOI:
            pcs = new lac::ILPColoring<graph_type> (sg); break;
        case AlgorithmTypeEnum::LP_GUROBI:
            pcs = new lac::LPColoring<graph_type> (sg); break;
        case AlgorithmTypeEnum::MIS_GUROBI:
            pcs = new lac::MISColoring<graph_type> (sg); break;
#endif
#if LEMONCBC == 1
        case AlgorithmTypeEnum::ILP_CBC:
            pcs = new lac::ILPColoringLemonCbc<graph_type> (sg); break;
#endif
#if CSDP == 1
        case AlgorithmTypeEnum::SDP_CSDP:
            pcs = new lac::SDPColoringCsdp<graph_type> (sg); break;
#endif
        case AlgorithmTypeEnum::BACKTRACK:
            pcs = new lac::BacktrackColoring<graph_type> (sg);
            break;
        default: mplAssertMsg(0, "unknown algorithm type");
    }
    pcs->stitch_weight(0.1);
    pcs->color_num(m_db->color_num());
    pcs->threads(1); // we use parallel at higher level 

    return pcs;
}

void SimpleMPL::vdd_biconnected_component()
{

}

/// given a graph, solve coloring 
/// contain nested call for itself 
uint32_t SimpleMPL::solve_graph_coloring(uint32_t comp_id, SimpleMPL::graph_type const& dg, 
        std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        uint32_t simplify_strategy, std::vector<int8_t>& vColor, std::set<vertex_descriptor> vdd_set)
{
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (dg, m_db->color_num());
	gs.set_isVDDGND(vdd_set);
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
#ifdef DEBUG
    if (comp_id == 130)
        mplPrint(kDEBUG, "stop\n"); 
#endif
    // set max merge level, actually it only works when MERGE_SUBK4 is on 
    if (m_db->color_num() == 3)
        gs.max_merge_level(3);
    else if (m_db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works better 
        gs.max_merge_level(2);
    gs.simplify(simplify_strategy);
	// collect simplified information 
    std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

    // for debug, it does not affect normal run 
	if (comp_id == m_db->dbg_comp_id() && simplify_strategy != graph_simplification_type::MERGE_SUBK4)
		gs.write_simplified_graph_dot("graph_simpl");

	// in order to recover color from articulation points 
	// we have to record all components and mappings 
	// but graph is not necessary 
	std::vector<std::vector<int8_t> > mSubColor (gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig (gs.num_component());
	double acc_obj_value = 0;
	std::cout << "number of components : " << gs.num_component() << std::endl;
	
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		std::cout << std::endl;
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		
		vSubColor.assign(num_vertices(sg), -1);

		for (std::map<uint32_t, uint32_t>::iterator it = dgGlobal2Local.begin(); it != dgGlobal2Local.end(); it++)
		{
			if (std::find(vSimpl2Orig.begin(), vSimpl2Orig.end(), it->second) != vSimpl2Orig.end())
			{
				dgCompId[it->first] = sub_comp_id + 1;
			}
		}

		// solve coloring 
		typedef lac::Coloring<graph_type> coloring_solver_type;
		coloring_solver_type* pcs = create_coloring_solver(sg);

		// set precolored vertices 
        boost::graph_traits<graph_type>::vertex_iterator vi, vie;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = vColor[vSimpl2Orig[v]];
			if (color >= 0 && color < m_db->color_num())
            {
				pcs->precolor(v, color);
                vSubColor[v] = color; // necessary for 2nd trial
            }
		}
        // 1st trial 
		
		double obj_value1 = (*pcs)(); // solve coloring 
	/*	
#ifdef DEBUG
        mplPrint(kDEBUG, "comp_id = %u, %lu vertices, obj_value1 = %g\n", comp_id, num_vertices(sg), obj_value1); 
#endif
        // 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
        double obj_value2 = std::numeric_limits<double>::max();
	
#ifndef DEBUG_NONINTEGERS
		std::set<vertex_descriptor> s_vdd_set;
		for (uint32_t i = 0; i < vSimpl2Orig.size(); i++)
		{
			if (vdd_set.find(vSimpl2Orig[i]) != vdd_set.end())
			{
				s_vdd_set.insert(i);
			}
		}
        // very restrict condition to determin whether perform MERGE_SUBK4 or not 
        if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::LP_GUROBI || m_db->algo() == AlgorithmTypeEnum::SDP_CSDP)
                && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
            obj_value2 = solve_graph_coloring(comp_id, sg, itBgn, pattern_cnt, graph_simplification_type::MERGE_SUBK4, vSubColor, s_vdd_set); // call again 
#endif

        // choose smaller objective value 
        if (obj_value1 < obj_value2)
        {
            acc_obj_value += obj_value1;

            // collect coloring results from simplified graph 
            for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
            {
                vertex_descriptor v = *vi;
                int8_t color = pcs->color(v);
				// if (color < 0) color = m_db->color_num();
                mplAssert(color >= 0 && color < m_db->color_num());
                vSubColor[v] = color;
            }
        }
        else // no need to update vSubColor, as it is already updated by sub call 
            acc_obj_value += obj_value2;
		delete pcs;
		*/
	}

#ifdef _DGOUT
	return 1;
#endif
	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 
	gs.recover(vColor, mSubColor, mSimpl2Orig);

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
    RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();

    return acc_obj_value;
}

void SimpleMPL::construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt, 
        SimpleMPL::graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor,std::set<vertex_descriptor>& vdd_set, bool flag) const
{
	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		vColor[i] = m_db->vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
		if (isVDDGND[v])
			vdd_set.insert(i);
	}

	// edges 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
		///< add conflict edges
		for (std::vector<uint32_t>::const_iterator it = m_mAdjVertex[v].begin(); it != m_mAdjVertex[v].end(); ++it)
		{
			uint32_t u = *it;
#ifdef DEBUG
			mplAssert(mGlobal2Local.count(u));
#endif
			uint32_t j = mGlobal2Local[u];
			if (i < j) // avoid duplicate 
			{
				std::pair<edge_descriptor, bool> e = edge(i, j, dg);
				if (!e.second) // make sure no duplicate 
				{
					e = add_edge(i, j, dg);
					mplAssert(e.second);
                    boost::put(boost::edge_weight, dg, e.first, 1);
				}
			}
		} // for, end of conflict edges
		
		///< add stitch edges
		if (flag == true)
		{
			for (std::vector<uint32_t>::const_iterator it = StitchRelation[v].begin(); it != StitchRelation[v].end(); ++it)
			{
				uint32_t s = *it;

				mplAssert(mGlobal2Local.count(s));
				uint32_t j = mGlobal2Local[s];
				if (i < j)
				{
					std::pair<edge_descriptor, bool> e = edge(i, j, dg);
					if (!e.second)
					{
						e = add_edge(i, j, dg);
						mplAssert(e.second);
						// for stitch, edge_weight is negative.
						boost::put(boost::edge_weight, dg, e.first, -1);
					}
				}
			}
		}
	}
}

uint32_t SimpleMPL::solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	if (itBgn == itEnd) return 0;
#ifdef DEBUG
    // check order 
	for (std::vector<uint32_t>::const_iterator it = itBgn+1; it != itEnd; ++it)
	{
		uint32_t v1 = *(it-1), v2 = *it;
		mplAssert(m_vCompId[v1] == m_vCompId[v2]);
	}
#endif

    uint32_t acc_obj_value = std::numeric_limits<uint32_t>::max();
    // if current pattern does not contain uncolored patterns, directly calculate conflicts 
    if (check_uncolored(itBgn, itEnd)) 
        acc_obj_value = coloring_component(itBgn, itEnd, comp_id);
#ifdef _DGOUT
	return 1;
#endif
	// update global color density map 
	// if parallelization is enabled, there will be uncertainty in the density map 
	// because the density is being updated while being read 
	for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
        int8_t color = m_db->vPatternBbox[*it]->color();
		uint32_t& color_density = m_vColorDensity[color];
#ifdef _OPENMP
#pragma omp atomic
#endif
		color_density += 1;
	}

	uint32_t component_conflict_num = conflict_num(itBgn, itEnd);
    // only valid under no stitch 
    // if (acc_obj_value != std::numeric_limits<uint32_t>::max())
    //    mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, (uint32_t)(itEnd-itBgn), component_conflict_num);

	return component_conflict_num;
}

uint32_t SimpleMPL::coloring_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	// construct a graph for current component 
	uint32_t pattern_cnt = itEnd-itBgn;

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...\n", comp_id, pattern_cnt);

	// decomposition graph 
    // must allocate memory here 
	graph_type dg (pattern_cnt);
	std::vector<int8_t> vColor (pattern_cnt, -1); // coloring results 
	map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

	bool flag;
    // construct decomposition graph for component 
	if (m_db->use_stitch())
		flag = true;
	else
		flag = false;

	std::set<vertex_descriptor> vdd_set;
	construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, flag);

	std::map<uint32_t, uint32_t>().swap(dgGlobal2Local);

	dgGlobal2Local.insert(mGlobal2Local.begin(), mGlobal2Local.end());

    // for debug, it does not affect normal run 
	if (comp_id == m_db->dbg_comp_id())
        write_graph(dg, "graph_init");

	// graph simplification 
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
    uint32_t simplify_strategy = graph_simplification_type::NONE;
	// keep the order of simplification 
	if (m_db->simplify_level() > 1)
         simplify_strategy |= graph_simplification_type::HIDE_SMALL_DEGREE;
	if (m_db->simplify_level() > 2)
        simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

    // solve graph coloring 
    uint32_t acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);

#ifdef _DGOUT
	return 1;
#endif

#ifdef DEBUG
	/*
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		mplAssert(vColor[i] >= 0 && vColor[i] < m_db->color_num());
	*/
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		int8_t color = vColor[i];
		if (color < 0) color = m_db->color_num();
        m_db->set_color(v, color);
	}

    return acc_obj_value;
}

uint32_t SimpleMPL::conflict_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd) const
{
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	uint32_t cnt = 0;
    for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
		uint32_t v = *it;
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				int8_t color2 = vPatternBbox[u]->color();
				if (color2 >= 0 && color2 < m_db->color_num())
				{
					if (color1 == color2) ++cnt;
				}
				else ++cnt; // uncolored vertex is counted as conflict 
			}
		}
		else ++cnt; // uncolored vertex is counted as conflict 
	}
	// conflicts will be counted twice 
	return (cnt>>1);
}

uint32_t SimpleMPL::conflict_num() const
{
	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 <= m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
                if (v < u) // avoid duplicate 
                {
                    int8_t color2 = vPatternBbox[u]->color();
                    if (color2 >= 0 && color2 <= m_db->color_num())
                    {
                        if (color1 == color2) 
                            m_vConflict.push_back(std::make_pair(v, u));
                    }
                    else // uncolored vertex is counted as conflict 
                        mplAssertMsg(0, "uncolored vertex %u = %d", u, color2);
                }
			}
		}
        else // uncolored vertex is counted as conflict 
            mplAssertMsg(0, "uncolored vertex %u = %d", v, color1);
	}
	return m_vConflict.size();
}

bool SimpleMPL::check_uncolored(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd) const
{
    for (; itBgn != itEnd; ++itBgn)
        if (m_db->vPatternBbox[*itBgn]->color() < 0)
            return true;
    return false;
}

void SimpleMPL::write_graph(SimpleMPL::graph_type& g, std::string const& filename) const
{
    // in order to make the .gv file readable by boost graphviz reader 
    // I dump it with boost graphviz writer 
    boost::dynamic_properties dp;
    dp.property("id", boost::get(boost::vertex_index, g));
    dp.property("node_id", boost::get(boost::vertex_index, g));
    dp.property("label", boost::get(boost::vertex_index, g));
    // somehow edge properties need mutable graph_type& 
    dp.property("weight", boost::get(boost::edge_weight, g));
    dp.property("label", boost::get(boost::edge_weight, g));
    std::ofstream out ((filename+".gv").c_str());
    boost::write_graphviz_dp(out, g, dp, string("id"));
    out.close();
    la::graphviz2pdf(filename);
}

void SimpleMPL::print_welcome() const
{
  mplPrint(kNONE, "\n\n");
  mplPrint(kNONE, "=======================================================================\n");
  mplPrint(kNONE, "                      OpenMPL - Version 1.1                          \n");
  mplPrint(kNONE, "                                by                                   \n");  
  mplPrint(kNONE, "                Yibo Lin, Bei Yu, Qi Sun and  David Z. Pan           \n");
  mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
  mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
  mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
  mplPrint(kNONE, "            Contact Authors:  {yibolin,dpan}@cerc.utexas.edu         \n");
  mplPrint(kNONE, "                              {byu, qsun}@cse.cuhk.edu.hk            \n");
  mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE

