/*************************************************************************
> File Name: SimpleMPL.cpp
> Author: Yibo Lin, Qi Sun
> Mail: yibolin@utexas.edu, qsun@cse.cuhk.edu.hk
> Created Time: Wed May 20 22:38:50 2015
************************************************************************/

/*
solve_component : 
	mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);
LayoutDBPolygon:
	set_color()
*/
#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"
#include <sstream>
#include <time.h>
#include <stack>
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

void SimpleMPL::run(int32_t argc, char** argv)
{
	this->reset(false);
	this->read_cmd(argc, argv);
	this->read_gds();
	this->solve();
#ifdef COMPONENTS
	return;
#endif
	if(m_db->gen_stitch())
		return;
	this->report();
	this->write_gds();
	return;
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

void SimpleMPL::read_cmd(int32_t argc, char** argv)
{
	// a little bit tricky here 
	// in order to support run-time switch of layoutdb_type
	// we need to construct a dummy ControlParameter for CmdParser
	// then construct actual layoutdb_type according to the option of CmdParser

	// construct a dummy ControlParameter
	ControlParameter tmpParms;
	// read command 
	CmdParser cmd(tmpParms);
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
	boost::timer::auto_cpu_timer timer(buf);
	mplPrint(kINFO, "Reading input file %s\n", m_db->input_gds().c_str());
	// read input gds file 
	GdsReader reader(*m_db);
	mplAssertMsg(reader(m_db->input_gds()), "failed to read %s", m_db->input_gds().c_str());
	// must call initialize after reading 
	m_db->initialize_data();
	// report data 
	m_db->report_data();
}

void SimpleMPL::gen_proj_target()
{
	uint32_t count = 0;
	proj_target.resize(m_db->vPatternBbox.size(), false);
	std::vector<uint32_t> vBookmark(m_comp_cnt);
	// std::cout << "==== After Projection vBookmark ====" << std::endl;
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}
	/*
	std::cout << "bookmark : " << std::endl;
	for(uint32_t i = 0; i < m_comp_cnt; i++)
	{
		std::cout << "comp " << i << " starts at " << vBookmark[i] << std::endl;
	}
	*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		if (check_uncolored(itBgn, itEnd))
		{
			uint32_t pattern_cnt = itEnd - itBgn;
			graph_type dg(pattern_cnt);
			std::vector<int8_t> vColor(pattern_cnt, -1); // coloring results 
			map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

			construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, false);

			// graph simplification 
			typedef lac::GraphSimplification<graph_type> graph_simplification_type;
			uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
			// simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

			graph_simplification_type gs(dg, m_db->color_num());
			gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 

			gs.simplify(simplify_strategy);
			std::vector<bool> vSim2OriTotal(mGlobal2Local.size(), false);

			for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
			{
				graph_type sg;
				std::vector<vertex_descriptor> vSimpl2Orig;
				gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
				for(std::vector<vertex_descriptor>::iterator it = vSimpl2Orig.begin(), ite = vSimpl2Orig.end(); it != ite; ++it)
					vSim2OriTotal[*it] = true;
			}

			std::map<uint32_t, uint32_t>::iterator it, ite;
			for(it = mGlobal2Local.begin(), ite = mGlobal2Local.end(); it != ite; ++it)
			{
				if(vSim2OriTotal[it->second])
				{
					proj_target[it->first] = true;
					count ++;
					rectangle_pointer_type const& pPattern = m_db->vPatternBbox[it->first];
					//std::cout <<"proj_target " << it->first << " : " << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl;
				}
			}
			
		}
	}
	std::cout << "number of projection : " << count << std::endl;
}


void SimpleMPL::write_gds()
{
	char buf[256];
	mplSPrint(kINFO, buf, "writing output file takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);
	if (m_db->output_gds().empty())
	{
		mplPrint(kWARN, "Output file not specified, no file generated\n");
		return;
	}
	// write output gds file 
	GdsWriter writer;
	if(m_db->stitch())
	{
		std::vector<std::vector<uint32_t> > Final_Stitches;
		uint32_t final_stitch_num = stitch_num(Final_Stitches);
		mplPrint(kINFO, "Total %u final stitches.\n", final_stitch_num);
		mplPrint(kINFO, "Write output gds file: %s\n", (m_db->output_gds() + "_" + limbo::to_string(m_db->color_num()) + "_coloring_with_stitch.gds").c_str());
		writer(m_db->output_gds() + "_" + limbo::to_string(m_db->color_num()) + "_coloring_with_stitch.gds", *m_db, m_vConflict, Final_Stitches, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
	else
	{	
		mplPrint(kINFO, "Write output gds file: %s\n", (m_db->output_gds() + "_" + limbo::to_string(m_db->color_num()) + "_coloring_no_stitch_result.gds").c_str());
		writer(m_db->output_gds() + "_" + limbo::to_string(m_db->color_num()) + "_coloring_result.gds", *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
}

void SimpleMPL::solve()
{
	// skip if no uncolored layer
	if (m_db->parms.sUncolorLayer.empty())
		return;
	
	if (m_db->vPatternBbox.empty())
	{
		mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}

	clock_t cons_start = clock();
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

	clock_t cons_end = clock();
	mplPrint(kINFO, "construct_graph takes  %f.\n", (double)(cons_end - cons_start)/CLOCKS_PER_SEC);

	if(m_db->stitch() || m_db->gen_stitch())
	{
		gen_proj_target();
		runProjection();
	}	
	
	if(m_db->gen_stitch())
		return;

	char buf[256];
	mplSPrint(kINFO, buf, "stitch coloring takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);

	// create bookmark to index the starting position of each component

	std::vector<uint32_t> vBookmark(m_comp_cnt);
	// std::cout << "==== After Projection vBookmark ====" << std::endl;
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

	mplPrint(kINFO, "Solving %u independent components...\n", m_comp_cnt);
	// thread number controled by user option 
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
#ifdef DEBUG
		//if (comp_id != 130)
		//    continue; 
#endif
		// construct a component 
		clock_t comp_start = clock();
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		// solve component 
		// pass iterators to save memory 
		this->solve_component(itBgn, itEnd, comp_id);
		clock_t comp_end = clock();
		// mplPrint(kINFO, "Component %d takes %f seconds.\n", comp_id, (double)(comp_end - comp_start)/CLOCKS_PER_SEC);
	}
#ifdef COMPONENTS
	return ;
#endif

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

void SimpleMPL::report() const
{
	mplPrint(kINFO, "Total conflict number = %u\n", conflict_num());
	for (int32_t i = 0, ie = m_db->color_num(); i != ie; ++i)
		mplPrint(kINFO, "Color %d density = %u\n", i, m_vColorDensity[i]);
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
	edge_num = edge_num >> 1;

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
		rectangle_type rect(*pPattern);
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
				{
					/*
					std::cout << "distance : " << distance << " " <<"coloring_distance : " << m_db->coloring_distance   << std::endl;
					std::cout << gtl::xl(*pAdjPattern) << " " << gtl::yl(*pAdjPattern) << " " << gtl::xh(*pAdjPattern) << " " << gtl::yh(*pAdjPattern) << std::endl;
					std::cout << gtl::xl(*pPattern) << " " << gtl::yl(*pPattern) << " " << gtl::xh(*pPattern) << " " << gtl::yh(*pPattern) << std::endl;
					*/
					vAdjVertex.push_back(pAdjPattern->pattern_id());
				}
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
		set<uint32_t> sAdjVertex(vAdjVertex.begin(), vAdjVertex.end());
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
	std::vector<uint32_t> vTmpOrder(m_vVertexOrder.size());
	for (uint32_t v = 0; v != vertex_num; ++v)
		vTmpOrder[m_vVertexOrder[v]] = v;
	std::swap(m_vVertexOrder, vTmpOrder);

#ifdef DEBUG
	// check ordered 
	for (uint32_t v = 1; v != vertex_num; ++v)
		mplAssert(m_vCompId[m_vVertexOrder[v - 1]] <= m_vCompId[m_vVertexOrder[v]]);
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
		pcs = new lac::ILPColoring<graph_type>(sg); break;
	case AlgorithmTypeEnum::LP_GUROBI:
		pcs = new lac::LPColoring<graph_type>(sg); break;
	case AlgorithmTypeEnum::MIS_GUROBI:
		pcs = new lac::MISColoring<graph_type>(sg); break;
#endif
#if LEMONCBC == 1
	case AlgorithmTypeEnum::ILP_CBC:
		pcs = new lac::ILPColoringLemonCbc<graph_type>(sg); break;
#endif
#if CSDP == 1
	case AlgorithmTypeEnum::SDP_CSDP:
		pcs = new lac::SDPColoringCsdp<graph_type>(sg); break;
#endif
	case AlgorithmTypeEnum::BACKTRACK:
		pcs = new lac::BacktrackColoring<graph_type>(sg);
		break;
	default: mplAssertMsg(0, "unknown algorithm type");
	}
	pcs->stitch_weight(0.1);
	pcs->color_num(m_db->color_num());
	pcs->threads(1); // we use parallel at higher level 

	// std::cout << "\ncreate coloring solver down.\n";
	return pcs;
}

/// given a graph, solve coloring 
/// contain nested call for itself 
uint32_t SimpleMPL::solve_graph_coloring(uint32_t comp_id, SimpleMPL::graph_type const& dg,
	std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt,
	uint32_t simplify_strategy, std::vector<int8_t>& vColor) const
{
	// mplPrint(kINFO, "solve graph coloring.");
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs(dg, m_db->color_num());
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
	std::vector<std::vector<int8_t> > mSubColor(gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig(gs.num_component());
	double acc_obj_value = 0;
	// std::cout << "solve_graph_coloring !" << std::endl;

	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		clock_t sub_comp_start = clock();
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);

		vSubColor.assign(num_vertices(sg), -1);
		// std::cout << "solve subcomponent " << sub_comp_id << " has " << vSubColor.size() << " nodes." << std::endl;
/*
		std::cout << "\n============================ subcomponent : " << sub_comp_id << std::endl;
		for(std::vector<vertex_descriptor>::iterator it = vSimpl2Orig.begin(); it != vSimpl2Orig.end(); ++it)
		{
			rectangle_pointer_type const& pPattern = m_db->vPatternBbox[(*it) + (*itBgn)];
			std::cout <<"simplification " << *it << " : " << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl;
		}
		std::cout << "==============================" << std::endl;

		if(vSubColor.size() >= 80)
		{
			std::string filename = m_db->output_gds() + "g_" + limbo::to_string(comp_id) + "_" + limbo::to_string(sub_comp_id);
			write_graph(sg, filename);
		}
	
		for (uint32_t i = 0; i != pattern_cnt; ++i)
		{
			uint32_t const& v = *(itBgn + i);
			m_db->vPatternBbox[v]->color(5);
			m_db->set_color(v, 5);
		}
		std::string intermediate_name = m_db->output_gds() + limbo::to_string(comp_id) + "_" + limbo::to_string(sub_comp_id) + ".gds";
		GdsWriter writer;
		std::cout << "here should output gds file. It has " << pattern_cnt << " patterns." << std::endl;
		mplPrint(kINFO, "Write output component gds file: %s\n", intermediate_name.c_str());
		writer.write_intermediate(intermediate_name, m_db->polyrect_patterns(), 100, m_db->strname, m_db->unit*1e+6);
*/
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
/*
#ifdef QDEBUG
		int abutting_count = 0;
		for(uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
		{
			m_db->vPatternBbox[i]->color(0);
			m_db->set_color(i, 0);
		}
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			abutting_count ++ ;
			m_db->vPatternBbox[vSimpl2Orig[v]]->color(abutting_count % 5 + 1);
			m_db->set_color(vSimpl2Orig[v], abutting_count % 5 + 1);
		}
		std::string intermediate_name = m_db->output_gds() + "_comp_" + limbo::to_string(comp_id) + "_sub_" + limbo::to_string(sub_comp_id) + "_abutting.gds";
		GdsWriter writer;
		// mplPrint(kINFO, "Write output component gds file: %s\n", intermediate_name.c_str());
		writer.write_intermediate(intermediate_name, m_db->polyrect_patterns(), 100, m_db->strname, m_db->unit*1e+6);
#endif
*/
		// 1st trial 
		// std::cout << "now start solving : " << std::endl;
		double obj_value1 = (*pcs)(); // solve coloring 
		// std::cout << "solving done.\n" << std::endl;
#ifdef DEBUG
		mplPrint(kDEBUG, "comp_id = %u, %lu vertices, obj_value1 = %g\n", comp_id, num_vertices(sg), obj_value1);
#endif
		// 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
		double obj_value2 = std::numeric_limits<double>::max();
/*
#ifndef DEBUG_NONINTEGERS
		// very restrict condition to determine whether perform MERGE_SUBK4 or not 
		if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::LP_GUROBI || m_db->algo() == AlgorithmTypeEnum::SDP_CSDP)
			&& (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
			obj_value2 = solve_graph_coloring(comp_id, sg, itBgn, pattern_cnt, graph_simplification_type::MERGE_SUBK4, vSubColor); // call again 
#endif
*/		
#ifdef COMPONENTS
		write_graph(sg, "c_" + limbo::to_string(comp_id) + "_" + limbo::to_string(sub_comp_id));
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
				mplAssert(color >= 0 && color < m_db->color_num());
				vSubColor[v] = color;
			}
		}
		else // no need to update vSubColor, as it is already updated by sub call 
			acc_obj_value += obj_value2;

#ifdef QDEBUG
		for(uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
		{
			m_db->vPatternBbox[i]->color(0);
			m_db->set_color(i, 0);
		}
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = pcs->color(v);
			m_db->vPatternBbox[vSimpl2Orig[v]]->color(color + 1);
			m_db->set_color(vSimpl2Orig[v], color + 1);
		}
		std::string intermediate_name = m_db->output_gds() + "_comp_" + limbo::to_string(comp_id) + "_sub_" + limbo::to_string(sub_comp_id) + "_coloring.gds";
		GdsWriter writer;
		// mplPrint(kINFO, "Write output component gds file: %s\n", intermediate_name.c_str());
		writer.write_intermediate(intermediate_name, m_db->polyrect_patterns(), 100, m_db->strname, m_db->unit*1e+6);
#endif

		clock_t sub_comp_end = clock();
		mplPrint(kINFO, "Comp_%d_subcomp_%d has %d nodes, takes %fs\n\n\n", comp_id, sub_comp_id, vSubColor.size(), (double)(sub_comp_end - sub_comp_start)/CLOCKS_PER_SEC);

		delete pcs;
	}

#ifdef COMPONENTS
	return 1;
#endif
	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 
	clock_t recover_start = clock();
	gs.recover(vColor, mSubColor, mSimpl2Orig);

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
	RecoverHiddenVertex(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();
	// RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();
	clock_t recover_end = clock();
	// mplPrint(kINFO, "\nComp_%d recovery takes %f seconds.\n", comp_id, (double)(recover_end - recover_end)/CLOCKS_PER_SEC);
	return acc_obj_value;
}

void SimpleMPL::construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt,
	SimpleMPL::graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor, bool flag) const
{
	
	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn + i);
		vColor[i] = m_db->vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
	}
	uint32_t stitch_count = 0;
	// edges 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn + i);
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
					/*
					rectangle_pointer_type const & tempA = m_db->vPatternBbox[v];
					rectangle_pointer_type const & tempB = m_db->vPatternBbox[u];
					std::cout << "add edge : " << std::endl;
					std::cout << "left " << v << " : " << gtl::xl(*tempA) << ", " << gtl::yl(*tempA) << ", " << gtl::xh(*tempA) << ", " << gtl::yh(*tempA) << std::endl;
					std::cout << "right " << u << " : " << gtl::xl(*tempB) << ", " << gtl::yl(*tempB) << ", " << gtl::xh(*tempB) << ", " << gtl::yh(*tempB) << std::endl;
					*/
					boost::put(boost::edge_weight, dg, e.first, 1);
				}
			}
		}
		if (flag && m_db->stitch())
		{
			for(std::vector<uint32_t>::const_iterator it = StitchRelation[v].begin(); it != StitchRelation[v].end(); ++it)
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
						stitch_count += 1;
/*
						rectangle_pointer_type tempA = m_db->vPatternBbox[v];
						std::cout << v << " : " << gtl::xl(*tempA) << ", " << gtl::yl(*tempA) << ", " << gtl::xh(*tempA) << ", " << gtl::yh(*tempA) << std::endl;
						std::cout << "inserts stitch with\n";
						rectangle_pointer_type tempB = m_db->vPatternBbox[s];
						std::cout << s << " : " << gtl::xl(*tempB) << ", " << gtl::yl(*tempB) << ", " << gtl::xh(*tempB) << ", " << gtl::yh(*tempB) << std::endl << std::endl;
*/
					}
				}
			}
		}
	}
	// mplPrint(kINFO, "%u stitches inserted.\n", stitch_count);
}

uint32_t SimpleMPL::solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{

	if (itBgn == itEnd) 
	{
		return 0; 
	}

#ifdef DEBUG
	// check order
	/*
	for (std::vector<uint32_t>::const_iterator it = itBgn + 1; it != itEnd; it++)
	{
		uint32_t v1 = *(it - 1), v2 = *it;
		mplAssert(m_vCompId[v1] == m_vCompId(v2);
	}
	*/
#endif
	/*
	std::cout << "In solve component : "  << comp_id << std::endl; 
	for(std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; it++)
	{
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[*it];
		std::cout << pPattern->pattern_id() << " : " << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl;
	}
	*/
	uint32_t acc_obj_value = std::numeric_limits<uint32_t>::max();
	// if current pattern does not contain uncolored patterns, directly calculate conflicts 
	if (check_uncolored(itBgn, itEnd))
		acc_obj_value = coloring_component(itBgn, itEnd, comp_id);
	else
		std::cout << "precolored." << std::endl;

#ifdef COMPONENTS
	return acc_obj_value;
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
/*
	// only valid under no stitch 
	if(!m_db->stitch())
	{
		if (acc_obj_value != std::numeric_limits<uint32_t>::max())
			mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);
	}
	else {
		uint32_t component_stitch_num = stitch_num(itBgn, itEnd);
		mplPrint(kINFO, "Component %u has %u patterns...%u final stitches\n\n", comp_id, (uint32_t)(itEnd - itBgn), component_stitch_num);
	}
*/
	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, (uint32_t)(itEnd - itBgn), component_conflict_num);		

	return component_conflict_num;
}

uint32_t SimpleMPL::coloring_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	// construct a graph for current component 
	uint32_t pattern_cnt = itEnd - itBgn;

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...\n", comp_id, pattern_cnt);


	// decomposition graph 
	// must allocate memory here 
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1); // coloring results 
	map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

	// mplPrint(kINFO, "In component %u\n", comp_id);
										   // construct decomposition graph for component
	construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, true);

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
	uint32_t acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor);

#ifdef COMPONENTS
	return acc_obj_value;
#endif

#ifdef DEBUG
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		mplAssert(vColor[i] >= 0 && vColor[i] < m_db->color_num());
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn + i);
		m_db->set_color(v, vColor[i]);
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
	return (cnt >> 1);
}

uint32_t SimpleMPL::conflict_num() const
{
	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				if (v < u) // avoid duplicate 
				{
					int8_t color2 = vPatternBbox[u]->color();
					if (color2 >= 0 && color2 < m_db->color_num())
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
	// std::ofstream out(("benchout/" + filename + ".gv").c_str());
	std::ofstream out((filename + ".gv").c_str());
	boost::write_graphviz_dp(out, g, dp, string("id"));
	out.close();
	la::graphviz2pdf(filename);
}

void SimpleMPL::runProjection()
{
	char buf[256];
	mplSPrint(kINFO, buf, "stitch generation takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);

	// initialization
	uint32_t vertex_num = m_db->vPatternBbox.size();
	std::vector<uint32_t> new_vCompId_vec;
	std::vector<rectangle_pointer_type> rect_vec = m_db->polyrect_patterns();
	std::vector<uint32_t> poly_rect_begin = m_db->polyrectBgnId();
	mplAssertMsg(poly_rect_begin.size() == vertex_num, "polyrectBgnId.size() must be equal to polyrect_patterns.size()");

	// generate poly_rect_end, which stores the end rectangle index when querying from parent polygon id
	std::vector<uint32_t> poly_rect_end;
	poly_rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		poly_rect_end[i] = poly_rect_begin[i + 1] - 1;
	poly_rect_end[vertex_num - 1] = rect_vec.size() - 1;

/*
#ifdef QDEBUG
	for(uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
	{
		rectangle_pointer_type pPattern = m_db->vPatternBbox[i];
		uint32_t start_idx = poly_rect_begin[pPattern->pattern_id()];
		uint32_t end_idx = poly_rect_end[pPattern->pattern_id()];
		int count = 0;
		for(uint32_t j = start_idx; j <= end_idx; j++)
		{
			rect_vec[j]->color(count);
			count ++;
		}
	}

	GdsWriter writer;
	writer.write_intermediate(m_db->output_gds() + "_abutting.gds", rect_vec, 1, m_db->strname, m_db->unit*1e+6);
#endif
*/
	std::vector<rectangle_pointer_type> new_rect_vec;		// store the newly-generated rectangles
	std::vector<uint32_t> rect_to_parent;					// map from rectangles to its parent polygon
	std::vector<std::vector<uint32_t> >().swap(ori2new);
	ori2new.resize(vertex_num);

	// store new vertex order
	std::vector<uint32_t> new_vertex_order;

	uint32_t new_polygon_id = -1;
	uint32_t new_rect_id = -1;

	clock_t projection_start = clock();
	double reconstruct_polygon_total = 0;
	double project_total = 0;

	for (uint32_t ver = 0; ver < vertex_num; ver++)
	{
		uint32_t v = m_vVertexOrder[ver];
		uint32_t comp_id = m_vCompId[v];
		
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[v];
		uint32_t pid = pPattern->pattern_id();
		// polygon v's neighbor polygons
		std::vector<uint32_t> & nei_Vec = m_mAdjVertex[pid];
		uint32_t start_idx = poly_rect_begin[pid];
		uint32_t end_idx = poly_rect_end[pid];

		if(proj_target[pid])
		{
			// std::cout << pid << " generate stitches : " << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl; 
			// use poss_nei_vec to obtain all the possible neighbor rectangles from neighbor polygons
			std::vector<rectangle_pointer_type> poss_nei_vec;
			//std::cout << "=====" << pid << " has neighbors : \n";
			for (std::vector<uint32_t>::iterator it = nei_Vec.begin(); it != nei_Vec.end(); it++)
			{
				uint32_t s_idx = poly_rect_begin[*it];
				uint32_t e_idx = poly_rect_end[*it];
				for (uint32_t a = s_idx; a <= e_idx; a++)
				{
					//std::cout << gtl::xl(*rect_vec[a]) << " " << gtl::yl(*rect_vec[a]) << " " << gtl::xh(*rect_vec[a]) << " " << gtl::yh(*rect_vec[a]) << std::endl;
					poss_nei_vec.push_back(rect_vec[a]);
				}
			}

			std::vector<std::pair<rectangle_pointer_type, uint32_t> > poly_split; 
			std::vector<uint32_t> new_polygon_id_list;
			// traverse all the rectangles in polygon v, to generate the intersections
			for (uint32_t j = start_idx; j <= end_idx; j++)
			{
				rectangle_type rect(*rect_vec[j]);
				// the generated split patterns
				std::vector<rectangle_pointer_type> split;
				clock_t each_pro = clock();
				projection(rect, split, poss_nei_vec);
				clock_t each_pend = clock();
				project_total += (double)(each_pend - each_pro)/CLOCKS_PER_SEC;

				for(std::vector<rectangle_pointer_type>::iterator it = split.begin(); it!=split.end(); it++)
					poly_split.push_back(std::make_pair(*it, rect.pattern_id()));
			}
			uint32_t pivot = new_polygon_id;

			std::vector<std::vector<uint32_t> > stitch_list;
			clock_t each_re = clock();
			reconstruct_polygon(new_polygon_id, new_polygon_id_list, poly_split, stitch_list);
			clock_t each_end = clock();
			reconstruct_polygon_total += (double)(each_end - each_re)/CLOCKS_PER_SEC;

			assert(new_polygon_id_list.size() == poly_split.size());

			StitchRelation.insert(StitchRelation.end(), stitch_list.begin(), stitch_list.end());

			for(uint32_t i = 0; i < poly_split.size(); i++)
			{
				poly_split[i].first->pattern_id(++new_rect_id);
				if(m_db->gen_stitch())
					poly_split[i].first->color(new_polygon_id_list[i]%7);
				new_rect_vec.push_back(poly_split[i].first);
				rect_to_parent.push_back(new_polygon_id_list[i]);

				if(pivot != new_polygon_id_list[i])
				{
					pivot = new_polygon_id_list[i];
					new2ori.push_back(pid);
					ori2new[pid].push_back(pivot);
					new_vertex_order.push_back(pivot);
					new_vCompId_vec.push_back(comp_id);
				}
			}
		}
		else
		{
			new_polygon_id += 1;
			for (uint32_t j = start_idx; j <= end_idx; j++)
			{
				rectangle_pointer_type rect = rect_vec[j];
				rect->pattern_id(++new_rect_id);
				if(m_db->gen_stitch())
					rect->color(new_polygon_id%7);
				new_rect_vec.push_back(rect);
				rect_to_parent.push_back(new_polygon_id);
			}
			new2ori.push_back(pid);
			new_vertex_order.push_back(new_polygon_id);
			new_vCompId_vec.push_back(comp_id);
			ori2new[pid].push_back(new_polygon_id);
			StitchRelation.push_back(std::vector<uint32_t>());
		}
		
	}
	clock_t projection_end = clock();
	mplPrint(kINFO, "projection takes : %f seconds\n", (double)(projection_end - projection_start)/CLOCKS_PER_SEC );
	mplPrint(kINFO, "reconstruct_polygon_total : %f seconds\n", reconstruct_polygon_total);
/*
	std::cout << "\n\n===========================\n";
	for(uint32_t i = 0; i < new_rect_vec.size(); i++)
	{
		rectangle_pointer_type pPattern = new_rect_vec[i];
		std::cout << "polygon " << rect_to_parent[i] << " rect " << pPattern->pattern_id() << " : ";
		std::cout << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl << std::endl;
	}
*/	

	// update information in m_db;
	m_db->refresh(new_rect_vec, rect_to_parent);

	// Now, update the adjacency list, we still need original adjacency list
	std::vector<std::set<uint32_t> > new_mAdjVertex;
	new_mAdjVertex.resize(m_db->vPatternBbox.size());
	uint32_t edge_num = 0;

	// generate new index informaiton
	std::vector<uint32_t> new_poly_rect_begin = m_db->polyrectBgnId();
	std::vector<uint32_t> new_poly_rect_end;
	// update vertex_num
	vertex_num = new_poly_rect_begin.size();
	assert(vertex_num == m_db->vPatternBbox.size());

	new_poly_rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		new_poly_rect_end[i] = new_poly_rect_begin[i + 1] - 1;
	new_poly_rect_end[vertex_num - 1] = new_rect_vec.size() - 1;

	clock_t new_relation_start = clock();
	// traverse all the original polygons to get the original neighbor list
	for (uint32_t i = 0, ie = m_mAdjVertex.size(); i < ie; i++)
	{
		// list all the possible polygon neighbors
		std::vector<uint32_t> poss_nei;
		bool split_flag = false;
		for (uint32_t j = 0, je = m_mAdjVertex[i].size(); j < je; j++)
		{
			if(!proj_target[m_mAdjVertex[i][j]] && !proj_target[i])
				new_mAdjVertex[ori2new[i].front()].insert(ori2new[m_mAdjVertex[i][j]].front());
			else
				split_flag = true;
			for (std::vector<uint32_t>::iterator it = ori2new[m_mAdjVertex[i][j]].begin(); it != ori2new[m_mAdjVertex[i][j]].end(); it++)
				poss_nei.push_back(*it);

		}
		if(split_flag)
		{
			// traverse all the newly-generated polygons in current polygon
			for (std::vector<uint32_t>::iterator it = ori2new[i].begin(); it != ori2new[i].end(); it++)
			{
				// std::cout << *it << " : " << gtl::xl(*m_db->vPatternBbox[*it]) << ", " << gtl::yl(*m_db->vPatternBbox[*it]) << ", " << gtl::xh(*m_db->vPatternBbox[*it]) << ", " << gtl::yh(*m_db->vPatternBbox[*it]) << std::endl;
				// std::cout << "\n\n=========== new polygon " << *it << " possible neighbors: " << std::endl;
				uint32_t start_idx = new_poly_rect_begin[*it];
				uint32_t end_idx = new_poly_rect_end[*it];
				// traverse all rectangles in current newly-generated polygon.
				for(uint32_t now_rect = start_idx; now_rect <= end_idx; now_rect++)
				{
					// traverse all possible newly-generated neighbor polygons
					for(std::vector<uint32_t>::iterator nei_poly = poss_nei.begin(); nei_poly != poss_nei.end(); nei_poly++)
					{
						uint32_t nei_start_idx = new_poly_rect_begin[*nei_poly];
						uint32_t nei_end_idx = new_poly_rect_end[*nei_poly];
						// traverse all rectangles in current newly-generated neighbor polygon
						for(uint32_t nei_rect = nei_start_idx; nei_rect <= nei_end_idx; nei_rect++)
						{
							uint32_t nei_rect_pid = new_rect_vec[nei_rect]->pattern_id();
							coordinate_difference distance = boost::geometry::distance(*new_rect_vec[now_rect], *new_rect_vec[nei_rect_pid]);
							if(distance < m_db->coloring_distance)
							{
								new_mAdjVertex[*it].insert(*nei_poly);
								// if the rectangles are close to each other, it means corresponding polygons are neighbors
								break;
							}
						}
					}
				}
			}
		}
	}
	clock_t new_relation_end = clock();
	mplPrint(kINFO, "Generate new relationships takes %f seconds.\n", (double)(new_relation_end - new_relation_start)/CLOCKS_PER_SEC );

	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex.resize(new_mAdjVertex.size());
	for(uint32_t i = 0; i < new_mAdjVertex.size(); i++)
	{
		for(std::set<uint32_t>::iterator it = new_mAdjVertex[i].begin(); it != new_mAdjVertex[i].end(); it++)
		{
			edge_num++;
			m_mAdjVertex[i].push_back(*it);
		}
	}
	
	if(m_db->gen_stitch())
	{
		std::vector<std::pair<uint32_t, uint32_t> > stitch_pair;
		uint32_t stitch_number = 0;
		for(uint32_t i = 0; i < StitchRelation.size(); i++)
		{
			stitch_number += StitchRelation[i].size();
			for(uint32_t j = 0; j < StitchRelation[i].size(); j++)
			{
				stitch_pair.push_back(std::make_pair(i, StitchRelation[i][j]));
			}
		}

		GdsWriter writer;
		mplPrint(kINFO, "Stitches number: %u\n", stitch_number / 2);
		mplPrint(kINFO, "AdjVertex size: %u\n", m_mAdjVertex.size());
		mplPrint(kINFO, "Write output gds file: %s\n", (m_db->output_gds() + "_gen_stitch.gds").c_str());

		writer(m_db->output_gds() + "_gen_stitch.gds", *m_db, stitch_pair, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
/*
	std::cout << "\n\n\n======== stitch relationships ========== \n\n";
	for(uint32_t i = 0; i < StitchRelation.size(); i++)
	{
		rectangle_pointer_type tempA = m_db->vPatternBbox[i];
		for(uint32_t j = 0; j < StitchRelation[i].size(); j++)
		{
			rectangle_pointer_type tempB = m_db->vPatternBbox[StitchRelation[i][j]];
			std::cout << i << " : " << gtl::xl(*tempA) << ", " << gtl::yl(*tempA) << ", " << gtl::xh(*tempA) << ", " << gtl::yh(*tempA);
			std::cout << "\nhas stitch with\n" ;
			std::cout << StitchRelation[i][j] << " : " << gtl::xl(*tempB) << ", " << gtl::yl(*tempB) << ", " << gtl::xh(*tempB) << ", " << gtl::yh(*tempB) << std::endl << std::endl;
		}
	}
#ifdef QDEBUG
	std::cout << "\n\n\n========= conflict relationships ==========\n\n";
	for(uint32_t i = 0; i < m_mAdjVertex.size(); i++)
	{
		std::cout << i << " conflicts with : " ;
		for(uint32_t j = 0; j < m_mAdjVertex[i].size(); j++)
			std::cout << m_mAdjVertex[i][j] << " ";
		std::cout << std::endl;
	}
#endif
*/
	std::vector<uint32_t>().swap(m_vCompId);
	m_vCompId.swap(new_vCompId_vec);

	std::vector<uint32_t>().swap(m_vVertexOrder);
	m_vVertexOrder.swap(new_vertex_order);
/*
	for(uint32_t i = 0; i < m_vVertexOrder.size(); ++i)
	{
		std::cout << "order " << i << " is " << m_vVertexOrder[i];
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[m_vVertexOrder[i]];
		std::cout << " polygon " << pPattern->pattern_id() << " : " << gtl::xl(*pPattern) << ", " << gtl::yl(*pPattern) << ", " << gtl::xh(*pPattern) << ", " << gtl::yh(*pPattern) << std::endl;
	}
*/
	mplPrint(kINFO, "After stitch insertion, %u vertices, %u edges\n", m_db->vPatternBbox.size(), edge_num);
	return;
}

bool SimpleMPL::whetherHorizontal(rectangle_type temp)
{
	double xl = gtl::xl(temp);
	double yl = gtl::yl(temp);
	double xh = gtl::xh(temp);
	double yh = gtl::yh(temp);
	return (xh - xl) > (yh - yl);
}

LayoutDB::rectangle_type SimpleMPL::interSectionRect(rectangle_type rect1, rectangle_type rect2)
{
	coordinate_type xl = std::max(gtl::xl(rect1), gtl::xl(rect2));
	coordinate_type yl = std::max(gtl::yl(rect1), gtl::yl(rect2));
	coordinate_type xh = std::min(gtl::xh(rect1), gtl::xh(rect2));
	coordinate_type yh = std::min(gtl::yh(rect1), gtl::yh(rect2));
	rectangle_pointer_type output = new rectangle_type(xl, yl, xh, yh);
	return *output;
}

void SimpleMPL::projection(rectangle_type & pRect, std::vector<rectangle_pointer_type>& split, std::vector<rectangle_pointer_type> nei_Vec)
{
	bool hor = whetherHorizontal(pRect);
	std::vector<coordinate_type> vstitches;
	coordinate_difference width = gtl::xh(pRect) - gtl::xl(pRect);
	if(width <= 5000)
	{
		std::set<coordinate_type> vset;
		std::vector<coordinate_type> vPossibleStitches;
		if (hor)
		{
			vset.insert(gtl::xl(pRect));
			vset.insert(gtl::xh(pRect));
		}
		else {
			vset.insert(gtl::yl(pRect));
			vset.insert(gtl::yh(pRect));
		}
		std::vector<rectangle_type> vInterSect;
		// generate intersections
		for (std::vector<rectangle_pointer_type>::iterator it = nei_Vec.begin(); it != nei_Vec.end(); it++)
		{
			coordinate_difference distance = boost::geometry::distance(pRect, *(*it));
			if(distance >= m_db->coloring_distance) continue;
			rectangle_type extendPattern(*(*it));
			gtl::bloat(extendPattern, gtl::HORIZONTAL, m_db->coloring_distance);
			gtl::bloat(extendPattern, gtl::VERTICAL, m_db->coloring_distance);

			rectangle_type temp = interSectionRect(pRect, extendPattern);
			temp.pattern_id(extendPattern.pattern_id());
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

		// uint32_t nei_num = nei_Vec.size();
		
		GenerateStitchPosition_Bei(pRect, vInterSect, vPossibleStitches, vstitches);
		//GenerateStitchPosition_Jian(pRect, vInterSect, vPossibleStitches, nei_num, vstitches);

		// check the stitch positions' legalities
		// if the position is very colse to the rectangle's boundary, it's illegal.
		coordinate_type lower_boundary;
		coordinate_type upper_boundary;
		if (hor)
		{
			lower_boundary = gtl::xl(pRect);
			upper_boundary = gtl::xh(pRect);
		}
		else
		{
			lower_boundary = gtl::yl(pRect);
			upper_boundary = gtl::yh(pRect);
		}
		coordinate_type threshold = 0;
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

	}
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
					vstitches[j], gtl::yl(pRect), vstitches[j+1], gtl::yh(pRect));
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

void SimpleMPL::GenerateStitchPosition_Jian(const rectangle_type pRect, std::vector<rectangle_type> vInterSect, std::vector<coordinate_type> & vPossibleStitches, 
	uint32_t nei_num, std::vector<coordinate_type> & vstitches)
{
	bool ishor = whetherHorizontal(pRect);
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, std::set<uint32_t> > > vStages;
	std::set<uint32_t> temp;
	for (uint32_t i = 1; i < vPossibleStitches.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(vPossibleStitches[i - 1], vPossibleStitches[i]), temp));
	// calculate the times every stage covered by all the intersections.
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		coordinate_type lower = vStages[i].first.first;
		coordinate_type upper = vStages[i].first.second;
		for (uint32_t j = 0; j < vInterSect.size(); j++)
		{
			if (ishor)
			{
				if (lower < gtl::xl(vInterSect[j])) continue;
				if (upper > gtl::xh(vInterSect[j])) continue;
			}
			else
			{
				if (lower < gtl::yl(vInterSect[j])) continue;
				if (upper > gtl::yh(vInterSect[j])) continue;
			}
			vStages[i].second.insert(vInterSect[j].pattern_id());
		}
	}

	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		std::set<uint32_t> twosideSet;
		for (uint32_t j = 0; j <= i; j++)
			twosideSet.insert(vStages[j].second.begin(), vStages[j].second.end());
		if (twosideSet.size()  >= nei_num)
			continue;
		std::set<uint32_t>().swap(twosideSet);
		for (uint32_t j = i; j < vStages.size(); j++)
			twosideSet.insert(vStages[i].second.begin(), vStages[i].second.end());
		if (twosideSet.size() >= nei_num)
			continue;

		//std::cout << "Insert one stitch" << std::endl;
		vstitches.push_back((vStages[i].first.first + vStages[i].first.second) / 2);
	}
	sort(vstitches.begin(), vstitches.end());

	return;
}

void SimpleMPL::GenerateStitchPosition_Bei(const rectangle_type pRect, std::vector<rectangle_type> vInterSect, std::vector<coordinate_type> & vPossibleStitches,
	std::vector<coordinate_type> & vstitches)
{
	bool ishor = whetherHorizontal(pRect);
	coordinate_type lower, upper;
	if (ishor) {
		lower = gtl::xl(pRect);
		upper = gtl::xh(pRect);
	}
	else {
		lower = gtl::yl(pRect);
		upper = gtl::yh(pRect);
	}
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t > > vStages;
	for (uint32_t i = 1; i < vPossibleStitches.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(vPossibleStitches[i - 1], vPossibleStitches[i]), 0));
	// calculate the times every stage covered by all the intersections.
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
			}
			else
			{
				if (left < gtl::yl(vInterSect[j])) continue;
				if (right > gtl::yh(vInterSect[j])) continue;
			}
			vStages[i].second++;
		}
	}
	
	if (vStages.size() <= 0) return;
	if (vStages[0].second != 0)
		vStages.insert(vStages.begin(), std::make_pair(std::make_pair(lower, lower), 0));
	if (vStages[vStages.size() - 1].second != 0)
		vStages.push_back(std::make_pair(std::make_pair(upper, upper), 0));

	std::vector<uint32_t> vZeroIds;
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		if (vStages[i].second > 0) continue;
		vZeroIds.push_back(i);
	}

	std::vector<coordinate_type>().swap(vstitches);
	for (uint32_t i = 0; i < vZeroIds.size() - 1; i++)
	{
		uint32_t pos1 = vZeroIds[i];
		uint32_t pos2 = vZeroIds[i + 1];

		// since ((lower, lower), 0) has been added into vStages, so pos1 must be 0.
		if (i == 0) mplAssertMsg(0 == pos1, "pos1 %d doesn't equal to 0", pos1);
		// remove the useless stitches
		else if (pos1 == 2)
		{
			bool find = false;
			if (vStages.size() < 5) find = true;
			else if (i != 1) find = true;
			else if (1 != vStages[1].second) find = true;
			else if (1 != vStages[3].second) find = true;
			else if (0 != vStages[4].second) find = true;
			if (find == false)
				continue;
			coordinate_type position = (vStages[2].first.first + vStages[2].first.second) / 2;
			vstitches.push_back(position);
		}
		else if (pos1 == vStages.size() - 3)
		{
		//	std::cout << "i = " << i << "\t" << " vZeroIds.size() = " << vZeroIds.size() << std::endl;
			mplAssert(i == vZeroIds.size() - 2);
			bool find = false;
			uint32_t zsize = vZeroIds.size();
			if (vStages.size() < 5) find = true;
			else if (i != zsize - 2) find = true;
			else if (1 != vStages[pos1 + 1].second) find = true;
			else if (1 != vStages[pos1 - 1].second) find = true;
			else if (0 != vStages[pos1 - 2].second) find = true;
			if (find == false) continue;
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}
		else
		{
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}

		// search lost stitch in vStages[pos1 --> pos2]
		double maxValue = 0.9;
		uint32_t posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			uint32_t mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
			uint32_t diff = std::abs(static_cast<int>(vStages[i + 1].second - vStages[i - 1].second));
			double value = (double)mind + (double)diff * 0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			uint32_t position = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
			vstitches.push_back(position);
		}
	}
	sort(vstitches.begin(), vstitches.end());
	return;
}

void SimpleMPL::reconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t> & new_polygon_id_list, std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list, std::vector<std::vector<uint32_t> >& stitch_list)
{
	uint32_t start = polygon_id + 1;
	/*
	char buf[256];
	mplSPrint(kINFO, buf, " takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);
	*/
	std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_temp;
	std::vector<uint32_t> new_polygon_id_temp;

	uint32_t rect_list_size = rect_list.size();
	std::vector<bool> visited(rect_list_size, false);

	std::vector<uint32_t>().swap(new_polygon_id_list);
	new_polygon_id_list.resize(rect_list_size, std::numeric_limits<uint32_t>::max());

	for(uint32_t i = 0; i< rect_list_size; i++)
	{
		if(new_polygon_id_list[i] == std::numeric_limits<uint32_t>::max())
		{
			polygon_id += 1;
			new_polygon_id_list[i]=polygon_id;
		}
		visited[i] = true;
		for(uint32_t j = i+1; j < rect_list_size; j++)
		{
			if(visited[j] == false && rect_list[i].second != rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				visited[j] = true;
				new_polygon_id_list[j] = new_polygon_id_list[i];
			}
		}
	}

	stitch_list.resize( polygon_id - start + 1);
	for(uint32_t i = 0; i < rect_list.size(); i++)
	{
		for(uint32_t j = i + 1; j < rect_list.size(); j++)
		{
			if(rect_list[i].second == rect_list[j].second)
			{
				stitch_list[new_polygon_id_list[i] - start].push_back(new_polygon_id_list[j]);
				stitch_list[new_polygon_id_list[j] - start].push_back(new_polygon_id_list[i]);
			}
		}
	}
	for(uint32_t i = start, ie = polygon_id; i <= ie; ++i)
	{
		for(uint32_t j = 0; j < new_polygon_id_list.size(); ++j)
		{
			if(new_polygon_id_list[j] == i)
			{
				new_polygon_id_temp.push_back(i);
				rect_temp.push_back(rect_list[j]);
			}
		}
	}
	rect_list.swap(rect_temp);
	new_polygon_id_list.swap(new_polygon_id_temp);
}

uint32_t SimpleMPL::stitch_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd) const
{
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	uint32_t cnt = 0;
	for(std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
		uint32_t v = *it;
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for(std::vector<uint32_t>::const_iterator itAdj = StitchRelation[v].begin(); itAdj != StitchRelation[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				int8_t color2 = vPatternBbox[u]->color();
				if (color2 >= 0 && color2 < m_db->color_num())
				{
					if (color1 != color2 ) ++cnt;
				}
				else
					++cnt;
			}
		}
		else ++cnt;
	}
	return (cnt >> 1);
}

uint32_t SimpleMPL::stitch_num(std::vector<std::vector<uint32_t> >& Final_Stitches) const
{
	// std::cout << "\n\n=========== final stitch ===============\n";
	uint32_t cnt = 0;
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	Final_Stitches.resize(vPatternBbox.size());
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = StitchRelation[v].begin(); itAdj != StitchRelation[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				if (v < u) // avoid duplicate 
				{
					int8_t color2 = vPatternBbox[u]->color();
					if (color2 >= 0 && color2 < m_db->color_num())
					{
						if (color1 != color2)
						{
							Final_Stitches[v].push_back(u);
							Final_Stitches[u].push_back(v);							
							++cnt;
						}
					}
					else // uncolored vertex is counted as conflict 
						mplAssertMsg(0, "uncolored vertex %u = %d", u, color2);
				}
			}
		}
		else // uncolored vertex is counted as conflict 
			mplAssertMsg(0, "uncolored vertex %u = %d", v, color1);
	}
	return cnt;
}


void SimpleMPL::print_welcome() const
{
	mplPrint(kNONE, "\n\n");
	mplPrint(kNONE, "=======================================================================\n");
	mplPrint(kNONE, "                        OpenMPL - Version 1.1                        \n");
	mplPrint(kNONE, "                                by                                   \n");
	mplPrint(kNONE, "                Yibo Lin, Bei Yu, Qi Sun and David Z. Pan            \n");
	mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
	mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
	mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
	mplPrint(kNONE, "            Contact Authors:  {yibolin, dpan}@cerc.utexas.edu        \n");
	mplPrint(kNONE, "                              {byu, qsun}@cse.cuhk.edu.hk            \n");
	mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE

