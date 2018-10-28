/*************************************************************************
> File Name: SimpleMPL.cpp
> Author: Yibo Lin, Qi Sun
> Mail: yibolin@utexas.edu, qsun@cse.cuhk.edu.hk
> Created Time: Wed May 20 22:38:50 2015
************************************************************************/

#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"

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
	if (m_db->gen_stitch())
	{
		this->runProjection();
		this->write_gds();
		return;
	}
	//else {
		this->solve();
		this->report();
		this->write_gds();
	//}
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
	mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
	writer(m_db->output_gds(), *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
}

void SimpleMPL::solve()
{
	// skip if no uncolored layer
	if (m_db->parms.sUncolorLayer.empty())
		return;
	
	char buf[256];
	mplSPrint(kINFO, buf, "coloring takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);
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
	std::vector<uint32_t> vBookmark(m_comp_cnt);
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
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
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

void SimpleMPL::report() const
{
	mplPrint(kINFO, "Conflict number = %u\n", conflict_num());
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

	return pcs;
}


/// given a graph, solve coloring 
/// contain nested call for itself 
uint32_t SimpleMPL::solve_graph_coloring(uint32_t comp_id, SimpleMPL::graph_type const& dg,
	std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt,
	uint32_t simplify_strategy, std::vector<int8_t>& vColor) const
{
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
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);

		vSubColor.assign(num_vertices(sg), -1);

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
#ifdef DEBUG
		mplPrint(kDEBUG, "comp_id = %u, %lu vertices, obj_value1 = %g\n", comp_id, num_vertices(sg), obj_value1);
#endif
		// 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
		double obj_value2 = std::numeric_limits<double>::max();
#ifndef DEBUG_NONINTEGERS
		// very restrict condition to determine whether perform MERGE_SUBK4 or not 
		if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::LP_GUROBI || m_db->algo() == AlgorithmTypeEnum::SDP_CSDP)
			&& (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
			obj_value2 = solve_graph_coloring(comp_id, sg, itBgn, pattern_cnt, graph_simplification_type::MERGE_SUBK4, vSubColor); // call again 
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
		delete pcs;
	}

	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 
	gs.recover(vColor, mSubColor, mSimpl2Orig);

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
	RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();

	return acc_obj_value;
}

void SimpleMPL::construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt,
	SimpleMPL::graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor) const
{
	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn + i);
		vColor[i] = m_db->vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
	}

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
					boost::put(boost::edge_weight, dg, e.first, 1);
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
	for (std::vector<uint32_t>::const_iterator it = itBgn + 1; it != itEnd; it++)
	{
		uint32_t v1 = *(it - 1), v2 = *it;
		mplAssert(m_vCompId[v1] == m_vCompId(v2);
	}
#endif

	uint32_t acc_obj_value = std::numeric_limits<uint32_t>::max();
	// if current pattern does not contain uncolored patterns, directly calculate conflicts 
	if (check_uncolored(itBgn, itEnd))
		acc_obj_value = coloring_component(itBgn, itEnd, comp_id);

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
	if (acc_obj_value != std::numeric_limits<uint32_t>::max())
		mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);

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

										   // construct decomposition graph for component 
	construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor);

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
	std::ofstream out((filename + ".gv").c_str());
	boost::write_graphviz_dp(out, g, dp, string("id"));
	out.close();
	la::graphviz2pdf(filename);
}

void SimpleMPL::runProjection()
{
	if (m_db->vPatternBbox.empty())
	{
		mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}
#ifdef QDEBUG
	std::cout << "============= All original patterns ============" << std::endl;
	for (uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
	{
		rectangle_pointer_type temp = m_db->vPatternBbox[i];
		std::cout << i <<"\t" ; 
		std::cout << temp->pattern_id() << " -- " << gtl::xl(*temp) << "  " << gtl::yl(*temp) << "  " << gtl::xh(*temp) << "  " << gtl::yh(*temp) << std::endl;
	}
	std::cout << "================================================\n\n\n " << std::endl;
#endif

	uint32_t vertex_num = m_db->vPatternBbox.size();

	m_mAdjVertex.resize(vertex_num);
	
	uint32_t edge_num = construct_graph_from_distance(vertex_num);
	edge_num = edge_num >> 1;

	mplPrint(kINFO, "%u vertices, %u edges\n", vertex_num, edge_num);

	uint32_t num_new_pattern;

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:num_new_pattern)
#endif
	std::vector<rectangle_pointer_type> rect_vec = m_db->polyrect_patterns();
	std::vector<uint32_t> poly_rect_begin = m_db->polyrectBgnId();
	assert(poly_rect_begin.size() == vertex_num);

	// generate poly_rect_end, which stores the end index when querying from parent polygon id
	std::vector<uint32_t> poly_rect_end;
	poly_rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		poly_rect_end[i] = poly_rect_begin[i + 1] - 1;
	poly_rect_end[vertex_num - 1] = rect_vec.size() - 1;


	// generate m_mSplitMappingBbox, which stores all the newly-generated rectangles corresponding to original rectangles
	uint32_t rect_num = rect_vec.size();
	std::vector<std::vector<rectangle_pointer_type> > m_mSplitPatternBbox;
	m_mSplitPatternBbox.resize(rect_num);

	for (uint32_t v = 0; v < vertex_num; v++)
	{
		// polygon v
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[v];
		// polygon v's neighbor polygons
		std::vector<uint32_t> & nei_Vec = m_mAdjVertex[v];

		// use poss_nei_vec to obtain all the possible neighbor rectangles from neighbor polygons
		std::vector<rectangle_pointer_type> poss_nei_vec;
		for (std::vector<uint32_t>::iterator it = nei_Vec.begin(); it != nei_Vec.end(); it++)
		{
			uint32_t s_idx = poly_rect_begin[*it];
			uint32_t e_idx = poly_rect_end[*it];
			for (uint32_t a = s_idx; a <= e_idx; a++)
				poss_nei_vec.push_back(rect_vec[a]);
		}

		uint32_t start_idx = poly_rect_begin[v];
		uint32_t end_idx = poly_rect_end[v];
		
		// traverse all the rectangles in polygon v, to generate the intersections
		for (uint32_t j = start_idx; j <= end_idx; j++)
		{
			rectangle_type rect(*rect_vec[j]);
			// the generated split patterns
			std::vector<rectangle_pointer_type>& split = m_mSplitPatternBbox[j];
			projection(rect, split, poss_nei_vec);
			split.swap(split);
			num_new_pattern += split.size();
		}
	}

	// std::cout << "vertex_num now : " << vertex_num << "\n\n\n" << std::endl;

	std::vector<rectangle_pointer_type>().swap(m_db->vPatternBbox);
	
#ifdef DEPRECATED
	std::vector<uint32_t>().swap(new2ori);
	SplitMapping.resize(vertex_num);
	new2ori.resize(num_new_pattern);
	std::vector<uint32_t>().swap(m_vVertexOrder);
#endif

	uint32_t pattern_id = 0;
	for (uint32_t v = 0; v < rect_num; v++)
	{
		for (uint32_t j = 0; j < m_mSplitPatternBbox[v].size(); j++)
		{	
			// SplitMapping[v].push_back(pattern_id);
			// new2ori[pattern_id] = v;
			// m_vVertexOrder.push_back(pattern_id);
			m_mSplitPatternBbox[v][j]->pattern_id(pattern_id);
			m_db->vPatternBbox.push_back(m_mSplitPatternBbox[v][j]);
			pattern_id++;
		}
	}
#ifdef QDEBUG
	// std::cout << "============= All new patterns ============" << std::endl;
	//std::cout << "total patterns : " << num_new_pattern <<std::endl;
	/*
	for (uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
	{
		rectangle_pointer_type temp = m_db->vPatternBbox[i];
		std::cout << i << " : " <<  temp->pattern_id() << " -- " << gtl::xl(*temp) << "  " << gtl::yl(*temp) << "  " << gtl::xh(*temp) << "  " << gtl::yh(*temp) << std::endl;
	}
	*/
#endif

	mplPrint(kINFO, "Now it has %u new patterns.\n", num_new_pattern);
	m_vConflict.clear();
	m_mAdjVertex.clear();

#ifdef DEPRECATD
	std::vector<std::pair<uint32_t, uint32_t> >().swap(m_vConflict);
	std::vector<std::vector<uint32_t> > new_mAdjVertex;
	new_mAdjVertex.resize(pattern_id);

	adj4NewPatterns(m_mSplitPatternBbox, new_mAdjVertex);

	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex.resize(new_mAdjVertex.size());
	std::copy(m_mAdjVertex.begin(), m_mAdjVertex.end(), new_mAdjVertex.begin());

	std::vector<rectangle_pointer_type>().swap(m_db->vPatternBbox);
	for (uint32_t v = 0; v < vertex_num; v++)
	{
		for (uint32_t j = 0; j < m_mSplitPatternBbox[v].size(); j++)
			m_db->vPatternBbox.push_back(m_mSplitPatternBbox[v][j]);
	}
	std::cout << "end of runprojection!"<<std::endl;

	m_vCompId.resize(m_db->vPatternBbox.size(), std::numeric_limits<uint32_t>::max());
#endif
	

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

	uint32_t nei_num = nei_Vec.size();

	std::vector<coordinate_type> vstitches;
	
	GenerateStitchPosition_Bei(pRect, vInterSect, vPossibleStitches, nei_num, vstitches);
	//GenerateStitchPosition_Jian(pRect, vInterSect, vPossibleStitches, nei_num, vstitches);

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
#ifdef QDEBUG

	for(uint32_t i = 0;i < vStages.size(); i++)
	{
		std::cout << vStages[i].first.first << " -- " << vStages[i].first.second << "  :  " << vStages[i].second.size()<< std::endl;
	}
#endif
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		std::set<uint32_t> twosideSet;
		for (uint32_t j = 0; j <= i; j++)
			twosideSet.insert(vStages[j].second.begin(), vStages[j].second.end());
		if (twosideSet.size() >= 3 )// >= nei_num)
			continue;
		std::set<uint32_t>().swap(twosideSet);
		for (uint32_t j = i; j < vStages.size(); j++)
			twosideSet.insert(vStages[i].second.begin(), vStages[i].second.end());
		if (twosideSet.size() >= 3)// nei_num)
			continue;

		//std::cout << "Insert one stitch" << std::endl;
		vstitches.push_back((vStages[i].first.first + vStages[i].first.second) / 2);
	}
	sort(vstitches.begin(), vstitches.end());
#ifdef QDEBUG
	/*
	std::cout << "stitch position : " << std::endl;
	for(uint32_t i = 0; i < vstitches.size(); i++)
		std::cout << vstitches[i] << " ";
	std::cout << std::endl;
*/
#endif
	return;
}

void SimpleMPL::GenerateStitchPosition_Bei(const rectangle_type pRect, std::vector<rectangle_type> vInterSect, std::vector<coordinate_type> & vPossibleStitches,
	uint32_t nei_num, std::vector<coordinate_type> & vstitches)
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
#ifdef QDEBUG
		//	std::cout << "i = " << i << " \t" << "pos1 = " << pos1 << " \t" << "pos2 = " << pos2 << std::endl;
#endif
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

void SimpleMPL::adj4NewPatterns(std::vector<std::vector<rectangle_pointer_type> > & m_mSplitPatternBbox, std::vector<std::vector<uint32_t> > & new_mAdjVertex)
{
	uint32_t edge_num = 0;
	uint32_t vertex_num = 0;
	std::cout <<"in adj4" << std::endl;
	for (uint32_t i = 0; i < m_mSplitPatternBbox.size(); i++)
	{
		std::vector<rectangle_pointer_type> &  box_vec = m_mSplitPatternBbox[i];
		std::vector<uint32_t> & original_nel_vec = m_mAdjVertex[i];
		/*
		for(uint32_t j = 0; j < box_vec.size(); j++)
			std::cout << box_vec[j]->pattern_id() << " ";
		std::cout << std::endl;
		for(uint32_t j = 0; j < original_nel_vec.size(); j++)
			std::cout << "original nei : " << original_nel_vec[j] << std::endl;
		*/
		for(std::vector<rectangle_pointer_type>::iterator new_Pattern = box_vec.begin(); new_Pattern < box_vec.end(); new_Pattern++)
		{
			vertex_num ++;
		//	std::cout << "===== child " << (*new_Pattern)->pattern_id() << " ====="<<std::endl;
			for(std::vector<uint32_t>::iterator original_nei = original_nel_vec.begin(); original_nei != original_nel_vec.end(); original_nei++)
			{
		//		std::cout << "now for its original nei " << *original_nei << "  total " << m_mSplitPatternBbox.size()<< std::endl;
				std::vector<rectangle_pointer_type> & ori_nei_children_vec = m_mSplitPatternBbox[*original_nei];
				for(std::vector<rectangle_pointer_type>::iterator child = ori_nei_children_vec.begin(); child!=ori_nei_children_vec.end(); child++)
				{
					/*
					std::cout << (*new_Pattern)->pattern_id() ;
					std::cout << " -- ori  " << *original_nei ;
					std::cout << " -- child  " <<  (*child)->pattern_id() << std::endl;
					std::cout << "new_Pattern " <<  gtl::xl(**new_Pattern) << " " << gtl::yl(**new_Pattern) << " " << gtl::xh(**new_Pattern) << " " << gtl::yh(**new_Pattern) << std::endl;
					std::cout << "child : " << gtl::xl(**child) << " " << gtl::yl(**child) << " " << gtl::xh(**child) << " " << gtl::yh(**child) << std::endl;
					*/
					coordinate_difference distance = boost::geometry::distance(**new_Pattern, **child);
					// std::cout << " distance " << distance <<std::endl;
					if(distance < m_db->coloring_distance)
					{
						edge_num ++;
						new_mAdjVertex[(*new_Pattern)->pattern_id()].push_back((*child)->pattern_id());
					}
				}
			}
		}
	}
	mplPrint(kINFO, "%u vertices, %u edges\n", vertex_num, edge_num);
	std::cout <<"end of  adj4" << std::endl;
	return;
}

void SimpleMPL::print_welcome() const
{
	mplPrint(kNONE, "\n\n");
	mplPrint(kNONE, "=======================================================================\n");
	mplPrint(kNONE, "                        OpenMPL - Version 1.1                        \n");
	mplPrint(kNONE, "                                by                                   \n");
	mplPrint(kNONE, "                Yibo Lin, Bei Yu, Qi Sun and  David Z. Pan           \n");
	mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
	mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
	mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
	mplPrint(kNONE, "            Contact Authors:  {yibolin, dpan}@cerc.utexas.edu        \n");
	mplPrint(kNONE, "                              {byu, qsun}@cse.cuhk.edu.hk            \n");
	mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE
