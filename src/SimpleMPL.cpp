/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"

#include <stack>
#include <boost/graph/graphviz.hpp>
#include <boost/timer/timer.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>

// only valid when gurobi is available 
#if GUROBI == 1
#include <limbo/algorithms/coloring/ILPColoring.h>
#include <limbo/algorithms/coloring/LPColoring.h>
#endif
// only valid when lemon cbc api is available 
#if LEMONCBC == 1
#include <limbo/algorithms/coloring/ILPColoringLemonCbc.h>
#endif
#include <limbo/algorithms/coloring/BacktrackColoring.h>

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
	this->report();
	this->write_gds();
}
void SimpleMPL::reset(bool init)
{
    // release memory and set to initial value 
    if (!init && m_db) delete m_db;
    m_db = NULL;
    m_comp_cnt = 0;
    std::vector<uint32_t>().swap(m_vVertexOrder);
    std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
    std::vector<uint32_t>().swap(m_vCompId);
    std::vector<uint32_t>().swap(m_vColorDensity);
    std::vector<std::pair<uint32_t, uint32_t> >().swap(m_vConflict);
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
	boost::timer::auto_cpu_timer timer (buf);
	if (m_db->vPatternBbox.empty())
	{
        mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}

	this->construct_graph();
	this->connected_component();

	// create bookmark to index the starting position of each component 
	std::vector<uint32_t> vBookmark (m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i-1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

	mplPrint(kINFO, "Solving %u independent components...\n", m_comp_cnt);
	// thread number controled by user option 
#ifdef _OPENMP
#pragma omp parallel num_threads (m_db->thread_num())
#endif
	{
#ifdef _OPENMP
#pragma omp for
#endif
		for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
		{
			// construct a component 
			std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+vBookmark[comp_id];
			std::vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+vBookmark[comp_id+1] : m_vVertexOrder.end();
#ifdef DEBUG
//					if (comp_id != 9941)
//						continue;
#endif
			// solve component 
			// pass iterators to save memory 
			this->solve_component(itBgn, itEnd, comp_id);
		}
	}
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
	{
		for (uint32_t v = 0; v != vertex_num; ++v)
		{
			rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
			std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];

			// find patterns connected with pPattern 
			// query tPatternBbox in m_db
			std::vector<rectangle_pointer_type> vAdjPattern;
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
			edge_num += vAdjVertex.size();
		}
	}
	else // construct from conflict edges in hPath
	{
		// at the same time, estimate a coloring distance from conflict edges 
		m_db->coloring_distance = 0;
		for (set<int32_t>::const_iterator itLayer = m_db->parms.sPathLayer.begin(), itLayerE = m_db->parms.sPathLayer.end(); itLayer != itLayerE; ++itLayer)
		{
			if (!m_db->hPath.count(*itLayer)) continue;

			std::vector<path_type> const& vPath = m_db->hPath[*itLayer];

			for (std::vector<path_type>::const_iterator itPath = vPath.begin(); itPath != vPath.end(); ++itPath)
			{
				path_type const& path = *itPath;

				point_type const& p1 = path.low();
				point_type const& p2 = path.high();
				std::vector<rectangle_pointer_type> vPattern1;
				m_db->tPatternBbox.query(bgi::contains(p1), std::back_inserter(vPattern1));
				std::vector<rectangle_pointer_type> vPattern2;
				m_db->tPatternBbox.query(bgi::contains(p2), std::back_inserter(vPattern2));
#ifdef DEBUG
				if (vPattern1.size() > 1)
				{
					for (uint32_t i = 0; i != vPattern1.size(); ++i)
						mplPrint(kDEBUG, "%d, %d, %d, %d, %d\n", vPattern1[i]->layer(), gtl::xl(*vPattern1[i]), gtl::yl(*vPattern1[i]), gtl::xh(*vPattern1[i]), gtl::yh(*vPattern1[i]));
				}
				if (vPattern2.size() > 1)
				{
					for (uint32_t i = 0; i != vPattern2.size(); ++i)
						mplPrint(kDEBUG, "%d, %d, %d, %d, %d\n", vPattern2[i]->layer(), gtl::xl(*vPattern2[i]), gtl::yl(*vPattern2[i]), gtl::xh(*vPattern2[i]), gtl::yh(*vPattern2[i]));
				}
#endif
				mplAssert(vPattern1.size() == 1);
				mplAssert(vPattern2.size() == 1);
				rectangle_pointer_type const& pPattern1 = vPattern1.front();
				rectangle_pointer_type const& pPattern2 = vPattern2.front();

				// if the input gds file contains duplicate edges 
				// we will have duplicate edges here 
				m_mAdjVertex[pPattern1->pattern_id()].push_back(pPattern2->pattern_id());
				m_mAdjVertex[pPattern2->pattern_id()].push_back(pPattern1->pattern_id());

				// estimate coloring distance 
				coordinate_difference distance = gtl::length(path);
				m_db->coloring_distance = std::max(distance, m_db->coloring_distance);
			}
		}
		// eliminate all duplicates in adjacency list 
		for (uint32_t v = 0; v != vertex_num; ++v)
		{
			std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];
			set<uint32_t> sAdjVertex (vAdjVertex.begin(), vAdjVertex.end());
			vAdjVertex.resize(sAdjVertex.size()); // shrink to fit 
			vAdjVertex.assign(sAdjVertex.begin(), sAdjVertex.end()); 
			edge_num += vAdjVertex.size();
		}
		m_db->parms.coloring_distance_nm = m_db->coloring_distance*(m_db->unit*1e+9);
		mplPrint(kINFO, "Estimated coloring distance from conflict edges = %lld (%g nm)\n", m_db->coloring_distance, m_db->coloring_distance_nm());
	}
	m_vCompId.resize(m_vVertexOrder.size(), std::numeric_limits<uint32_t>::max());
	m_vColorDensity.assign(m_db->color_num(), 0);
	m_vConflict.clear();
    edge_num = edge_num>>1;

	// report statistics 
	mplPrint(kINFO, "%u vertices, %u edges\n", vertex_num, edge_num);
}

void SimpleMPL::connected_component()
{
	mplPrint(kINFO, "Computing connected components...\n");
	uint32_t comp_id = 0;
	uint32_t order_id = 0; // position in an ordered pattern array with grouped components 

	uint32_t vertex_num = m_vVertexOrder.size();
	// here we employ the assumption that the graph data structure alreays has vertices labeled with 0~vertex_num-1
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

// do not use setS, it does not compile for subgraph
// do not use custom property tags, it does not compile for most utilities
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
        boost::property<boost::vertex_index_t, uint32_t, boost::property<boost::vertex_color_t, int> >, 
        boost::property<boost::edge_index_t, uint32_t, boost::property<boost::edge_weight_t, int> >,
        boost::property<boost::graph_name_t, std::string> > graph_type;
//typedef property<vertex_index_t, uint32_t> VertexId;
//typedef property<edge_index_t, uint32_t> EdgeID;
typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor; 
typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;

///// helper functions for SimpleMPL::solve_component
///// since we may want to call solve_graph_coloring() multiple times 
///// it is better to wrap it as an independent function 
/////
///// TO DO: current way to wrap functions are not clean enough, too many parameters 
///// find a cleaner way to wrap functions

/// create coloring solver pointer according to algorithm type
limbo::algorithms::coloring::Coloring<graph_type>* create_coloring_solver(graph_type const& sg, AlgorithmType const& algo, int32_t color_num)
{
    typedef limbo::algorithms::coloring::Coloring<graph_type> coloring_solver_type;
    coloring_solver_type* pcs = NULL;
    switch (algo.get())
    {
#if GUROBI == 1
        case AlgorithmTypeEnum::ILP_GURBOI:
            pcs = new limbo::algorithms::coloring::ILPColoring<graph_type> (sg); break;
        case AlgorithmTypeEnum::LP_GUROBI:
            pcs = new limbo::algorithms::coloring::LPColoring<graph_type> (sg); break;
#endif
#if LEMONCBC == 1
        case AlgorithmTypeEnum::ILP_CBC:
            pcs = new limbo::algorithms::coloring::ILPColoringLemonCbc<graph_type> (sg); break;
#endif
        case AlgorithmTypeEnum::BACKTRACK:
            pcs = new limbo::algorithms::coloring::BacktrackColoring<graph_type> (sg);
            break;
        default: mplAssertMsg(0, "unknown algorithm type");
    }
    pcs->stitch_weight(0.1);
    pcs->color_num(color_num);
    pcs->threads(1); // we use parallel at higher level 

    return pcs;
}
/// recover color of vertices simplified by HIDE_SMALL_DEGREE
/// consider density balance 
void recover_hide_vertex_colors(graph_type const& dg, SimpleMPL::layoutdb_type const* db, 
        const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, 
        std::vector<int8_t>& vColor, std::vector<uint32_t>& vColorDensity, std::stack<vertex_descriptor>& vHiddenVertices)
{
    typedef SimpleMPL::coordinate_type coordinate_type;
    typedef SimpleMPL::coordinate_difference coordinate_difference;

	uint32_t pattern_cnt = itEnd-itBgn;

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
	while (!vHiddenVertices.empty())
	{
		vertex_descriptor v = vHiddenVertices.top();
		vHiddenVertices.pop();

		// find available colors 
        std::vector<char> vUnusedColor (db->color_num(), true);
        boost::graph_traits<graph_type>::adjacency_iterator vi, vie;
		for (tie(vi, vie) = adjacent_vertices(v, dg); vi != vie; ++vi)
		{
			vertex_descriptor u = *vi;
			if (vColor[u] >= 0)
			{
				mplAssert(vColor[u] < db->color_num());
				vUnusedColor[vColor[u]] = false;
			}
		}

		// find the nearest distance of each color 
		// search all patterns in the component 
		// TO DO: further speedup is possible to search a local window 
		std::vector<coordinate_difference> vDist (db->color_num(), std::numeric_limits<coordinate_difference>::max());
		for (uint32_t u = 0; u != pattern_cnt; ++u)
		{
			if (v == u) continue;
			// skip uncolored vertices 
			if (vColor[u] < 0) continue; 
			// we consider euclidean distance
            // use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
			gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = db->euclidean_distance(*db->vPatternBbox[*(itBgn+v)], *db->vPatternBbox[*(itBgn+u)]);
#ifdef DEBUG
			mplAssert(vColor[u] < db->color_num() && distance >= 0);
#endif
			vDist[ vColor[u] ] = std::min(vDist[ vColor[u] ], distance);
		}

		// choose the color with largest distance 
		int8_t best_color = -1;
		double best_score = -std::numeric_limits<double>::max(); // negative max 
		for (int8_t i = 0; i != db->color_num(); ++i)
		{
			if (vUnusedColor[i])
			{
				double cur_score = (double)vDist[i]/(1.0+vColorDensity[i]);
				if (best_score < cur_score)
				{
					best_color = i;
					best_score = cur_score;
				}
			}
		}
		mplAssert(best_color >= 0 && best_color < db->color_num());
		vColor[v] = best_color;
	}
}
/// given a graph, solve coloring 
/// contain nested call for itself 
uint32_t solve_graph_coloring(uint32_t comp_id, graph_type const& dg, SimpleMPL::layoutdb_type const* db, 
        const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, 
        uint32_t simplify_strategy, std::vector<int8_t>& vColor, std::vector<uint32_t>& vColorDensity)
{
	typedef limbo::algorithms::coloring::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (dg, db->color_num());
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
    if (db->color_num() == 3)
        gs.max_merge_level(3);
    else if (db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works
        gs.max_merge_level(2);
    gs.simplify(simplify_strategy);
	// collect simplified information 
    std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

    // for debug, it does not affect normal run 
	if (comp_id == db->dbg_comp_id() && simplify_strategy != graph_simplification_type::MERGE_SUBK4)
		gs.write_simplified_graph_dot("graph_simpl");

	// in order to recover color from articulation points 
	// we have to record all components and mappings 
	// but graph is not necessary 
	std::vector<std::vector<int8_t> > mSubColor (gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig (gs.num_component());
	double acc_obj_value = 0;
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);

		vSubColor.assign(num_vertices(sg), -1);

		// solve coloring 
		typedef limbo::algorithms::coloring::Coloring<graph_type> coloring_solver_type;
		coloring_solver_type* pcs = create_coloring_solver(sg, db->algo(), db->color_num());

		// set precolored vertices 
        boost::graph_traits<graph_type>::vertex_iterator vi, vie;
		for (tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = vColor[vSimpl2Orig[v]];
			if (color >= 0 && color < db->color_num())
            {
				pcs->precolor(v, color);
                vSubColor[v] = color; // necessary for 2nd trial
            }
		}
        // 1st trial 
		double obj_value1 = (*pcs)(); // solve coloring 
        // 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
        double obj_value2 = std::numeric_limits<double>::max();
        // very restric condition to determin whether perform MERGE_SUBK4 or not 
        if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && db->algo() == AlgorithmTypeEnum::LP_GUROBI 
                && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
            obj_value2 = solve_graph_coloring(comp_id, sg, db, itBgn, itEnd, graph_simplification_type::MERGE_SUBK4, vSubColor, vColorDensity); // call again 

        if (obj_value1 < obj_value2)
        {
            acc_obj_value += obj_value1;

            // collect coloring results from simplified graph 
            for (tie(vi, vie) = vertices(sg); vi != vie; ++vi)
            {
                vertex_descriptor v = *vi;
                int8_t color = pcs->color(v);
                mplAssert(color >= 0 && color < db->color_num());
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
    recover_hide_vertex_colors(dg, db, itBgn, itEnd, vColor, vColorDensity, vHiddenVertices);

    return acc_obj_value;
}

uint32_t SimpleMPL::solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
//{{{
{
	if (itBgn == itEnd) return 0;
	std::vector<rectangle_pointer_type>& vPatternBbox = m_db->vPatternBbox;
#ifdef DEBUG
	for (std::vector<uint32_t>::const_iterator it = itBgn+1; it != itEnd; ++it)
	{
		uint32_t v1 = *(it-1), v2 = *it;
		mplAssert(m_vCompId[v1] == m_vCompId[v2]);
	}
#endif

	// construct a graph for current component 
	uint32_t pattern_cnt = itEnd-itBgn;

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...\n", comp_id, pattern_cnt);

	// decomposition graph 
	graph_type dg (pattern_cnt);
	std::vector<int8_t> vColor (pattern_cnt, -1); // coloring results 
	map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		vColor[i] = vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
	}

	// edges 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
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

    // for debug, it does not affect normal run 
	if (comp_id == m_db->dbg_comp_id())
	{
		boost::dynamic_properties dp;
		dp.property("id", boost::get(boost::vertex_index, dg));
		dp.property("node_id", boost::get(boost::vertex_index, dg));
		dp.property("label", boost::get(boost::vertex_index, dg));
		dp.property("weight", boost::get(boost::edge_weight, dg));
		dp.property("label", boost::get(boost::edge_weight, dg));
        std::ofstream out ("graph_init.gv");
		boost::write_graphviz_dp(out, dg, dp, string("id"));
		out.close();
		::system("dot -Tpdf graph_init.gv -o graph_init.pdf");
	}
	// graph simplification 
	typedef limbo::algorithms::coloring::GraphSimplification<graph_type> graph_simplification_type;
    uint32_t simplify_strategy = graph_simplification_type::NONE;
	// keep the order of simplification 
	if (m_db->simplify_level() > 0)
        simplify_strategy |= graph_simplification_type::HIDE_SMALL_DEGREE;
	if (m_db->simplify_level() > 1)
        simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

    // solve graph coloring 
    uint32_t acc_obj_value = solve_graph_coloring(comp_id, dg, m_db, itBgn, itEnd, simplify_strategy, vColor, m_vColorDensity);

#ifdef DEBUG
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		mplAssert(vColor[i] >= 0 && vColor[i] < m_db->color_num());
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
        m_db->set_color(v, vColor[i]);
	}

	// update global color density map 
	// if parallelization is enabled, there will be uncertainty in the density map 
	// because the density is being updated while being read 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t& color_density = m_vColorDensity[vColor[i]];
#ifdef _OPENMP
#pragma omp atomic
#endif
		++color_density;
	}

	uint32_t component_conflict_num = conflict_num(itBgn, itEnd);
    // only valid under no stitch 
	mplAssert(acc_obj_value == component_conflict_num);

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, pattern_cnt, component_conflict_num);

	return component_conflict_num;
}
//}}}


uint32_t SimpleMPL::conflict_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd) const
{
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	uint32_t pattern_cnt = itEnd-itBgn;
	uint32_t cnt = 0;
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
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

void SimpleMPL::print_welcome() const
{
  mplPrint(kNONE, "\n\n");
  mplPrint(kNONE, "=======================================================================\n");
  mplPrint(kNONE, "                      SimpleMPL - Version 1.0                        \n");
  mplPrint(kNONE, "                                by                                   \n");  
  mplPrint(kNONE, "                   Yibo Lin, Bei Yu, and  David Z. Pan               \n");
  mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
  mplPrint(kNONE, "                         Copyright (c) 2015                          \n");
  mplPrint(kNONE, "            Contact Authors:  {yibolin,bei,dpan}@cerc.utexas.edu     \n");
  mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE

