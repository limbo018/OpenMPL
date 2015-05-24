/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL.h"
#include <stack>
#include <deque>
#include <boost/graph/graphviz.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>
#include <limbo/algorithms/coloring/LPColoring.h>
#include <limbo/algorithms/coloring/ILPColoring.h>

namespace SimpleMPL {

using std::cout;
using std::endl;
using std::stack;
using std::deque;

void SimpleMPL::run(int argc, char** argv)
{
	this->read_cmd(argc, argv);
	this->read_gds();
	this->solve();
	this->report();
	this->write_gds();
}
void SimpleMPL::read_cmd(int argc, char** argv)
{
	// read command 
	CmdParser<coordinate_type> cmd (m_db);
	assert_msg(cmd(argc, argv), "failed to parse command");
}
void SimpleMPL::read_gds()
{
	// read input gds file 
	GdsReader<coordinate_type> reader (m_db);
	assert_msg(reader(m_db.input_gds), "failed to read " << m_db.input_gds);
}
void SimpleMPL::write_gds()
{
	// write output gds file 
	GdsWriter<coordinate_type> writer;
	writer(m_db.output_gds, m_db);
}
void SimpleMPL::solve()
{
	this->construct_graph();
	this->connected_component();

	// create bookmark to index the starting position of each component 
	vector<uint32_t> vBookmark (m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i-1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

#pragma omp parallel for
	for (uint32_t comp_id = 0; comp_id != m_comp_cnt; ++comp_id)
	{
		// construct a component 
		vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+vBookmark[comp_id];
		vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+vBookmark[comp_id+1] : m_vVertexOrder.end();
		// solve component 
		// pass iterators to save memory 
		this->solve_component(itBgn, itEnd);
	}
}
void SimpleMPL::report() const 
{
	cout << "(I) Conflict number = " << conflict_num() << endl;
}

void SimpleMPL::construct_graph()
{
	// construct vertices 
	// assume vertices start from 0 and end with vertex_num-1
	uint32_t vertex_num = m_db.vPattern.size();
	m_vVertexOrder.resize(vertex_num, std::numeric_limits<uint32_t>::max()); 

	// construct edges 
	m_mAdjVertex.resize(vertex_num);
	if (m_db.hPath.empty()) // construct from distance 
	{
		for (uint32_t v = 0; v != vertex_num; ++v)
		{
			rectangle_pointer_type const& pPattern = m_db.vPattern[v];
			vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];

			// find patterns connected with pPattern 
			// query tPattern in m_db
			vector<rectangle_pointer_type> vAdjPattern;
			rectangle_type rect (*pPattern);
			// bloat pPattern with minimum coloring distance 
			gtl::bloat(rect, gtl::HORIZONTAL, m_db.coloring_distance);
			gtl::bloat(rect, gtl::VERTICAL, m_db.coloring_distance);
			for (rtree_type::const_query_iterator itq = m_db.tPattern.qbegin(bgi::intersects(rect));
					itq != m_db.tPattern.qend(); ++itq)
			{
				rectangle_pointer_type const& pAdjPattern = *itq;
				if (pAdjPattern != pPattern) // skip pPattern itself 
				{
					assert(pAdjPattern->pattern_id() != pPattern->pattern_id());
					// we consider euclidean distance
					gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = gtl::euclidean_distance(*pAdjPattern, *pPattern);
					if (distance < m_db.coloring_distance)
						vAdjVertex.push_back(pAdjPattern->pattern_id());
				}
			}
#if 0
			m_db.tPattern.query(bgi::intersects(rect), std::back_inserter(vAdjPattern));
			vAdjVertex.reserve(vAdjPattern.size()); // reserve enough memory 
			for (vector<rectangle_pointer_type>::iterator it = vAdjPattern.begin(); 
					it != vAdjPattern.end(); ++it)
			{
				// add pAdjPattern to stack 
				rectangle_pointer_type const& pAdjPattern = *it;
				if (pAdjPattern != pPattern) // skip pPattern itself 
				{
					assert(pAdjPattern->pattern_id() != pPattern->pattern_id());
					// we consider euclidean distance
					gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = gtl::euclidean_distance(*pAdjPattern, *pPattern);
					if (distance < m_db.coloring_distance)
						vAdjVertex.push_back(pAdjPattern->pattern_id());
				}
			}
#endif
			vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
		}
	}
	else // construct from conflict edges in hPath
	{
		// at the same time, estimate a coloring distance from conflict edges 
		m_db.coloring_distance = 0;
		for (set<int32_t>::const_iterator itLayer = m_db.sPathLayer.begin(); itLayer != m_db.sPathLayer.end(); ++itLayer)
		{
			if (!m_db.hPath.count(*itLayer)) continue;

			vector<path_type> const& vPath = m_db.hPath[*itLayer];

			for (vector<path_type>::const_iterator itPath = vPath.begin(); itPath != vPath.end(); ++itPath)
			{
				path_type const& path = *itPath;

				point_type const& p1 = path.low();
				point_type const& p2 = path.high();
				vector<rectangle_pointer_type> vPattern1;
				m_db.tPattern.query(bgi::contains(p1), std::back_inserter(vPattern1));
				vector<rectangle_pointer_type> vPattern2;
				m_db.tPattern.query(bgi::contains(p2), std::back_inserter(vPattern2));
#ifdef DEBUG
				if (vPattern1.size() > 1)
				{
					for (uint32_t i = 0; i != vPattern1.size(); ++i)
						printf("%d, %d, %d, %d, %d\n", vPattern1[i]->layer(), gtl::xl(*vPattern1[i]), gtl::yl(*vPattern1[i]), gtl::xh(*vPattern1[i]), gtl::yh(*vPattern1[i]));
				}
				if (vPattern2.size() > 1)
				{
					for (uint32_t i = 0; i != vPattern2.size(); ++i)
						printf("%d, %d, %d, %d, %d\n", vPattern2[i]->layer(), gtl::xl(*vPattern2[i]), gtl::yl(*vPattern2[i]), gtl::xh(*vPattern2[i]), gtl::yh(*vPattern2[i]));
				}
#endif
				assert(vPattern1.size() == 1);
				assert(vPattern2.size() == 1);
				rectangle_pointer_type const& pPattern1 = vPattern1.front();
				rectangle_pointer_type const& pPattern2 = vPattern2.front();

				// if the input gds file contains duplicate edges 
				// we will have duplicate edges here 
				m_mAdjVertex[pPattern1->pattern_id()].push_back(pPattern2->pattern_id());
				m_mAdjVertex[pPattern2->pattern_id()].push_back(pPattern1->pattern_id());

				// estimate coloring distance 
				coordinate_difference distance = gtl::length(path);
				m_db.coloring_distance = std::max(distance, m_db.coloring_distance);
			}
		}
		// eliminate all duplicates in adjacency list 
		for (uint32_t v = 0; v != vertex_num; ++v)
		{
			vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];
			set<uint32_t> sAdjVertex (vAdjVertex.begin(), vAdjVertex.end());
			vAdjVertex.resize(sAdjVertex.size()); // shrink to fit 
			vAdjVertex.assign(sAdjVertex.begin(), sAdjVertex.end()); 
		}
	}
	m_vCompId.resize(m_vVertexOrder.size(), std::numeric_limits<uint32_t>::max());
}

void SimpleMPL::connected_component()
{
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
		assert(m_vCompId[v] != std::numeric_limits<uint32_t>::max()); 
#endif

#if 0
	// reorder with order_id 
	// to save memory, use an approach very like swap which takes O(nlogn)
	for (uint32_t v = 0, order = m_vVertexOrder[v]; v != vertex_num; )
	{
		if (order != std::numeric_limits<uint32_t>::max() && v != order)
		{
			uint32_t tmp_v = order;
			uint32_t tmp_order = m_vVertexOrder[tmp_v];
			m_vVertexOrder[order] = v; // move v to position order 
			m_vVertexOrder[v] = std::numeric_limits<uint32_t>::max(); // set position v to empty
			v = tmp_v;
			order = tmp_order;
			if (order == std::numeric_limits<uint32_t>::max()) // if we replaced an empty position, go back to 0 
			{v = 0; order = m_vVertexOrder[0];}
			else assert(v != order);
		}
		else {++v; order = m_vVertexOrder[v];} 
	}
#else 
	vector<uint32_t> vTmpOrder (m_vVertexOrder.size());
	for (uint32_t v = 0; v != vertex_num; ++v)
		vTmpOrder[m_vVertexOrder[v]] = v;
	// it may be better to use move in c++11 
	std::copy(vTmpOrder.begin(), vTmpOrder.end(), m_vVertexOrder.begin());
#endif

#ifdef DEBUG
	// check ordered 
	for (uint32_t v = 1; v != vertex_num; ++v)
		assert(m_vCompId[m_vVertexOrder[v-1]] <= m_vCompId[m_vVertexOrder[v]]);
#endif
}

void SimpleMPL::depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& order_id)
{
	stack<uint32_t> vStack; 
	vStack.push(source);

	while (!vStack.empty())
	{
		uint32_t current = vStack.top();
		vStack.pop();
		if (m_vCompId[current] == std::numeric_limits<uint32_t>::max()) // not visited 
		{
			m_vCompId[current] = comp_id; // set visited 
			m_vVertexOrder[current] = order_id++; // update position 

			for (vector<uint32_t>::const_iterator it = m_mAdjVertex[current].begin(); it != m_mAdjVertex[current].end(); ++it)
			{
				uint32_t const& child = *it;
				assert(current != child);
				vStack.push(child);
			}
		}
	}
}

void SimpleMPL::solve_component(const vector<uint32_t>::const_iterator itBgn, const vector<uint32_t>::const_iterator itEnd)
{
	if (itBgn == itEnd) return;
	vector<rectangle_pointer_type>& vPattern = m_db.vPattern;
#ifdef DEBUG
	for (vector<uint32_t>::const_iterator it = itBgn+1; it != itEnd; ++it)
	{
		uint32_t v1 = *(it-1), v2 = *it;
		assert(m_vCompId[v1] == m_vCompId[v2]);
	}
#endif

	// construct a graph for current component 
	using namespace boost;
	// do not use setS, it does not compile for subgraph
	// do not use custom property tags, it does not compile for most utilities
	typedef adjacency_list<vecS, vecS, undirectedS, 
			property<vertex_index_t, uint32_t, property<vertex_color_t, int> >, 
			property<edge_index_t, uint32_t, property<edge_weight_t, int> >,
			property<graph_name_t, string> > graph_type;
	typedef property<vertex_index_t, uint32_t> VertexId;
	typedef property<edge_index_t, uint32_t> EdgeID;
	typedef typename graph_traits<graph_type>::vertex_descriptor vertex_descriptor; 
	typedef typename graph_traits<graph_type>::edge_descriptor edge_descriptor;
	typedef typename property_map<graph_type, edge_weight_t>::type edge_weight_map_type;

	uint32_t pattern_cnt = itEnd-itBgn;
	// decomposition graph 
	graph_type dg (pattern_cnt);
	vector<int8_t> vColor (pattern_cnt, -1); // coloring results 
	map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		vColor[i] = vPattern[v]->color();
		mGlobal2Local[v] = i;
	}

	// edges 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
		for (vector<uint32_t>::const_iterator it = m_mAdjVertex[v].begin(); it != m_mAdjVertex[v].end(); ++it)
		{
			uint32_t u = *it;
#ifdef DEBUG
			assert(mGlobal2Local.count(u));
#endif
			uint32_t j = mGlobal2Local[u];
			if (i < j) // avoid duplicate 
			{
				pair<edge_descriptor, bool> e = edge(i, j, dg);
				if (!e.second) // make sure no duplicate 
				{
					e = add_edge(i, j, dg);
					assert(e.second);
					put(edge_weight, dg, e.first, 1);
				}
			}
		}
	}

#ifdef DEBUG
	if (0)
	{
		boost::dynamic_properties dp;
		dp.property("id", boost::get(boost::vertex_index, dg));
		dp.property("node_id", boost::get(boost::vertex_index, dg));
		dp.property("label", boost::get(boost::vertex_index, dg));
		dp.property("weight", boost::get(boost::edge_weight, dg));
		dp.property("label", boost::get(boost::edge_weight, dg));
		ofstream out ("graph_init.gv");
		boost::write_graphviz_dp(out, dg, dp, string("id"));
		out.close();
		system("dot -Tpdf graph_init.gv -o graph_init.pdf");
	}
#endif
	// graph simplification 
	typedef limbo::algorithms::coloring::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (dg);
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
	// keep the order of simplification 
	gs.hide_small_degree(m_db.color_num); // hide vertices with degree smaller than color_num
	if (m_db.color_num == 3)
		gs.merge_subK4(); // merge sub-K4 structure 
	// collect simplified information 
	vector<graph_simplification_type::vertex_status_type> const& vStatus = gs.status();
	vector<vector<vertex_descriptor> > const& vChildren = gs.children();
	stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();
	pair<graph_type, map<vertex_descriptor, vertex_descriptor> > sg = gs.simplified_graph(); // simplified graph and vertex mapping

#ifdef DEBUG
	//gs.write_simplified_graph_dot("graph_simpl");
#endif

	// solve coloring 
#if 0
	typedef limbo::algorithms::coloring::LPColoring<graph_type> coloring_solver_type;
	coloring_solver_type cs (sg.first);
	cs.stitchWeight(0.1);
	cs.conflictCost(false);
	cs.stitchMode(false);
	cs.roundingScheme(coloring_solver_type::DIRECT_ILP);
	// THREE or FOUR 
	if (m_db.color_num == 3)
		cs.colorNum(coloring_solver_type::THREE);
	else 
		cs.colorNum(coloring_solver_type::FOUR);
	cs(); // solve coloring 

	// collect coloring results from simplified graph 
	{
		graph_traits<graph_type>::vertex_iterator vi, vie;
		for (tie(vi, vie) = vertices(sg.first); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = cs.cg_vertexColor(v);
			assert(color >= 0 && color < m_db.color_num);
			vColor[ sg.second[v] ] = color;
		}
	}
#else 
	typedef limbo::algorithms::coloring::ILPColoring<graph_type> coloring_solver_type;
	coloring_solver_type cs (sg.first);
	cs.stitch_weight(0.1);
	cs.color_num(m_db.color_num);
	// set precolored vertices 
	graph_traits<graph_type>::vertex_iterator vi, vie;
	for (tie(vi, vie) = vertices(sg.first); vi != vie; ++vi)
	{
		vertex_descriptor v = *vi;
		int8_t color = vColor[sg.second[v]];
		if (color >= 0 && color < m_db.color_num)
			cs.precolor(v, color);
	}
	cs(); // solve coloring 

	// collect coloring results from simplified graph 
	for (tie(vi, vie) = vertices(sg.first); vi != vie; ++vi)
	{
		vertex_descriptor v = *vi;
		int8_t color = cs.color(v);
		assert(color >= 0 && color < m_db.color_num);
		vColor[ sg.second[v] ] = color;
	}
#endif
	
	// recover colors for simplified vertices with balanced assignment 
	// recover merged vertices 
	for (uint32_t v = 0; v != vStatus.size(); ++v)
	{
		if (vStatus[v] == graph_simplification_type::GOOD)
		{
			assert(vColor[v] >= 0 && vColor[v] < m_db.color_num);
			for (uint32_t j = 0; j != vChildren[v].size(); ++j)
			{
				vertex_descriptor u = vChildren[v][j];
				if (v != u) 
					vColor[u] = vColor[v];
			}
		}
	}
	// recover hidden vertices with local balanced density control 
	while (!vHiddenVertices.empty())
	{
		vertex_descriptor v = vHiddenVertices.top();
		vHiddenVertices.pop();

		// find available colors 
		deque<bool> vUnusedColor (m_db.color_num, true);
		graph_traits<graph_type>::adjacency_iterator vi, vie;
		for (tie(vi, vie) = adjacent_vertices(v, dg); vi != vie; ++vi)
		{
			vertex_descriptor u = *vi;
			if (vColor[u] >= 0)
			{
				assert(vColor[u] < m_db.color_num);
				vUnusedColor[u] = false;
			}
		}

		// find the nearest distance of each color 
		// search all patterns in the component 
		// TO DO: further speedup is possible to search a local window 
		vector<coordinate_difference> vDist (m_db.color_num, std::numeric_limits<coordinate_difference>::max());
		for (uint32_t u = 0; u != pattern_cnt; ++u)
		{
			if (v == u) continue;
			// skip uncolored vertices 
			if (vColor[u] < 0) continue; 
			// we consider euclidean distance
			gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = gtl::euclidean_distance(*vPattern[*(itBgn+v)], *vPattern[*(itBgn+u)]);
			vDist[ vColor[u] ] = std::min(vDist[ vColor[u] ], distance);
		}

		// choose the color with largest distance 
		int8_t best_color = -1;
		coordinate_difference best_dist = std::numeric_limits<coordinate_difference>::min();
		for (int8_t i = 0; i != m_db.color_num; ++i)
		{
			if (vUnusedColor[i])
			{
				if (best_dist < vDist[i])
				{
					best_color = i;
					best_dist = vDist[i];
				}
			}
		}
		assert(best_color >= 0 && best_color < m_db.color_num);
		vColor[v] = best_color;
	}

#ifdef DEBUG
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		assert(vColor[i] >= 0 && vColor[i] < m_db.color_num);
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		if (vPattern[v]->color() >= 0 && vPattern[v]->color() < m_db.color_num) // check precolored pattern 
			assert(vPattern[v]->color() == vColor[i]);
		else // assign color to uncolored pattern 
			vPattern[v]->color(vColor[i]);
	}
}

uint32_t SimpleMPL::conflict_num(const vector<uint32_t>::const_iterator itBgn, const vector<uint32_t>::const_iterator itEnd) const
{
	vector<rectangle_pointer_type> const& vPattern = m_db.vPattern;
	uint32_t pattern_cnt = itEnd-itBgn;
	uint32_t cnt = 0;
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
		int8_t color1 = vPattern[v]->color();
		if (color1 >= 0 && color1 < m_db.color_num)
		{
			for (vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				int8_t color2 = vPattern[u]->color();
				if (color2 >= 0 && color2 < m_db.color_num)
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
	vector<rectangle_pointer_type> const& vPattern = m_db.vPattern;
	uint32_t cnt = 0;
	for (uint32_t v = 0; v != vPattern.size(); ++v)
	{
		int8_t color1 = vPattern[v]->color();
		if (color1 >= 0 && color1 < m_db.color_num)
		{
			for (vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				int8_t color2 = vPattern[u]->color();
				if (color2 >= 0 && color2 < m_db.color_num)
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

} // namespace SimpleMPL

