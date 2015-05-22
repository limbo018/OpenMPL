/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL"
#include <stack>
#include <limbo/algorithms/coloring/GraphSimplification.h>
#include <limbo/algorithms/coloring/LPColoring.h>

namespace SimpleMPL {

using std::cout;
using std::endl;
using std::stack;

void SimpleMPL::run(int argc, char** argv)
{
	this->read_cmd(argc, argv);
	this->read_gds();
	this->write_gds();
	this->solve();
}
void SimpleMPL::read_cmd(int argc, char** argv)
{
	// read command 
	CmdParser cmd (m_db);
	assert(cmd(argc, argv));
}
void SimpleMPL::read_gds()
{
	// read input gds file 
	GdsReader reader (m_db);
	assert(reader(m_db.input_gds));
}
void SimpleMPL::write_gds()
{
	// write output gds file 
	GdsWriter writer;
	assert(writer.m_db.output_gds, m_db);
}
void SimpleMPL::solve()
{
	this->connected_component();

	// create bookmark to index the starting position of each component 
	vector<uint32_t> vBookmark (m_comp_cnt);
	for (size_t i = 0, comp_id = 0; i != m_db.vPattern.size(); ++i)
	{
		if (i == 0 || m_db.vPattern[i-1]->comp_id() != m_db.vPattern[i]->comp_id())
			vBookmark[m_db.vPattern[i]->comp_id()] = i;
	}

#pragma omp parallel for
	for (uint32_t comp_id = 0; comp_id != m_comp_cnt; ++comp_id)
	{
		// construct a component 
		vector<rectangle_pointer_type>::iterator itBgn = m_db.vPattern.begin()+vBookmark[comp_id];
		vector<rectangle_pointer_type>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_db.vPattern.begin()+vBookmark[comp_id+1] : m_db.vPattern.end();
		// solve component 
		// pass iterators to save memory 
		this->solve_component(itBgn, itEnd);
	}
}

void SimpleMPL::connected_component()
{
	uint32_t comp_id = 0;
	uint32_t pattern_id = 0; // position in an ordered pattern array with grouped components 
	for (vector<rectangle_pointer_type>::iterator it = m_db.vPattern.begin();
			it != m_db.vPattern.end(); ++it)
	{
		rectangle_pointer_type pPattern = *it;

		if (pPattern->comp_id() != std::numeric_limits<uint32_t>::max()) // not visited 
		{
			depth_first_search(pPattern, comp_id, pattern_id);
			comp_id += 1;
		}
	}
	// record maximum number of connected components
	m_comp_cnt = comp_id;

#ifdef DEBUG 
	// check visited 
	for (vector<rectangle_pointer_type>::iterator it = m_db.vPattern.begin();
			it != m_db.vPattern.end(); ++it)
	{
		rectangle_pointer_type pPattern = *it;
		assert(pPattern->comp_id() != std::numeric_limits<uint32_t>::max()); 
	}
#endif

	// reorder with pattern_id 
	// to save memory, use swap which takes O(nlogn)
	// if memory is not a problem, constructing a new array only takes O(n)
	for (size_t i = 0; i != m_db.vPattern.size(); )
	{
		rectangle_pointer_type& pPattern = m_db.vPattern[i];
		if (pPattern->pattern_id() != i) // not in the position desired 
			std::swap(pPattern, m_db.vPattern[pPattern->pattern_id()]);
		else ++i;
	}

#ifdef DEBUG
	// check ordered 
	for (vector<rectangle_pointer_type>::iterator it = m_db.vPattern.begin()+1;
			it != m_db.vPattern.end(); ++it)
	{
		assert((*(it-1))->pattern_id() < (*it)->pattern_id());
		assert((*(it-1))->comp_id() <= (*it)->comp_id());
	}
#endif
}

void SimpleMPL::depth_first_search(SimpleMPL::rectangle_pointer_type source, uint32_t comp_id, uint32_t& pattern_id)
{
	stack<rectangle_pointer_type> vStack; 
	vStack.push(source);

	while (!vStack.empty())
	{
		rectangle_pointer_type current = vStack.top();
		vStack.pop();
		if (current->comp_id() != std::numeric_limits<uint32_t>::max()) // not visited 
		{
			current->comp_id(comp_id); // set visited 
			current->pattern_id(pattern_id++); // update position 

			if (m_db.sPathLayer.empty()) // conflict edges from coloring_distance
			{
				// find patterns connected with current 
				// query tPattern in m_db
				vector<rectangle_pointer_type> vChildren;
				rectangle_type rect (*current);
				// bloat current with minimum coloring distance 
				gtl::bloat(rect, gtl::HORIZONTAL, m_db.coloring_distance);
				gtl::bloat(rect, gtl::VERTICAL, m_db.coloring_distance);
				m_db.tPattern.query(bgi::intersects(rect), std::back_inserter(vChildren));
				for (vector<rectangle_pointer_type>::iterator it = vChildren.begin(); 
						it != vChildren.end(); ++it)
				{
					// add child to stack 
					rectangle_pointer_type child = *it;
					if (child != current) // skip current itself 
					{
						assert(child->id() != current->id());
						// we consider eucliean distance
						gtl::coordinate_traits<coordinate_type>::coordinate_distance distance = gtl::eucliean_distance(*child, *current);
						if (distance < m_db.coloring_distance)
							vStack.push(child);
					}
				}
			}
			else // conflict edges from path layers 
			{
				// need rtree of path 
			}
		}
	}
}

void SimpleMPL::solve_component(const vector<rectangle_pointer_type>::const_iterator itBgn, const vector<rectangle_pointer_type>::const_iterator itEnd)
{
	if (itBgn == itEnd) return;
#ifdef DEBUG
	for (vector<rectangle_pointer_type>::const_iterator it = itBgn+1; it != itEnd; ++it)
		assert((*it)->comp_id() == (*(it-1))->comp_id());
#endif

	// construct a graph for current component 
	using namespace boost;
	// do not use setS, it does not compile for subgraph
	// do not use custom property tags, it does not compile for most utilities
	typedef adjacency_list<vecS, vecS, undirectedS, 
			property<vertex_index_t, std::size_t, property<vertex_color_t, int> >, 
			property<edge_index_t, std::size_t, property<edge_weight_t, int> >,
			property<graph_name_t, string> > graph_type;
	typedef property<vertex_index_t, std::size_t> VertexId;
	typedef property<edge_index_t, std::size_t> EdgeID;
	typedef typename graph_traits<graph_type>::vertex_descriptor vertex_descriptor; 
	typedef typename graph_traits<graph_type>::edge_descriptor edge_descriptor;
	typedef typename property_map<graph_type, edge_weight_t>::type edge_weight_map_type;

	size_t pattern_cnt = itEnd-itBgn;
	// decomposition graph 
	graph_type dg (pattern_cnt);
	vector<int8_t> vColor (pattern_cnt, -1); // coloring results 

	// precolored patterns 
	for (size_t i = 0; i != pattern_cnt; ++i)
		vColor[i] = (itBgn+i)->color();

	if (pattern_cnt < 30) // when there are very few patterns, construct graph by direct searching 
	{
		for (size_t i = 0; i != pattern_cnt; ++i)
			for (size_t j = i+1; j != pattern_cnt; ++j)
			{
				// we consider eucliean distance
				gtl::coordinate_traits<coordinate_type>::coordinate_distance distance = gtl::eucliean_distance(*(itBgn+i), *(itBgn+j));
				if (distance < m_db.coloring_distance)
				{
					pair<edge_descriptor, bool> e = add_edge(i, j, dg);
					assert(e.second);
					put(edge_weight, dg, e.first, 1);
				}
			}
	}
	else // build local rtree for fast searching  
	{
		rtree_type tPattern (itBgn, itEnd); 
		map<rectangle_pointer_type, size_t> mPattern2Idx; // map pattern to index in vPattern
		for (size_t i = 0; i != pattern_cnt; ++i)
			assert(mPattern2Idx.insert(std::make_pair((itBgn+i), i)).second);

		for (size_t i = 0; i != pattern_cnt; ++i)
		{
			rectangle_pointer_type pPattern1 = *(itBgn+i);

			vector<rectangle_pointer_type> vPattern2;
			rectangle_type rect (*pPattern1);
			// bloat current with minimum coloring distance 
			gtl::bloat(rect, gtl::HORIZONTAL, m_db.coloring_distance);
			gtl::bloat(rect, gtl::VERTICAL, m_db.coloring_distance);
			tPattern.query(bgi::intersects(rect), std::back_inserter(vPattern2));
			for (vector<rectangle_pointer_type>::iterator it = vPattern2.begin(); 
					it != vPattern2.end(); ++it)
			{
				rectangle_pointer_type pPattern2 = *it;
				if (pPattern2 != pPattern1) // skip pPattern1 itself 
				{
					assert(pPattern2->id() != pPattern1->id());
					// we consider eucliean distance
					gtl::coordinate_traits<coordinate_type>::coordinate_distance distance = gtl::eucliean_distance(*pPattern1, *pPattern2);
					if (distance < m_db.coloring_distance)
					{
						size_t j = mPattern2Idx[pPattern2];
						pair<edge_descriptor, bool> e = edge(i, j, dg);
						if (!e.second) // avoid duplicate edges 
						{
							e = add_edge(i, j, dg);
							assert(e.second);
							put(edge_weight, dg, e.first, 1);
						}
					}
				}
			}
		}
	}

	// graph simplification 
	typedef limbo::algorithms::coloring::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (g);
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
	// keep the order of simplification 
	gs.hide_small_degree(m_db.color_num); // hide vertices with degree smaller than color_num
	if (m_db.color_num == 3)
		gs.merge_subK4(); // merge sub-K4 structure 
	// collect simplified information 
	vector<graph_simplification_type::vertex_status_type> vStatus const& = gs.status();
	vector<vector<vertex_descriptor> > const& vChildren = gs.children();
	stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();
	pair<graph_type, map<vertex_descriptor, vertex_descriptor> > sg = gs.simplified_graph(); // simplified graph and vertex mapping

	// solve coloring 
	typedef limbo::algorithms::coloring::LPColoring coloring_solver_type;
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
	
	// recover colors for simplified vertices with balanced assignment 
	// recover merged vertices 
	for (size_t v = 0; v != vStatus.size(); ++v)
	{
		if (vStatus[v] == graph_simplification_type::GOOD)
		{
			assert(vColor[v] >= 0 && vColor[v] < m_db.color_num);
			for (size_t j = 0; j != vChildren[v].size(); ++j)
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
		vector<coordinate_distance> vDist (m_db.color_num, std::numeric_limits<coordinate_distance>::max());
		for (size_t u = 0; u != pattern_cnt; ++u)
		{
			if (v == u) continue;
			// skip uncolored vertices 
			if (vColor[u] < 0) continue; 
			// we consider eucliean distance
			gtl::coordinate_traits<coordinate_type>::coordinate_distance distance = gtl::eucliean_distance(*(itBgn+v), *(itBgn+u));
			vDist[ vColor[u] ] = std::min(vDist[ vColor[u] ], distance);
		}

		// choose the color with largest distance 
		int8_t best_color = -1;
		coordinate_distance best_dist = std::numeric_limits<coordinate_distance>::min();
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
}

} // namespace SimpleMPL

