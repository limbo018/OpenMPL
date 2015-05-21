/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL"
#include <stack>
#include <limbo/algorithms/coloring/GraphSimplification.h>

namespace SimpleMPL {

using std::cout;
using std::endl;
using std::stack;

void SimpleMPL::connected_component()
{
	uint32_t comp_id = 0;
	for (vector<rectangle_pointer_type>::iterator it = m_db.vPattern.begin();
			it != m_db.vPattern.end(); ++it)
	{
		rectangle_pointer_type pPattern = *it;

		if (pPattern->comp_id() != std::numeric_limits<uint32_t>::max()) // not visited 
		{
			depth_first_search(pPattern, comp_id);
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
}

void SimpleMPL::depth_first_search(SimpleMPL::rectangle_pointer_type source, uint32_t comp_id)
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
	}
}

void SimpleMPL::solve_component(vector<SimpleMPL::rectangle_pointer_type> const& vPattern)
{
#ifdef DEBUG
	for (vector<rectangle_pointer_type>::const_iterator it = vPattern.begin()+1; it != vPattern.end(); ++it)
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

	// decomposition graph 
	graph_type dg (vPattern.size());

	if (vPattern.size() < 30) // when there are very few patterns, construct graph by direct searching 
	{
		for (size_t i = 0; i != vPattern.size(); ++i)
			for (size_t j = i+1; j != vPattern.size(); ++j)
			{
				// we consider eucliean distance
				gtl::coordinate_traits<coordinate_type>::coordinate_distance distance = gtl::eucliean_distance(*vPattern[i], *vPattern[j]);
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
		rtree_type tPattern (vPattern.begin(), vPattern.end()); 
		map<rectangle_pointer_type, size_t> mPattern2Idx; // map pattern to index in vPattern
		for (size_t i = 0; i != vPattern.size(); ++i)
			assert(mPattern2Idx.insert(std::make_pair(vPattern[i], i)).second);

		for (size_t i = 0; i != vPattern.size(); ++i)
		{
			rectangle_pointer_type pPattern1 = *vPattern[i];

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
	limbo::algorithms::coloring::GraphSimplification<graph_type> gs (g);
	gs.hide_small_degree(m_db.color_num);
	gs.merge_subK4();
	// solve coloring 
	// recover colors for simplified vertices with balanced assignment 
}

} // namespace SimpleMPL

