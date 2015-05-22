/*************************************************************************
    > File Name: SimpleMPL.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:21:16 2015
 ************************************************************************/

#ifndef SIMPLEMPL_H
#define SIMPLEMPL_H

#include <iostream>
#include "io.h"

namespace SimpleMPL {

class SimpleMPL
{
	public:
		typedef int32_t coordinate_type;
		typedef LayoutDB<coordinate_type> layoutdb_type;
		typedef typename layoutdb_type::point_type point_type;
		typedef typename layoutdb_type::rectangle_type rectangle_type;
		typedef typename layoutdb_type::polygon_type polygon_type;
		typedef typename layoutdb_type::polygon_pointer_type polygon_pointer_type;
		typedef typename layoutdb_type::rectangle_pointer_type rectangle_pointer_type;
		typedef typename layoutdb_type::path_type path_type;

#if 0
		// do not use setS, it does not compile for subgraph
		// do not use custom property tags, it does not compile for most utilities
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
				boost::property<boost::vertex_index_t, std::size_t, boost::property<boost::vertex_color_t, int> >, 
				boost::property<boost::edge_index_t, std::size_t, boost::property<boost::edge_weight_t, int> >,
				boost::property<boost::graph_name_t, string> > graph_type;
#endif
		/// top api to solve decomposition
		void run(int32_t argc, char** argv);
		void read_cmd(int32_t argc, char** argv);
		void read_gds();
		void write_gds();
		/// solve decomposition
		void solve();
	protected:
		/// compute connected component 
		void connected_component();
		/// DFS for connected component computation
		void depth_first_search(rectangle_pointer_type source, uint32_t comp_id, uint32_t& pattern_id);
		/// solve a single component 
		void solve_component(const vector<rectangle_pointer_type>::const_iterator itBgn, const vector<rectangle_pointer_type>::const_iterator itEnd);

		layoutdb_type m_db; ///< layout database and user-defined options 
		uint32_t m_comp_cnt; ///< maximum number of connected components 
};

} // namespace SimpleMPL

#endif
