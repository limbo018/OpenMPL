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
		typedef layoutdb_type::coordinate_difference   coordinate_difference;
		typedef layoutdb_type::point_type              point_type;
		typedef layoutdb_type::rectangle_type          rectangle_type;
		typedef layoutdb_type::polygon_type            polygon_type;
		typedef layoutdb_type::polygon_pointer_type    polygon_pointer_type;
		typedef layoutdb_type::rectangle_pointer_type  rectangle_pointer_type;
		typedef layoutdb_type::path_type               path_type;
		typedef layoutdb_type::rtree_type              rtree_type;

		/// top api to solve decomposition
		void run(int32_t argc, char** argv);
		void read_cmd(int32_t argc, char** argv);
		void read_gds();
		void write_gds();
		/// solve decomposition
		void solve();
		/// report statistics 
		void report() const;
	protected:
		/// initialize graph from layoutdb_type
		void construct_graph();
		/// compute connected component 
		void connected_component();
		/// DFS for connected component computation
		void depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& pattern_id);
		/// solve a single component 
		uint32_t solve_component(const vector<uint32_t>::const_iterator itBgn, const vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);

		/// report conflict number for a component 
		uint32_t conflict_num(const vector<uint32_t>::const_iterator itBgn, const vector<uint32_t>::const_iterator itEnd) const;
		/// report conflict number for the whole layout 
		/// collect conflict patterns to m_vConflict
		uint32_t conflict_num() const;

		layoutdb_type m_db; ///< layout database and user-defined options 
		/// adjacency list data structure for a graph 
		vector<uint32_t> m_vVertexOrder; ///< vertex id 
		vector<vector<uint32_t> > m_mAdjVertex; ///< adjcency list 
		vector<uint32_t> m_vCompId; ///< independent component id 
		uint32_t m_comp_cnt; ///< maximum number of connected components 

		/// density balancing 
		vector<uint32_t> m_vColorDensity; ///< number of colors used so far 

		/// conflict report 
		mutable vector<pair<uint32_t, uint32_t> > m_vConflict; ///< conflict patterns  
};

} // namespace SimpleMPL

#endif
