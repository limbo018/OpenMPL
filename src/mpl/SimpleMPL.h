/*************************************************************************
    > File Name: SimpleMPL.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:21:16 2015
 ************************************************************************/

#ifndef SIMPLEMPL_SIMPLEMPL_H
#define SIMPLEMPL_SIMPLEMPL_H
#define LIWEI_BEFORE 0
#define LIWEI_DEBUG 1
#define _DGOUT
#include <iostream>
#include <stack>
#include <limbo/algorithms/coloring/Coloring.h>
#include "GdsiiIO.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace la = limbo::algorithms;
namespace lac = la::coloring;

class SimpleMPL
{
	public:
		typedef LayoutDB layoutdb_type;
        typedef layoutdb_type::coordinate_type coordinate_type;
		typedef layoutdb_type::coordinate_difference   coordinate_difference;
		typedef layoutdb_type::point_type              point_type;
		typedef layoutdb_type::rectangle_type          rectangle_type;
		typedef layoutdb_type::polygon_type            polygon_type;
		typedef layoutdb_type::polygon_pointer_type    polygon_pointer_type;
		typedef layoutdb_type::rectangle_pointer_type  rectangle_pointer_type;
		typedef layoutdb_type::path_type               path_type;
		typedef layoutdb_type::rtree_type              rtree_type;
        typedef layoutdb_type::graph_type              graph_type;
        typedef layoutdb_type::vertex_descriptor       vertex_descriptor;
        typedef layoutdb_type::edge_descriptor         edge_descriptor;

        /// default constructor 
        SimpleMPL();
        /// destructor 
        ~SimpleMPL();

		/// top api to solve decomposition
		void run(int32_t argc, char** argv);
		void read_cmd(int32_t argc, char** argv);
		void read_gds();
		void write_gds();
		/// solve decomposition
		void solve();
		/// report statistics 
		void report() const;
        /// print welcome information
        void print_welcome() const;
	protected:
		void lgSimplification();
		
		void lgSimplification(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);
		
		void dgSimplColoring();

		void dgSimplColoring(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::map<vertex_descriptor, std::set<uint32_t> >& m_ArtiPoint, std::vector<std::vector<vertex_descriptor> >& m_CompVertex);
		
		void projection();

		uint32_t merge_K4_coloring(graph_type const& dg, std::vector<int8_t>& vColor) const;

		void splitRectangle(rectangle_type& pRect, std::vector<rectangle_pointer_type>& split, std::vector<rectangle_pointer_type> nei_Vec);
		
		void GenerateStitchPosition(const rectangle_type pRect, std::vector<rectangle_type> vInterSect, std::vector<coordinate_type>& vPossibleStitches, std::vector<coordinate_type>& vstitches);
		
		void setVddGnd();
		void cal_boundaries();
		void reconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t>& new_polygon_id_list, std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list);
		
		void updateConflictRelation();
		
		/// merge vdd, and conduct biconnected component devision once more, and simplify the graph further.
		void vdd_biconnected_component();
		/// initialize graph from layoutdb_type
		void construct_graph();
        /// construct graph from coloring distance, set adjacency list m_mAdjVertex 
        /// pass \param vertex_num
        /// \return total number of edges 
        uint32_t construct_graph_from_distance(uint32_t vertex_num);
        /// construct graph from paths in gdsii file, set adjacency list m_mAdjVertex 
        /// pass \param vertex_num
        /// \return total number of edges 
        uint32_t construct_graph_from_paths(uint32_t vertex_num);
		/// compute connected component 
		void connected_component();
		/// DFS for connected component computation
		void depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& pattern_id);
		/// solve a single component 
        /// it wraps up solve_graph_coloring()
		uint32_t solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);
        /// kernel for coloring a component 
        uint32_t coloring_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);
        /// create solver and initialize 
        /// \parm sg is the simplified graph 
        /// \return a point of solver base type
        lac::Coloring<graph_type>* create_coloring_solver(graph_type const& sg) const;
        /// given a graph, solve coloring, contain nested call for itself 
        /// \param dg is decomposition graph before simplification
        uint32_t solve_graph_coloring(uint32_t comp_id, graph_type const& dg, 
                std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
                uint32_t simplify_strategy, std::vector<int8_t>& vColor, std::set<vertex_descriptor> vdd_set) ;
        /// given a component, construct graph, mapping from global index to local index, and set precolor 
        void construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt, 
                graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor, std::set<vertex_descriptor>& vdd_set, bool flag) const;

		/// report conflict number for a component 
		uint32_t conflict_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd) const;
		/// report conflict number for the whole layout 
		/// collect conflict patterns to m_vConflict
		uint32_t conflict_num() const;
        /// reset data members 
        /// \param init denote whether run in initialize mode 
        void reset(bool init);
        /// check whether a component contains non-colored patterns 
        /// \return false if all precolored 
        bool check_uncolored(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd) const;

        /// for debug 
        /// \param g is mutable because edge properties for boost::dynamic_properties need mutable graph 
        /// \param filename should not contain extension 
        void write_graph(graph_type& g, std::string const& filename) const;

        layoutdb_type* m_db; ///< pointer of layout database and user-defined options 
		/// adjacency list data structure for a graph 
		std::vector<uint32_t>					m_vVertexOrder;		///< vertex id, vertices in the same component are abutting,
		// In vPatternBbox, the vertices are in the order of vertex id, but in m_vVertexOrder, verteices are in the order of neighboring relationships.
		std::vector<std::vector<uint32_t> >		m_mAdjVertex;		///< adjcency list
		std::vector<uint32_t>					m_vCompId;			///< independent component id
		uint32_t								m_comp_cnt;			///< max# of connected components
		std::vector<uint32_t>					vBookmark;			///< a bookmark marking the start and end locations of different components
		/// density balancing 
		std::vector<uint32_t>					m_vColorDensity;	///< number of colors used so far 
		// std::vector<uint32_t>					IsVDDGND;


		/// conflict report 
		mutable std::vector<std::pair<uint32_t, uint32_t> > m_vConflict; ///< conflict patterns  

		///< in_DG : controls whether that node has been removed in lgSimplification, in other words, that's stitch candidate.
		///< also, only when a node is in DG, it will generate projections on other nodes.
		std::vector<bool>						in_DG;		///< store the nodes left after lgSimplification, which means these nodes should be in DGs.
		///< articulation_set : if it's articulation point, we won't split it.
		std::vector<bool>						articulation_vec; ///< store the global articulation polygons
		///< if it's VGGGND, we won't split it.
		std::vector<bool>						isVDDGND;			  ///< whether this node is VDDGND
		std::vector<uint32_t>					new2ori_polygon;
		std::vector<std::vector<uint32_t> >		ori2new_polygon;
		std::vector<std::vector<uint32_t> >		StitchRelation;

		std::vector<uint32_t>					dgCompId;
		std::map<uint32_t, uint32_t>			dgGlobal2Local;
		uint32_t								globalCompId;
		std::vector<std::vector<uint32_t> >		vdd_multi_comp;
};

SIMPLEMPL_END_NAMESPACE

#endif
