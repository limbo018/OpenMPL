/*************************************************************************
    > File Name: SimpleMPL.h
    > Author: Yibo Lin, Qi Sun
    > Mail: yibolin@utexas.edu, qsun@cse.cuhk.edu.hk
    > Created Time: Wed May 20 22:21:16 2015
 ************************************************************************/

#ifndef SIMPLEMPL_SIMPLEMPL_H
#define SIMPLEMPL_SIMPLEMPL_H

#include <iostream>
#include <stack>
#include <string>
#include <deque>
#include <limbo/algorithms/coloring/Coloring.h>
#include "GdsiiIO.h"
#include <limbo/geometry/api/GeoBoostPolygonApi.h>
#include <limbo/geometry/Geometry.h>
#include <boost/polygon/interval_data.hpp>

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
        void produce_graph_run(int32_t argc, char** argv);
        void read_cmd(int32_t argc, char** argv);
        void read_gds();
        void write_gds();
        /// solve decomposition
        void solve(std::string simplified_graph);
        /// report statistics 
        void report() const;
        /// print welcome information
        void print_welcome() const;

    protected:
        /// initialize graph from layoutdb_type
        void construct_graph();
        /// initialize graph from layoutdb_type
        /// This method is created by Qi Sun to output the adjacency matrix of the graph.
        void construct_graph_with_outputs(std::string vertex_file_name, std::string adjacency_file_name);
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
        uint32_t solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::string simplified_graph);
        /// Created by Qi Sun
        /// Store every component (this function was implemented long long ago, for other uses)
        void store_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);
        /// Created by Qi Sun
        /// Store every component for Dancinglinks
        void store_component_dlx(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id);
        /// kernel for coloring a component
        uint32_t coloring_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::string simplified_graph);
        /// create solver and initialize 
        /// \parm sg is the simplified graph 
        /// \return a point of solver base type
        lac::Coloring<graph_type>* create_coloring_solver(graph_type const& sg) const;
        /// given a graph, solve coloring, contain nested call for itself 
        /// \param dg is decomposition graph before simplification
        uint32_t solve_graph_coloring(uint32_t comp_id, graph_type const& dg, 
                std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
                uint32_t simplify_strategy, std::vector<int8_t>& vColor, std::string simplified_graph) const;
        /// given a component, construct graph, mapping from global index to local index, and set precolor 
        void construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt, 
                graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor) const;

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
        std::vector<uint32_t>               m_vVertexOrder; ///< vertex id
        std::vector<std::vector<uint32_t> > m_mAdjVertex;   ///< adjcency list
        std::vector<uint32_t>               m_vCompId;      ///< independent component id
        uint32_t                            m_comp_cnt;     ///< max# of connected components

        /// density balancing 
        std::vector<uint32_t> m_vColorDensity; ///< number of colors used so far 

        /// conflict report 
        mutable std::vector<std::pair<uint32_t, uint32_t> > m_vConflict; ///< conflict patterns  

        //*********************** Stitch Insertion ***********************//
    public:
        // void relation4rect();
        // I think the storage operation is a little complex, so I combine the storage operation and projection() into one.
        // The vBookmark can help index the starting position of each component.
        void runProjection(const std::vector<uint32_t> & vBookmark);
        // Used to generate projections
        // The reference variables are used to store the intermediate results.
        void projection(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd,
            uint32_t comp_id, std::vector<uint32_t> & new_CompId, std::vector<uint32_t> & new_VertexOrder, 
            std::vector<std::vector<uint32_t>> & new_AdjVertex, std:vector<rectangle_pointer_type> new_PatternBox);
        // Constructor, I am not sure whether the pitch is useful here.
        // But the default value is 0.
        SimpleMPL(double pitch);

    protected:
        // get the width of the rectangle
        coordinate_type getWidth(rectangle_pointer_type rect) { return gtl::xh(*rect) - gtl::xl(*rect); }
        // get the height of the rectangle
        coordinate_type getHeight(rectangle_pointer_type rect) { return gtl::yh(*rect) - gtl::yl(*rect); }


        // judge whether the rectangle is horizontal, which may affect the following stitch insertion strategy.
        bool whetherHorizontal (rectangle_pointer_type tmp);

        // Generate Stitch Insertion Points for DPL
        void StitchGenerateDPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, const coordinate_type left, const coordinate_type right, bool vddstitch = false);
	    // Generate Stitch Insertion Points for TPL
        void StitchGenerateTPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, const coordinate_type left, const coordinate_type right, bool vddstitch = false);

        std::vector<rectangle_type> m_final_rectangle;

        std::vector<std::vector<uint32_t>> m_Touch;     // This vector is for polygon type data, storing the touching rectangles.
        
        std::vector<std::vector<uint32_t>> m_Safe;      // This vector stores the safe patterns.
        std::vector<std::vector<rectangle_type>> m_interRect;

        std::vector<std::vector<uint32_t>>  m_StitchNeigh;  // This vector stores the patterns which belong to the same parent pattern. 

        double MINDPLDIST           = 0;    // minimual DPL distance
        double PITCH                = 0;    // still not sure whether use pitch in this problem
        bool dplstitch              = true; // whether use DPL method to generate stitch candidates, if false, then use TPL method.

        coordinate_type m_layout_left       = INT32_MAX;
        coordinate_type m_layout_right      = INT32_MIN;
        coordinate_type m_layout_top        = INT32_MIN;
        coordinate_type m_layout_down       = INT32_MAX;

};

SIMPLEMPL_END_NAMESPACE

#endif