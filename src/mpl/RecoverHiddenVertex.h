/*************************************************************************
    > File Name: RecoverHiddenVertex.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 08 Oct 2015 01:19:55 PM CDT
 ************************************************************************/

#ifndef SIMPLEMPL_RECOVERHIDDENVERTEX_H
#define SIMPLEMPL_RECOVERHIDDENVERTEX_H

#include <stack>
#include "LayoutDB.h"

SIMPLEMPL_BEGIN_NAMESPACE

/// ==========================================================
/// Function objects to recover colors of vertices that are 
/// hidden by graph simplification HIDE_SMALL_DEGREE process. 
/// The recovery can subject to various cost functions  
/// ==========================================================

/// base class for vertex recovery 
/// simply recover with first available color 
class RecoverHiddenVertex
{
    public:
        typedef LayoutDB layoutdb_type;
        typedef layoutdb_type::coordinate_type coordinate_type;
        typedef layoutdb_type::coordinate_difference   coordinate_difference;
        typedef layoutdb_type::point_type              point_type;
        typedef layoutdb_type::rectangle_type          rectangle_type;
        typedef layoutdb_type::rectangle_pointer_type  rectangle_pointer_type;
        typedef layoutdb_type::graph_type              graph_type;
        typedef layoutdb_type::vertex_descriptor       vertex_descriptor;
        typedef layoutdb_type::edge_descriptor         edge_descriptor;
        typedef layoutdb_type::graph_edge_type         graph_edge_type;

        /// constructor 
        RecoverHiddenVertex(graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
                std::vector<int8_t>& vColor, std::stack<vertex_descriptor>& vHiddenVertices, 
                std::vector<uint32_t> const& vColorDensity, layoutdb_type const& db);
        /// destructor
        virtual ~RecoverHiddenVertex() {}

        /// top api 
        virtual void operator()();
    protected:
        /// recover a single vertex 
        virtual void recover_vertex(vertex_descriptor v);
        /// find available colors for a single vertex 
        /// this function does not to be virtual function, because it always remains the same 
        void find_unused_colors(vertex_descriptor v);
        /// find best color for a single vertex from available colors 
        /// \return best color 
        virtual int8_t find_best_color(vertex_descriptor v);

        graph_type const& m_dg;
        std::vector<uint32_t>::const_iterator m_itBgn;
        uint32_t m_pattern_cnt;
        std::vector<int8_t>& m_vColor;
        std::stack<vertex_descriptor>& m_vHiddenVertices;
        std::vector<uint32_t> const& m_vColorDensity;
        layoutdb_type const& m_db;
        std::vector<char> m_vUnusedColor; ///< local variable to avoid frequent construction
        std::vector<char> m_vStitchColor; ///< local variable to store the stitch neighbors' colors
};

/// recover with distance heuristic 
/// try to maximize minimum distance 
class RecoverHiddenVertexDistance : public RecoverHiddenVertex
{
    public:
        typedef RecoverHiddenVertex base_type;

        /// constructor
        RecoverHiddenVertexDistance(graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
                std::vector<int8_t>& vColor, std::stack<vertex_descriptor>& vHiddenVertices, 
                std::vector<uint32_t> const& vColorDensity, layoutdb_type const& db); 

    protected:
        virtual int8_t find_best_color(vertex_descriptor v);

        std::vector<coordinate_difference> m_vDist; ///< minimum distance for different colors to current vertex 
};

/// recover with image contrast heuristic 
/// try to maximize image contrast score
class RecoverHiddenVertexImageContrast : public RecoverHiddenVertex
{
    public:
        typedef RecoverHiddenVertex base_type;

        /// constructor
        RecoverHiddenVertexImageContrast(graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
                std::vector<int8_t>& vColor, std::stack<vertex_descriptor>& vHiddenVertices, 
                std::vector<uint32_t> const& vColorDensity, layoutdb_type const& db); 

    protected:
        virtual int8_t find_best_color(vertex_descriptor v);
};

SIMPLEMPL_END_NAMESPACE

#endif
