/*************************************************************************
    > File Name: RecoverHiddenVertex.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 08 Oct 2015 01:35:35 PM CDT
 ************************************************************************/

#include "RecoverHiddenVertex.h"

SIMPLEMPL_BEGIN_NAMESPACE

RecoverHiddenVertex::RecoverHiddenVertex(RecoverHiddenVertex::graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        std::vector<int8_t>& vColor, std::stack<RecoverHiddenVertex::vertex_descriptor>& vHiddenVertices, 
        std::vector<uint32_t> const& vColorDensity, RecoverHiddenVertex::layoutdb_type const& db)
    : m_dg (dg)
    , m_itBgn (itBgn)
    , m_pattern_cnt (pattern_cnt)
    , m_vColor (vColor)
    , m_vHiddenVertices (vHiddenVertices)
    , m_vColorDensity (vColorDensity)
    , m_db (db)
    , m_vStitchColor(db.color_num())
    , m_vUnusedColor(db.color_num())
{
}

void RecoverHiddenVertex::operator()()
{
    std::cout << "in RecoverHiddenVerte.\n";
	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
	while (!m_vHiddenVertices.empty())
	{
		vertex_descriptor v = m_vHiddenVertices.top();
		m_vHiddenVertices.pop();

        // recover color for v 
        recover_vertex(v);
	}
}

void RecoverHiddenVertex::recover_vertex(RecoverHiddenVertex::vertex_descriptor v)
{
    rectangle_pointer_type temp = m_db.vPatternBbox[*(m_itBgn+v)];
    std::cout << "\nid " << temp->pattern_id() << " : " << std::endl;
    // find available colors 
    find_unused_colors(v);
    // find best color 
    int8_t best_color = find_best_color(v);
    // assign color 
    mplAssert(best_color >= 0 && best_color < m_db.color_num());

    std::cout << "location : " << gtl::xl(*temp) << ", " << gtl::yl(*temp) 
            << ", " << gtl::xh(*temp) << ", " << gtl::yh(*temp) << ", color " << +unsigned(best_color) << std::endl;
    m_vColor[v] = best_color;
}

void RecoverHiddenVertex::find_unused_colors(RecoverHiddenVertex::vertex_descriptor v)
{
     // find available colors 
    std::fill(m_vStitchColor.begin(), m_vStitchColor.end(), false);
    std::fill(m_vUnusedColor.begin(), m_vUnusedColor.end(), true);
    boost::graph_traits<graph_type>::adjacency_iterator vi, vie;
    int count = 0;
    for (boost::tie(vi, vie) = adjacent_vertices(v, m_dg); vi != vie; ++vi)
    {
        count ++;
        vertex_descriptor u = *vi;
        if (m_vColor[u] >= 0)
        {
            mplAssert(m_vColor[u] < m_db.color_num());
            // added by Qi Sun to sovle the stitch color
            std::pair<graph_edge_type, bool> e12 = boost::edge(v, u, m_dg);
            assert(e12.second);
            if (boost::get(boost::edge_weight, m_dg, e12.first) < 0) 
            {
                std::cout << "stitch : " << +unsigned(m_vColor[u]) << std::endl;
                m_vStitchColor[m_vColor[u]] = true;
            }
            else 
            {
                std::cout << "conflict : " << +unsigned(m_vColor[u]) << std::endl; 
                m_vUnusedColor[m_vColor[u]] = false;
            }
        }
    }
    std::cout << "deg : " << count << std::endl;
}

int8_t RecoverHiddenVertex::find_best_color(RecoverHiddenVertex::vertex_descriptor /*v*/) 
{
    // find the first available color 
    int8_t best_color = -1;
    for (int8_t i = 0; i != m_db.color_num(); ++i)
    {
        if (m_vUnusedColor[i] && m_vStitchColor[i])
        {
            best_color = i;
            return best_color;
        }
    }

    for(int8_t i = 0; i != m_db.color_num(); ++i)
    {
        if (m_vUnusedColor[i])
        {
            best_color = i;
            break;
        }
    }
    return best_color;
}

RecoverHiddenVertexDistance::RecoverHiddenVertexDistance(RecoverHiddenVertexDistance::graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        std::vector<int8_t>& vColor, std::stack<RecoverHiddenVertexDistance::vertex_descriptor>& vHiddenVertices, 
        std::vector<uint32_t> const& vColorDensity, RecoverHiddenVertex::layoutdb_type const& db) 
    : RecoverHiddenVertexDistance::base_type(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, vColorDensity, db) 
    , m_vDist(db.color_num())
{
}

int8_t RecoverHiddenVertexDistance::find_best_color(RecoverHiddenVertexDistance::vertex_descriptor v) 
{
    // find the nearest distance of each color 
    // search all patterns in the component 
    // TO DO: further speedup is possible to search a local window 
    std::fill(m_vDist.begin(), m_vDist.end(), std::numeric_limits<coordinate_difference>::max());
    for (uint32_t u = 0; u != m_pattern_cnt; ++u)
    {
        if (v == u) continue;
        // skip uncolored vertices 
        if (m_vColor[u] < 0) continue; 
        // we consider euclidean distance
        // use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
        gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = m_db.euclidean_distance(*m_db.vPatternBbox[*(m_itBgn+v)], *m_db.vPatternBbox[*(m_itBgn+u)]);
#ifdef DEBUG
        mplAssert(m_vColor[u] < m_db.color_num() && distance >= 0);
#endif
        m_vDist[ m_vColor[u] ] = std::min(m_vDist[ m_vColor[u] ], distance);
    }

    // choose the color with largest distance 
    int8_t best_color = -1;
    double best_score = -std::numeric_limits<double>::max(); // negative max 
    for (int8_t i = 0; i != m_db.color_num(); ++i)
    {
        if (m_vUnusedColor[i])
        {
            double cur_score = (double)m_vDist[i]/(1.0+m_vColorDensity[i]);
            if (best_score < cur_score)
            {
                best_color = i;
                best_score = cur_score;
            }
        }
    }
    return best_color;
}

RecoverHiddenVertexImageContrast::RecoverHiddenVertexImageContrast(RecoverHiddenVertexImageContrast::graph_type const& dg, std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        std::vector<int8_t>& vColor, std::stack<RecoverHiddenVertexImageContrast::vertex_descriptor>& vHiddenVertices, 
        std::vector<uint32_t> const& vColorDensity, RecoverHiddenVertex::layoutdb_type const& db) 
    : RecoverHiddenVertexImageContrast::base_type(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, vColorDensity, db) 
{
}

int8_t RecoverHiddenVertexImageContrast::find_best_color(RecoverHiddenVertexImageContrast::vertex_descriptor /*v*/)
{
    mplAssertMsg(0, "not ready yet"); 

    int8_t best_color = -1;
    double best_score = -std::numeric_limits<double>::max();
    for (int8_t i = 0; i != m_db.color_num(); ++i)
    {
        if (m_vUnusedColor[i])
        {
            // compute image contrast score here 
            double cur_score = 0;
            // update best_score and best_color
            if (cur_score > best_score)
            {
                best_color = i;
                best_score = cur_score;
            }
        }
    }
    return best_color;
}

SIMPLEMPL_END_NAMESPACE
