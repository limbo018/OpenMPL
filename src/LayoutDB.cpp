/*************************************************************************
    > File Name: LayoutDB.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:45:41 AM CDT
 ************************************************************************/

#include "LayoutDB.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;
namespace bg = boost::geometry;
namespace bgi = bg::index;

LayoutDB::LayoutDB() 
    : LayoutDB::base_type()
{
    initialize();
}
LayoutDB::LayoutDB(LayoutDB::coordinate_type xl, LayoutDB::coordinate_type yl, LayoutDB::coordinate_type xh, LayoutDB::coordinate_type yh) 
    : LayoutDB::base_type(xl, yl, xh, yh)
{
    initialize();
}
LayoutDB::LayoutDB(LayoutDB const& rhs) 
    : LayoutDB::base_type(rhs)
{
    copy(rhs);
}
LayoutDB::~LayoutDB()
{
    // recycle 
    for (std::vector<rectangle_pointer_type>::iterator it = vPatternBbox.begin(), ite = vPatternBbox.end(); it != ite; ++it)
        delete *it;
}
LayoutDB& LayoutDB::operator=(LayoutDB const& rhs)
{
    if (this != &rhs)
    {
        this->base_type::operator=(rhs);
        copy(rhs);
    }
    return *this;
}
void LayoutDB::initialize()
{
    strname                  = "TOPCELL";
    unit                     = 0.001;
    coloring_distance        = 0;
}
void LayoutDB::copy(LayoutDB const& rhs)
{
    hPath    = rhs.hPath;
    strname  = rhs.strname;
    unit     = rhs.unit;
    coloring_distance = rhs.coloring_distance;
    // options
    parms    = rhs.parms;
    // layout information
    tPatternBbox = rhs.tPatternBbox;
    vPatternBbox = rhs.vPatternBbox;
}
void LayoutDB::swap(LayoutDB& rhs)
{
    std::swap(hPath, rhs.hPath);
    std::swap(strname, rhs.strname);
    std::swap(unit, rhs.unit);
    std::swap(coloring_distance, rhs.coloring_distance);
    parms.swap(rhs.parms);
    std::swap(tPatternBbox, rhs.tPatternBbox);
    std::swap(vPatternBbox, rhs.vPatternBbox);
}
void LayoutDB::add(int32_t layer, std::vector<point_type> const& vPoint)
{
    // classify features 
    // sometimes path may be defined as boundary 
    if (parms.sPathLayer.count(layer))
        this->add_path(layer, vPoint);
    else 
        this->add_pattern(layer, vPoint);
}
void LayoutDB::add_path(int32_t layer, std::vector<point_type> const& vPoint)
{
    if (vPoint.size() < 2) return;
    else if (vPoint.size() == 4) // probably it is initialized from boundary 
    {
        gtl::coordinate_traits<coordinate_type>::coordinate_distance dist[] = {
            gtl::euclidean_distance(vPoint[0], vPoint[1]),
            gtl::euclidean_distance(vPoint[1], vPoint[2])
        };
        // simply check if one edge is much longer than the other and choose the longer one  
        int32_t offset = -1;
        if (dist[0] > 10*dist[1])
            offset = 0;
        else if (10*dist[0] < dist[1])
            offset = 1;
        if (offset >= 0)
        {
            path_type p (vPoint[offset], vPoint[offset+1]);
            if (hPath.count(layer))
                hPath[layer].push_back(p);
            else mplAssert(hPath.insert(make_pair(layer, std::vector<path_type>(1, p))).second);
            return;
        }
    }
}
void LayoutDB::initialize_data()
{
    // remember to set coloring_distance from coloring_distance_nm and unit 
    coloring_distance = (coordinate_difference)round(coloring_distance_nm()/(unit*1e+9));
}
void LayoutDB::report_data_kernel() const
{
    mplPrint(kINFO, "Total patterns # = %lu\n", vPatternBbox.size());
    mplPrint(kINFO, "Shape mode = %s\n", std::string(shape_mode()).c_str());
    mplPrint(kINFO, "Coloring distance = %lld db ( %g nm )\n", coloring_distance, coloring_distance_nm());
    mplPrint(kINFO, "Color num = %d\n", color_num());
    mplPrint(kINFO, "Simplification level = %u\n", simplify_level());
    mplPrint(kINFO, "Thread num = %u\n", thread_num());
    mplPrint(kINFO, "Uncolored layer # = %lu", parms.sUncolorLayer.size());
    if (!parms.sUncolorLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = parms.sUncolorLayer.begin(), ite = parms.sUncolorLayer.end(); it != ite; ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Precolored layer # = %lu", parms.sPrecolorLayer.size());
    if (!parms.sPrecolorLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = parms.sPrecolorLayer.begin(), ite = parms.sPrecolorLayer.end(); it != ite; ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Path layer # = %lu", parms.sPathLayer.size());
    if (!parms.sPathLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = parms.sPathLayer.begin(), ite = parms.sPathLayer.end(); it != ite; ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Algorithm = %s\n", std::string(algo()).c_str());
}
void LayoutDB::update_bbox(base_type const& bbox)
{
    if (vPatternBbox.empty())
        gtl::assign(*this, bbox);
    else 
        gtl::encompass(*this, bbox);
}
void LayoutDB::check_layer_and_color(int32_t layer, bool& pattern_layer_flag, int8_t& color) const
{
    pattern_layer_flag = true;
    color = -1;
    if (parms.sPrecolorLayer.count(layer)) 
    {
        // for precolored patterns 
        // set colors 
        color = layer-*parms.sPrecolorLayer.begin();
    }
    else if (parms.sUncolorLayer.count(layer))
    {
        // for uncolored patterns
    }
    else pattern_layer_flag = false;
}
/// remove overlapping patterns 
/// based on scanline approach, horizontal sort and then vertical sort 
/// O(nlogn)
void LayoutDB::remove_overlap(std::vector<LayoutDB::rectangle_pointer_type>& vTargetPattern)
{
    // do nothing if empty
    if (vTargetPattern.empty()) return;
    std::sort(vTargetPattern.begin(), vTargetPattern.end(), compare_rectangle_type());
    // only duplicate is removed so far
    // TO DO: remove overlapping patterns 
    boost::dynamic_bitset<uint32_t, std::allocator<uint32_t> > vValid (vTargetPattern.size()); // use a bit set to record validity 
    vValid.set(); // set all to 1
#ifdef DEBUG
    mplAssert(vValid[0] && vValid[vTargetPattern.size()-1]);
#endif
    uint32_t duplicate_cnt = 0;
    for (std::vector<rectangle_pointer_type>::iterator it1 = vTargetPattern.begin(), it2 = vTargetPattern.begin();
            it2 != vTargetPattern.end(); ++it2)
    {
        if (it2 != vTargetPattern.begin()) 
        {
            rectangle_pointer_type& pPattern1 = *it1;
            rectangle_pointer_type& pPattern2 = *it2;

            if (gtl::equivalence(*pPattern1, *pPattern2)) 
            {
#ifdef DEBUG
                mplPrint(kWARN, "%s duplicates with %s ignored %u\n", std::string(*pPattern1).c_str(), std::string(*pPattern2).c_str(), duplicate_cnt);
#endif
                vValid[pPattern2->pattern_id()] = false;
                duplicate_cnt += 1;
                continue;
            }
        }
        it1 = it2;
    }
    mplPrint(kINFO, "Ignored %u duplicate patterns\n", duplicate_cnt);

    // erase invalid patterns 
    uint32_t invalid_cnt = 0;
    for (uint32_t i = 0; i < vTargetPattern.size(); )
    {
        if (!vValid[vTargetPattern[i]->pattern_id()])
        {
            std::swap(vTargetPattern[i], vTargetPattern.back());
            delete vTargetPattern.back(); // recycle 
            vTargetPattern.pop_back();
            invalid_cnt += 1;
        }
        else ++i;
    }
    mplAssert(duplicate_cnt == invalid_cnt);
#ifdef DEBUG
    for (uint32_t i = 0; i != vTargetPattern.size(); ++i)
        mplAssert(vValid[vTargetPattern[i]->pattern_id()]);
#endif
    // update pattern_id 
    for (uint32_t i = 0; i != vTargetPattern.size(); ++i)
        vTargetPattern[i]->pattern_id(i);
}


SIMPLEMPL_END_NAMESPACE
