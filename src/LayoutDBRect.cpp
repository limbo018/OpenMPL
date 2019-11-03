/*************************************************************************
    > File Name: LayoutDBRect.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:48:22 AM CDT
 ************************************************************************/

#include "LayoutDBRect.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;
namespace bg = boost::geometry;
namespace bgi = bg::index;

//////////////// LayoutDBRect ////////////////
LayoutDBRect::LayoutDBRect() : LayoutDBRect::base_type() 
{
}
LayoutDBRect::LayoutDBRect(LayoutDBRect::coordinate_type xl, LayoutDBRect::coordinate_type yl, LayoutDBRect::coordinate_type xh, LayoutDBRect::coordinate_type yh) 
    : LayoutDBRect::base_type(xl, yl, xh, yh) 
{
}
LayoutDBRect::LayoutDBRect(LayoutDBRect const& rhs) : LayoutDBRect::base_type(rhs)
{
    copy(rhs);
}
LayoutDBRect::~LayoutDBRect()
{
}
LayoutDBRect& LayoutDBRect::operator=(LayoutDBRect const& rhs)
{
    if (this != &rhs)
    {
        this->base_type::operator=(rhs);
        copy(rhs);
    }
    return *this;
}
void LayoutDBRect::copy(LayoutDBRect const&)
{
}
void LayoutDBRect::add_pattern(int32_t layer, std::vector<point_type> const& vPoint)
{
    // collect patterns 
    bool pattern_layer_flag = true;
    int8_t color = -1;
    check_layer_and_color(layer, pattern_layer_flag, color);
    // skip layers not interseted 
    if (!pattern_layer_flag) return;

    mplAssert(vPoint.size() >= 4 && vPoint.size() < 6);

    rectangle_pointer_type pPattern(new rectangle_type());
    for (std::vector<point_type>::const_iterator it = vPoint.begin(); it != vPoint.end(); ++it)
    {
        if (it == vPoint.begin())
            gtl::set_points(*pPattern, *it, *it);
        else 
            gtl::encompass(*pPattern, *it);
    }
    pPattern->layer(layer);
    pPattern->pattern_id(vPatternBbox.size());

    // update layout boundary 
    update_bbox(*pPattern);

    pPattern->color(color);
    // collect pattern 
    // initialize rtree later will contribute to higher efficiency in runtime 
    vPatternBbox.push_back(pPattern);
}

void cal_bound(){
    std::vector<coordinate_type>().swap(boundaries);
	mplAssert(boundaries.empty());
	uint32_t vertex_num = vPatternBbox.size();
	coordinate_type tmp_bound[4] = {INT_MAX,0,INT_MAX ,0};


	for (uint32_t i = 0; i < vertex_num; i++)
    {
		coordinate_type test0 = gtl::xl(*(vPatternBbox[i]));
		coordinate_type test1 = gtl::xh(*(vPatternBbox[i]));
		coordinate_type test2 = gtl::yl(*(vPatternBbox[i]));
		coordinate_type test3 = gtl::yh(*(vPatternBbox[i]));
		tmp_bound[0] = std::min(test0,tmp_bound[0]);
		tmp_bound[1] = std::max(test1,tmp_bound[1]);
		tmp_bound[2] = std::min(test2,tmp_bound[2]);
		tmp_bound[3] = std::max(test3,tmp_bound[3]);
	}
	boundaries.push_back(tmp_bound[0]);
	boundaries.push_back(tmp_bound[1]);
	boundaries.push_back(tmp_bound[2]);
	boundaries.push_back(tmp_bound[3]);
	return;
}
/// call it to initialize rtree 
/// it should be faster than gradually insertion 
void LayoutDBRect::initialize_data()
{
    this->base_type::initialize_data();
    // I assume there may be duplicate in the input gds, but no overlapping patterns 
    // duplicates are removed with following function
    remove_overlap(vPatternBbox);
    // construction with packing algorithm 
    rtree_type tTmp (vPatternBbox.begin(), vPatternBbox.end());
    uint32_t vertex_num = vPatternBbox.size();
    for (uint32_t i = 0; i < vertex_num; i++)
        vRectBeginId.push_back(i);
    tPatternBbox.swap(tTmp);
}
void LayoutDBRect::set_color(uint32_t pattern_id, int8_t color)
{
    rectangle_pointer_type pPattern = vPatternBbox[pattern_id];
    if (pPattern->color() >= 0 && pPattern->color() < color_num()) // check precolored pattern 
        mplAssert(pPattern->color() == color);
    else // assign color to uncolored pattern 
        pPattern->color(color);
}
LayoutDBRect::point_type LayoutDBRect::get_point_closest_to_center(uint32_t pattern_id) const
{
    point_type resPoint; 
    gtl::center(resPoint, *vPatternBbox[pattern_id]); 
    return resPoint;
}
void LayoutDBRect::report_data() const 
{
    mplPrint(kINFO, "Input data for rectangle based layout...\n");
    this->base_type::report_data_kernel();
}

SIMPLEMPL_END_NAMESPACE
