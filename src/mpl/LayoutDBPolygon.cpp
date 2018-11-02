/*************************************************************************
    > File Name: LayoutDBPolygon.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:50:17 AM CDT
 ************************************************************************/

#include "LayoutDBPolygon.h"
#include <stack>

SIMPLEMPL_BEGIN_NAMESPACE

//////////////// LayoutDBPolygon ////////////////
LayoutDBPolygon::LayoutDBPolygon() : LayoutDBPolygon::base_type() 
{
}
LayoutDBPolygon::LayoutDBPolygon(LayoutDBPolygon::coordinate_type xl, LayoutDBPolygon::coordinate_type yl, LayoutDBPolygon::coordinate_type xh, LayoutDBPolygon::coordinate_type yh) 
    : LayoutDBPolygon::base_type(xl, yl, xh, yh) 
{
}
LayoutDBPolygon::LayoutDBPolygon(LayoutDBPolygon const& rhs) : LayoutDBPolygon::base_type(rhs)
{
    copy(rhs);
}
LayoutDBPolygon::~LayoutDBPolygon()
{
    // recycle 
    for (std::vector<rectangle_pointer_type>::iterator it = vPolyRectPattern.begin(), ite = vPolyRectPattern.end(); it != ite; ++it)
        delete *it;
}
LayoutDBPolygon& LayoutDBPolygon::operator=(LayoutDBPolygon const& rhs)
{
    if (this != &rhs)
    {
        this->base_type::operator=(rhs);
        copy(rhs);
    }
    return *this;
}
void LayoutDBPolygon::copy(LayoutDBPolygon const& rhs)
{
    vParentPolygonId = rhs.vParentPolygonId;
    vPolyRectPattern = rhs.vPolyRectPattern;
    vPolyRectBeginId = rhs.vPolyRectBeginId;
}
void LayoutDBPolygon::add_pattern(int32_t layer, std::vector<point_type> const& vPoint)
{
    // collect patterns 
    bool pattern_layer_flag = true;
    int8_t color = -1;
    check_layer_and_color(layer, pattern_layer_flag, color);
    // skip other layers 
    if (!pattern_layer_flag) return;

    mplAssert(vPoint.size() >= 4);

    polygon_set_type tmpPolygonSet (gtl::HORIZONTAL); 
    polygon_90_data<coordinate_type> tmpPolygon; 
    tmpPolygon.set(vPoint.begin(), vPoint.end());
    tmpPolygonSet.insert(tmpPolygon);

    std::vector<rectangle_data<coordinate_type> > vTmpRectangle;
    tmpPolygonSet.get_rectangles(vTmpRectangle);

    // update layout boundary in LayoutDBPolygon::initialize_data
    rectangle_data<coordinate_type> bbox;
    if (gtl::extents(bbox, tmpPolygonSet))
        update_bbox(bbox);

    // collect decomposed rectangles 
    for (std::vector<rectangle_data<coordinate_type> >::const_iterator it = vTmpRectangle.begin(), ite = vTmpRectangle.end(); it != ite; ++it)
    {
        rectangle_data<coordinate_type> const& rect = *it;
        rectangle_pointer_type pPattern = new rectangle_type(gtl::xl(rect), gtl::yl(rect), gtl::xh(rect), gtl::yh(rect));
        pPattern->layer(layer);
        pPattern->color(color);
        pPattern->pattern_id(vPolyRectPattern.size());

        // collect pattern 
        // initialize rtree later will contribute to higher efficiency in runtime 
        vPolyRectPattern.push_back(pPattern);
    }
}
/// call it to initialize rtree 
/// it should be faster than gradually insertion 
void LayoutDBPolygon::initialize_data()
{
    this->base_type::initialize_data();
    // I assume there may be duplicate in the input gds, but no overlapping patterns 
    // duplicates are removed with following function
    remove_overlap(vPolyRectPattern);
    // tPatternBbox is first used to store vPolyRectPattern
    // construction with packing algorithm 
    rtree_type tTmp1 (vPolyRectPattern.begin(), vPolyRectPattern.end());
    tPatternBbox.swap(tTmp1);
    // compute parent polygons for each rectangle by connected component algorithm 
    compute_parent_polygons();

    tPatternBbox.clear(); // I want to reduce peak memory, but not sure if it works or not
    // after initialization, tPatternBbox is used to store vPatternBbox
    // construction with packing algorithm 
    rtree_type tTmp2 (vPatternBbox.begin(), vPatternBbox.end());
    tPatternBbox.swap(tTmp2);
}
void LayoutDBPolygon::set_color(uint32_t pattern_id, int8_t color)
{
    rectangle_pointer_type pPattern = vPatternBbox[pattern_id];
#ifdef QDEBUG
//	std::cout << "set_color : " <<  pPattern->pattern_id()  << " has color " << +unsigned(pPattern->color()) << " and new color is " << +unsigned(color)<< std::endl;
#endif
	if (pPattern->color() >= 0 && pPattern->color() < color_num()) // check precolored pattern 
        mplAssert(pPattern->color() == color);
    else // assign color to uncolored pattern 
        pPattern->color(color);
    for (uint32_t i = vPolyRectBeginId[pattern_id], ie = vParentPolygonId.size(); i != ie && vParentPolygonId[i] == pattern_id; ++i)
        vPolyRectPattern[i]->color(color);
}
LayoutDBPolygon::point_type LayoutDBPolygon::get_point_closest_to_center(uint32_t pattern_id) const
{
    // given a pattern 
    // check every center of decomposed rectangles 
    // find the one closest to the center of bounding box of the pattern 
    const rectangle_pointer_type pPattern = vPatternBbox[pattern_id];
    point_type bboxCenter;
    gtl::center(bboxCenter, *pPattern);
    point_type resPoint; 
    coordinate_difference bestDist = std::numeric_limits<coordinate_difference>::max();
    uint32_t polyRectId1e = (pattern_id+1 == vPolyRectBeginId.size())? vPolyRectPattern.size() : vPolyRectBeginId[pattern_id+1];
    for (uint32_t polyRectId1 = vPolyRectBeginId[pattern_id]; polyRectId1 != polyRectId1e; ++polyRectId1)
    {
        point_type polyRectCenter;
        gtl::center(polyRectCenter, *vPolyRectPattern[polyRectId1]);
        coordinate_difference distance = gtl::euclidean_distance(bboxCenter, polyRectCenter);
        if (bestDist > distance)
        {
            resPoint = polyRectCenter;
            bestDist = distance;
        }
    }
    return resPoint;
}
void LayoutDBPolygon::compute_parent_polygons()
{
    vParentPolygonId.assign(vPolyRectPattern.size(), std::numeric_limits<uint32_t>::max());
    uint32_t polygon_id = 0;
    uint32_t order_id = 0; // order counter 
    std::vector<uint32_t> vOrderId (vPolyRectPattern.size(), 0); // necessary to have O(N) sorting during DFS
    for (uint32_t v = 0, ve = vPolyRectPattern.size(); v != ve; ++v)
    {
        if (vParentPolygonId[v] == std::numeric_limits<uint32_t>::max()) // not visited 
            depth_first_search(v, polygon_id++, order_id, vOrderId);
    }
    uint32_t num_polygons = polygon_id;

    // reorder vPolyRectPattern so that connected patterns are abutting 
    std::vector<rectangle_pointer_type> vTmpPolyRectPattern (vPolyRectPattern.size());
    for (uint32_t i = 0, ie = vPolyRectPattern.size(); i != ie; ++i)
        vTmpPolyRectPattern[vOrderId[i]] = vPolyRectPattern[i];
    vTmpPolyRectPattern.swap(vPolyRectPattern); 

    std::vector<uint32_t> vTmpParentPolygonId (vPolyRectPattern.size());
    for (uint32_t i = 0, ie = vPolyRectPattern.size(); i != ie; ++i)
        vTmpParentPolygonId[vOrderId[i]] = vParentPolygonId[i];
    vTmpParentPolygonId.swap(vParentPolygonId);

    // set pattern id 
    for (uint32_t i = 0, ie = vPolyRectPattern.size(); i != ie; ++i)
        vPolyRectPattern[i]->pattern_id(i);

#ifdef DEBUG
    mplAssert(vParentPolygonId.back()+1 == num_polygons);
#endif
    // compute vPatternBbox from vPolyRectPattern
    compute_parent_polygon_bboxes(num_polygons);
}
void LayoutDBPolygon::depth_first_search(uint32_t source, uint32_t polygon_id, uint32_t& order_id, std::vector<uint32_t>& vOrderId)
{
    std::stack<uint32_t> vStack; 
	vStack.push(source);

	while (!vStack.empty())
	{
		uint32_t current = vStack.top();
		vStack.pop();
		if (vParentPolygonId[current] == std::numeric_limits<uint32_t>::max()) // not visited 
		{
			vParentPolygonId[current] = polygon_id; // set visited 
            vOrderId[current] = order_id++; // set order 
            const rectangle_pointer_type pPattern = vPolyRectPattern[current];

			// find patterns connected with pPattern 
			// query tPatternBbox 
            std::vector<rectangle_pointer_type> vAdjPattern;
			rectangle_type rect (*pPattern);
			// slightly bloat pPattern 
			gtl::bloat(rect, gtl::HORIZONTAL, 1);
			gtl::bloat(rect, gtl::VERTICAL, 1);
			for (rtree_type::const_query_iterator itq = tPatternBbox.qbegin(bgi::intersects(rect));
					itq != tPatternBbox.qend(); ++itq)
			{
				rectangle_pointer_type const& pAdjPattern = *itq;
				if (pAdjPattern != pPattern) // skip pPattern itself 
				{
					mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
					// we consider euclidean distance
					gtl::coordinate_traits<coordinate_type>::coordinate_difference distance = gtl::euclidean_distance(*pAdjPattern, *pPattern);
					if (distance <= 0) // consider as neighbor
                        vStack.push(pAdjPattern->pattern_id());
				}
			}
		}
	}
}
void LayoutDBPolygon::compute_parent_polygon_bboxes(uint32_t num_polygons)
{
    // vPatternBbox is the bounding box of a polygon 
    vPatternBbox.assign(num_polygons, NULL);
    vPolyRectBeginId.assign(num_polygons, std::numeric_limits<uint32_t>::max());

    for (uint32_t i = 0, ie = vPolyRectPattern.size(); i != ie; ++i)
    {
        const rectangle_pointer_type pPolyRectPattern = vPolyRectPattern[i];
        uint32_t parentPolygonId = vParentPolygonId[i];
        rectangle_pointer_type& pPattern = vPatternBbox[parentPolygonId];
        if (!pPattern)
        {
            pPattern = new rectangle_type (*pPolyRectPattern);
            pPattern->pattern_id(parentPolygonId); // vPatternBbox.pattern_id() is different from vPolyRectPattern.pattern_id()
            vPolyRectBeginId[parentPolygonId] = i;
        }
        else 
        {
            mplAssert(pPattern->color() == pPolyRectPattern->color());
            gtl::encompass(*pPattern, *pPolyRectPattern);
        }
        mplAssert(vPolyRectBeginId[parentPolygonId] <= i);
    }
}
LayoutDBPolygon::coordinate_difference LayoutDBPolygon::euclidean_distance(rectangle_type const& r1, rectangle_type const& r2) const 
{
    uint32_t polygon_id1 = r1.pattern_id();
    uint32_t polygon_id2 = r2.pattern_id();
    uint32_t num_polyrects = vPolyRectPattern.size();

    coordinate_difference distance = std::numeric_limits<LayoutDBPolygon::coordinate_difference>::max();
    uint32_t polyRectId1e = (polygon_id1+1 == vPolyRectBeginId.size())? num_polyrects : vPolyRectBeginId[polygon_id1+1];
    uint32_t polyRectId2e = (polygon_id2+1 == vPolyRectBeginId.size())? num_polyrects : vPolyRectBeginId[polygon_id2+1];
    for (uint32_t polyRectId1 = vPolyRectBeginId[polygon_id1]; polyRectId1 != polyRectId1e; ++polyRectId1)
        for (uint32_t polyRectId2 = vPolyRectBeginId[polygon_id2]; polyRectId2 != polyRectId2e; ++polyRectId2)
            distance = std::min(distance, (coordinate_difference)gtl::euclidean_distance(*vPolyRectPattern[polyRectId1], *vPolyRectPattern[polyRectId2]));

    return distance;
}
void LayoutDBPolygon::report_data() const 
{
    mplPrint(kINFO, "Input data for polygon based layout...\n");
    this->report_data_kernel();
    this->base_type::report_data_kernel();
}
void LayoutDBPolygon::report_data_kernel() const 
{
    mplPrint(kINFO, "# polygon rectangles = %lu\n", vPolyRectPattern.size());
}

// more techniques may be needed to improve the storage performance and reduce the peak memory.
void LayoutDBPolygon::refresh(std::vector<rectangle_pointer_type>& new_rect_vec, std::vector<uint32_t>& rect_to_parent)
{
#ifdef QDEBUG
	std::cout << "===== REFRESH =====" << std::endl;
	std::cout << "rect_to_parent : " << std::endl;
	std::cout << "size : " << rect_to_parent.size() << std::endl;
	//for(uint32_t i = 0; i < rect_to_parent.size(); i++)
	//	std::cout << rect_to_parent[i] << std::endl;
	std::cout << "new_rect_vec : " << std::endl;
	std::cout << "size : " << new_rect_vec.size() << std::endl;
	//for(uint32_t j = 0; j < new_rect_vec.size(); j++)
	//	std::cout << new_rect_vec[j]->pattern_id() << std::endl;
	std::cout << "===================" <<  std::endl;
#endif
	vParentPolygonId.clear();
	vParentPolygonId.assign(rect_to_parent.begin(), rect_to_parent.end());

	vPolyRectPattern.clear();
	vPolyRectPattern.assign(new_rect_vec.begin(), new_rect_vec.end());

	// update vPolyRectBeginId
	std::vector<uint32_t>().swap(vPolyRectBeginId);
	uint32_t base = vParentPolygonId[0];
	vPolyRectBeginId.push_back(0);

	for (uint32_t i = 1, ie = vParentPolygonId.size(); i < ie; i++)
	{
		if (vParentPolygonId[i] != base)
		{
			base = vParentPolygonId[i];
			vPolyRectBeginId.push_back(i);
		}
	}
	
#ifdef QDEBUG
	std::cout << "vPolyRectBeginId : " << std::endl;
	for (uint32_t i = 1, ie = vPolyRectBeginId.size(); i < ie; i++)
		std::cout << i << " -- " << vPolyRectBeginId[i] << std::endl;
#endif
	uint32_t num_polygons = vPolyRectBeginId.size();
	std::vector<rectangle_pointer_type>().swap(vPatternBbox);

	// vPatternBbox is the bounding box of a polygon 
	vPatternBbox.assign(num_polygons, NULL);
	
	for (uint32_t i = 0, ie = vPolyRectPattern.size(); i != ie; ++i)
	{
		const rectangle_pointer_type &pPolyRectPattern = vPolyRectPattern[i];
		uint32_t parentPolygonId = vParentPolygonId[i];
		rectangle_pointer_type& pPattern = vPatternBbox[parentPolygonId];
		if (!pPattern)
		{
			pPattern = new rectangle_type(*pPolyRectPattern);
			pPattern->pattern_id(parentPolygonId); // vPatternBbox.pattern_id() is different from vPolyRectPattern.pattern_id()
			std::cout <<"=============================" << std::endl;
			std::cout << "pPattern : " << pPattern->pattern_id() << " - " << +unsigned(pPattern->color()) << "\n" << "pPolyRectPattern : " << pPolyRectPattern->pattern_id() << " - " << +unsigned(pPolyRectPattern->color()) << std::endl << std::endl; 
		}
		else
		{
			std::cout <<"=============================" << std::endl;
			std::cout << "pPattern : " << pPattern->pattern_id() << " - " << +unsigned(pPattern->color()) << "\n" << "pPolyRectPattern : " << pPolyRectPattern->pattern_id() << " - " << +unsigned(pPolyRectPattern->color()) << std::endl << std::endl;  
			mplAssert(pPattern->color() == pPolyRectPattern->color());
			gtl::encompass(*pPattern, *pPolyRectPattern);
		}
		mplAssert(vPolyRectBeginId[parentPolygonId] <= i);
	}

	// no need to store tPatternBbox again, since we don't need to compute the distance between two Polyons again,
	// Also, due to the stitch insertion, the distanaces in rtree are invalid. I'm not sure whether to delete this.
	tPatternBbox.clear(); // tPatternBbox is used to store vPatternBbox construction with packing algorithm 
	rtree_type tTmp(vPatternBbox.begin(), vPatternBbox.end());
	tPatternBbox.swap(tTmp);
}

SIMPLEMPL_END_NAMESPACE
