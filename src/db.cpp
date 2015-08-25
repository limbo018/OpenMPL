/*************************************************************************
    > File Name: db.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Mon 24 Aug 2015 03:43:08 PM CDT
 ************************************************************************/

#include "db.h"
#include <stack>

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
    for (std::vector<rectangle_pointer_type>::iterator it = vPattern.begin(), ite = vPattern.end(); it != ite; ++it)
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
    coloring_distance_nm     = 0;
    color_num                = 3;
    simplify_level           = 2;
    thread_num               = 1;
    verbose                  = false;
    input_gds                = "";
    output_gds               = "out.gds";
    algo                     = AlgorithmTypeEnum::BACKTRACK;
}
void LayoutDB::copy(LayoutDB const& rhs)
{
    hPath    = rhs.hPath;
    strname  = rhs.strname;
    unit     = rhs.unit;
    // options
    sUncolorLayer     = rhs.sUncolorLayer;
    sPrecolorLayer    = rhs.sPrecolorLayer;
    sPathLayer        = rhs.sPathLayer;
    coloring_distance = rhs.coloring_distance;
    coloring_distance_nm = rhs.coloring_distance_nm;
    color_num         = rhs.color_num;
    simplify_level    = rhs.simplify_level;
    thread_num        = rhs.thread_num;
    verbose           = rhs.verbose;
    input_gds         = rhs.input_gds;
    output_gds        = rhs.output_gds;
    algo              = rhs.algo;
    // layout information
    tPattern = rhs.tPattern;
    vPattern = rhs.vPattern;
}
void LayoutDB::add(int32_t layer, std::vector<point_type> const& vPoint)
{
    // classify features 
    // sometimes path may be defined as boundary 
    if (sPathLayer.count(layer))
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
    coloring_distance = (coordinate_difference)round(coloring_distance_nm/(unit*1e+9));
}
void LayoutDB::report_data() const
{
    mplPrint(kINFO, "Input data...\n");
    report_data_kernel();
}
void LayoutDB::report_data_kernel() const
{
    mplPrint(kINFO, "Total patterns # = %lu\n", vPattern.size());
    mplPrint(kINFO, "Coloring distance = %lld db ( %g nm )\n", coloring_distance, coloring_distance_nm);
    mplPrint(kINFO, "Color num = %d\n", color_num);
    mplPrint(kINFO, "Simplification level = %u\n", simplify_level);
    mplPrint(kINFO, "Thread num = %u\n", thread_num);
    mplPrint(kINFO, "Uncolored layer # = %lu", sUncolorLayer.size());
    if (!sUncolorLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = sUncolorLayer.begin(); it != sUncolorLayer.end(); ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Precolored layer # = %lu", sPrecolorLayer.size());
    if (!sPrecolorLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = sPrecolorLayer.begin(); it != sPrecolorLayer.end(); ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Path layer # = %lu", sPathLayer.size());
    if (!sPathLayer.empty())
    {
        mplPrint(kNONE, " ( ");
        for (std::set<int32_t>::const_iterator it = sPathLayer.begin(); it != sPathLayer.end(); ++it)
            mplPrint(kNONE, "%d ", *it);
        mplPrint(kNONE, ")");
    }
    mplPrint(kNONE, "\n");
    mplPrint(kINFO, "Algorithm = %s\n", std::string(algo).c_str());
}
void LayoutDB::update_bbox(base_type const& bbox)
{
    if (vPattern.empty())
        gtl::assign(*this, bbox);
    else 
        gtl::encompass(*this, bbox);
}
void LayoutDB::check_layer_and_color(int32_t layer, bool& pattern_layer_flag, int8_t& color) const
{
    pattern_layer_flag = true;
    color = -1;
    if (sPrecolorLayer.count(layer)) 
    {
        // for precolored patterns 
        // set colors 
        color = layer-*sPrecolorLayer.begin();
    }
    else if (sUncolorLayer.count(layer))
    {
        // for uncolored patterns
    }
    else pattern_layer_flag = false;
}
/// remove overlapping patterns 
/// based on scanline approach, horizontal sort and then vertical sort 
/// O(nlogn)
void LayoutDB::remove_overlap()
{
    std::sort(vPattern.begin(), vPattern.end(), compare_rectangle_type());
    // only duplicate is removed so far
    // TO DO: remove overlapping patterns 
    boost::dynamic_bitset<uint32_t, std::allocator<uint32_t> > vValid (vPattern.size()); // use a bit set to record validity 
    vValid.set(); // set all to 1
#ifdef DEBUG
    mplAssert(vValid[0] && vValid[vPattern.size()-1]);
#endif
    uint32_t duplicate_cnt = 0;
    for (std::vector<rectangle_pointer_type>::iterator it1 = vPattern.begin(), it2 = vPattern.begin();
            it2 != vPattern.end(); ++it2)
    {
        if (it2 != vPattern.begin()) 
        {
            rectangle_pointer_type& pPattern1 = *it1;
            rectangle_pointer_type& pPattern2 = *it2;

            if (gtl::equivalence(*pPattern1, *pPattern2)) 
            {
#ifdef DEBUG
                std::cout << "(W) " << *pPattern1 << " duplicates with " << *pPattern2 << " " << "ignored " << duplicate_cnt << std::endl;
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
    for (uint32_t i = 0; i < vPattern.size(); )
    {
        if (!vValid[vPattern[i]->pattern_id()])
        {
            std::swap(vPattern[i], vPattern.back());
            delete vPattern.back(); // recycle 
            vPattern.pop_back();
            invalid_cnt += 1;
        }
        else ++i;
    }
    mplAssert(duplicate_cnt == invalid_cnt);
#ifdef DEBUG
    for (uint32_t i = 0; i != vPattern.size(); ++i)
        mplAssert(vValid[vPattern[i]->pattern_id()]);
#endif
    // update pattern_id 
    for (uint32_t i = 0; i != vPattern.size(); ++i)
        vPattern[i]->pattern_id(i);
}


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
void LayoutDBRect::copy(LayoutDBRect const& rhs)
{
}
void LayoutDBRect::add_pattern(int32_t layer, std::vector<point_type> const& vPoint)
{
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
    pPattern->pattern_id(vPattern.size());

    // update layout boundary 
    update_bbox(*pPattern);

    // collect patterns 
    bool pattern_layer_flag = true;
    int8_t color = -1;
    check_layer_and_color(layer, pattern_layer_flag, color);

    if (pattern_layer_flag)
    {
        pPattern->color(color);
        // collect pattern 
        // initialize rtree later will contribute to higher efficiency in runtime 
        vPattern.push_back(pPattern);
    }
    else // recycle if it is not shared_ptr
        delete pPattern;
}
/// call it to initialize rtree 
/// it should be faster than gradually insertion 
void LayoutDBRect::initialize_data()
{
    this->base_type::initialize_data();
    // I assume there may be duplicate in the input gds, but no overlapping patterns 
    // duplicates are removed with following function
    remove_overlap();
    // construction with packing algorithm 
    rtree_type tTmp (vPattern.begin(), vPattern.end());
    tPattern.swap(tTmp);
}
void LayoutDBRect::report_data() const 
{
    mplPrint(kINFO, "Input data for rectangle based layout...\n");
    this->base_type::report_data_kernel();
}

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
}
void LayoutDBPolygon::add_pattern(int32_t layer, std::vector<point_type> const& vPoint)
{
    mplAssert(vPoint.size() >= 4);

    // collect patterns 
    bool pattern_layer_flag = true;
    int8_t color = -1;
    check_layer_and_color(layer, pattern_layer_flag, color);

    // skip other layers 
    if (!pattern_layer_flag)
        return;

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
        pPattern->pattern_id(vPattern.size());

        // collect pattern 
        // initialize rtree later will contribute to higher efficiency in runtime 
        vPattern.push_back(pPattern);
    }
}
/// call it to initialize rtree 
/// it should be faster than gradually insertion 
void LayoutDBPolygon::initialize_data()
{
    this->base_type::initialize_data();
    // I assume there may be duplicate in the input gds, but no overlapping patterns 
    // duplicates are removed with following function
    remove_overlap();
    // construction with packing algorithm 
    rtree_type tTmp (vPattern.begin(), vPattern.end());
    tPattern.swap(tTmp);
    // compute parent polygons for each rectangle by connected component algorithm 
    compute_parent_polygons();
}
void LayoutDBPolygon::compute_parent_polygons()
{
    vParentPolygonId.assign(vPattern.size(), std::numeric_limits<uint32_t>::max());
    uint32_t polygon_id = 0;
    for (uint32_t v = 0, ve = vPattern.size(); v != ve; ++v)
    {
        if (vParentPolygonId[v] == std::numeric_limits<uint32_t>::max()) // not visited 
            depth_first_search(v, polygon_id++);
    }
    num_polygons = polygon_id;
}
void LayoutDBPolygon::depth_first_search(uint32_t source, uint32_t polygon_id)
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
            const rectangle_pointer_type pPattern = vPattern[current];

			// find patterns connected with pPattern 
			// query tPattern 
            std::vector<rectangle_pointer_type> vAdjPattern;
			rectangle_type rect (*pPattern);
			// slightly bloat pPattern 
			gtl::bloat(rect, gtl::HORIZONTAL, 1);
			gtl::bloat(rect, gtl::VERTICAL, 1);
			for (rtree_type::const_query_iterator itq = tPattern.qbegin(bgi::intersects(rect));
					itq != tPattern.qend(); ++itq)
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
void LayoutDBPolygon::report_data() const 
{
    mplPrint(kINFO, "Input data for polygon based layout...\n");
    this->report_data_kernel();
    this->base_type::report_data_kernel();
}
void LayoutDBPolygon::report_data_kernel() const 
{
    mplPrint(kINFO, "# polygons = %u\n", num_polygons);
}

SIMPLEMPL_END_NAMESPACE

