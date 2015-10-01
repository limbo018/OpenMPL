/*************************************************************************
    > File Name: LayoutDBPolygon.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:49:21 AM CDT
 ************************************************************************/

#ifndef SIMPLEMPL_LAYOUTDBPOLYGON_H
#define SIMPLEMPL_LAYOUTDBPOLYGON_H

#include "LayoutDB.h"

SIMPLEMPL_BEGIN_NAMESPACE

/// polygon based layout 
struct LayoutDBPolygon : public LayoutDB
{
    typedef LayoutDB base_type;
	typedef base_type::coordinate_type coordinate_type;

    /// layout information
    /// base_type::vPatternBbox is used to store bounding boxes of polygons here 
    /// base_type::tPatternBbox is used to store poly rects during initialization, but for bounding boxes in vPatternBbox after that 
    std::vector<uint32_t> vParentPolygonId; ///< no need to actually store polygon patterns as we stored decomposed rectangles 
                                             ///< but we need to know which polygon a rectangle belongs to 
    std::vector<rectangle_pointer_type> vPolyRectPattern; ///< initial patterns decomposed from input polygons 
    std::vector<uint32_t> vPolyRectBeginId; ///< begin index in vPolyRectPattern when querying from parent polygon id 

	LayoutDBPolygon();
	LayoutDBPolygon(coordinate_type xl, coordinate_type yl, coordinate_type xh, coordinate_type yh);
	LayoutDBPolygon(LayoutDBPolygon const& rhs);
	virtual ~LayoutDBPolygon();
	LayoutDBPolygon& operator=(LayoutDBPolygon const& rhs);

	void copy(LayoutDBPolygon const& rhs);
	virtual void add_pattern(int32_t layer, std::vector<point_type> const& vPoint);
	/// call it to initialize rtree 
	/// it should be faster than gradually insertion 
	virtual void initialize_data();
    /// return poly rect patterns 
    virtual std::vector<rectangle_pointer_type> const& polyrect_patterns() const {return vPolyRectPattern;}
    /// set color for patterns 
    /// \param pattern_id is the index of vPatternBbox
    virtual void set_color(uint32_t pattern_id, int8_t color);
    /// mainly used for outputing edges 
    /// \return a point that is on the pattern and close to its center with given pattern id (polygon id for LayoutDBPolygon)
    /// default is to return the center of rectangle in vPatternBbox
    virtual point_type get_point_closest_to_center(uint32_t pattern_id) const; 
	virtual void report_data() const;
    void report_data_kernel() const;

    void compute_parent_polygons();
    void depth_first_search(uint32_t source, uint32_t polygon_id, uint32_t& order_id, std::vector<uint32_t>& vOrderId);
    void compute_parent_polygon_bboxes(uint32_t num_polygon);

    /// check vParentPolygonId for parent id 
    /// \return the euclidean distance of two patterns 
    virtual coordinate_difference euclidean_distance(rectangle_type const& r1, rectangle_type const& r2) const; 
};

SIMPLEMPL_END_NAMESPACE

#endif
