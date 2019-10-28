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

/// ================================================================
/// LayoutDBPolygon implements the support to rectilinear polygon
/// based layout. With every input polygon, I first decompose it to
/// rectangles and store them. After reading completed, DFS is used 
/// to find connected component and re-union rectangles into polygons 
/// by identifying the parent polygon id for each rectangle. 
/// So we actually only stores decomposed rectangles instead of polygons. 
/// For memory efficiency, the rtree stores poly rects during reading, 
/// but switches to bounding box of each polygon when reading is finished
/// (literally after LayoutDBPolygon::initialize_data()).
/// ================================================================

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
                                                        ///< during reading process, it has the same order as that in the input file 
                                                        ///< but reordered so that rectangles belong to the same polygon are abutting
    std::vector<uint32_t> vPolyRectBeginId; ///< begin index in vPolyRectPattern when querying from parent polygon id 

    /// default constructor 
	LayoutDBPolygon();
    /// constructor
	LayoutDBPolygon(coordinate_type xl, coordinate_type yl, coordinate_type xh, coordinate_type yh);
    /// copy constructor 
	LayoutDBPolygon(LayoutDBPolygon const& rhs);
    /// destructor 
	virtual ~LayoutDBPolygon();
    /// assignment 
	LayoutDBPolygon& operator=(LayoutDBPolygon const& rhs);

	void copy(LayoutDBPolygon const& rhs);
	virtual void add_pattern(int32_t layer, std::vector<point_type> const& vPoint);
	/// call it to initialize rtree 
	/// it should be faster than gradually insertion 
	virtual void initialize_data();
    // used to generate vPolyRectPattern from PolygonSet
    void generate_rect_vectors();
    /// return poly rect patterns 
    virtual std::vector<rectangle_pointer_type> const& polyrect_patterns() const {return vPolyRectPattern;}
    ///\return a vector indicating the parent poly id of each rect 
    virtual std::vector<uint32_t> const& ParentPolygonId() const {return vParentPolygonId; }
	virtual std::vector<uint32_t> const& PolyRectBgnLoc() const { return vPolyRectBeginId; }
    /// set color for patterns 
    /// \param pattern_id is the index of vPatternBbox
    virtual void set_color(uint32_t pattern_id, int8_t color);
    /// mainly used for outputing edges 
    /// \return a point that is on the pattern and close to its center with given pattern id (polygon id for LayoutDBPolygon)
    /// default is to return the center of rectangle in vPatternBbox
    virtual point_type get_point_closest_to_center(uint32_t pattern_id) const; 
	virtual void report_data() const;
    void report_data_kernel() const;

    /// compute parent polygons from vPolyRectPattern
    /// then reorder vPolyRectPattern so that connected polygons are abutting in the array  
    void compute_parent_polygons();
    /// DFS to find connected component 
    void depth_first_search(uint32_t source, uint32_t polygon_id, uint32_t& order_id, std::vector<uint32_t>& vOrderId);
    /// from vPolyRectPattern compute parent polygons and their bounding boxes 
    void compute_parent_polygon_bboxes(uint32_t num_polygon);

    /// check vParentPolygonId for parent id 
    /// \return the euclidean distance of two patterns 
    virtual coordinate_difference euclidean_distance(rectangle_type const& r1, rectangle_type const& r2) const; 

	virtual void refresh(std::vector<rectangle_pointer_type>& new_rect_vec, std::vector<uint32_t>& rect_to_parent);
};

SIMPLEMPL_END_NAMESPACE

#endif
