/*************************************************************************
    > File Name: LayoutDB.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:44:46 AM CDT
 ************************************************************************/

#ifndef SIMPLEMPL_LAYOUTDB_H
#define SIMPLEMPL_LAYOUTDB_H

#include <vector>
#include <map>
#include "shapes.h"
#include "parms.h"

SIMPLEMPL_BEGIN_NAMESPACE

/// base layout database class 
struct LayoutDB : public rectangle_data<int32_t>
{
    typedef int32_t coordinate_type;
	typedef gtl::coordinate_traits<coordinate_type>::manhattan_area_type area_type;
	typedef gtl::coordinate_traits<coordinate_type>::coordinate_difference coordinate_difference;
	typedef rectangle_data<coordinate_type> base_type;
	typedef point_data<coordinate_type> point_type;
	typedef segment_data<coordinate_type> path_type;
	typedef Rectangle<coordinate_type> rectangle_type;
	typedef Polygon<coordinate_type> polygon_type;
	typedef polygon_type* polygon_pointer_type;
	typedef rectangle_type* rectangle_pointer_type;
    // TO DO: experiment which hint is faster for rtree 
	//typedef bgi::rtree<rectangle_pointer_type, bgi::linear<16, 4> > rtree_type;
	typedef bgi::rtree<rectangle_pointer_type, bgi::rstar<16> > rtree_type;
	typedef polygon_90_set_data<coordinate_type> polygon_set_type;
    /// necessary for gtl::rectangle_traits
    using base_type::interval_type; 

	std::map<int32_t, std::vector<path_type> > hPath;    ///< path 
	std::string strname;                            ///< TOPCELL name, useful for dump out gds files 
	double unit;                               ///< keep output gdsii file has the same unit as input gdsii file 
	coordinate_difference coloring_distance;   ///< minimum coloring distance, set from coloring_distance_nm and unit

	/// options 
    ControlParameter parms; ///< control parameters from command line

	/// layout information 
    /// for rectangle-only layout, they denote patterns; otherwise, they denote decomposed rectangles from polygon patterns 
	rtree_type tPatternBbox;                       ///< rtree for components that intersects the LayoutDB
	std::vector<rectangle_pointer_type> vPatternBbox;   ///< uncolored and precolored patterns 

	struct compare_rectangle_type 
	{
		// by x and then by y
		bool operator() (rectangle_type const& r1, rectangle_type const& r2) const 
		{
			return gtl::xl(r1) < gtl::xl(r2) || (gtl::xl(r1) == gtl::xl(r2) && gtl::yl(r1) < gtl::yl(r2));
		}
		bool operator() (rectangle_pointer_type const& r1, rectangle_pointer_type const& r2) const 
		{return (*this)(*r1, *r2);}
	};

    LayoutDB();
	LayoutDB(coordinate_type xl, coordinate_type yl, coordinate_type xh, coordinate_type yh);
	LayoutDB(LayoutDB const& rhs);
	virtual ~LayoutDB();
	LayoutDB& operator=(LayoutDB const& rhs);

	void initialize();
	void copy(LayoutDB const& rhs);
    /// lib begins 
    virtual void begin_lib() {}
    /// cell begins 
    virtual void begin_str() {}
    /// lib ends 
    virtual void end_lib() {}
    /// str ends 
    virtual void end_str() {}
	virtual void add(int32_t layer, std::vector<point_type> const& vPoint);
    /// runtime disable this function instead of pure virtual function because we need to create LayoutDB object during initialization of SimpleMPL
	virtual void add_pattern(int32_t layer, std::vector<point_type> const& vPoint) = 0;
    /// add paths that indicating conflict edges to the layout 
	virtual void add_path(int32_t layer, std::vector<point_type> const& vPoint);
	/// call it to initialize rtree 
	/// it should be faster than gradually insertion 
	virtual void initialize_data();
    /// \return poly rect patterns 
    virtual std::vector<rectangle_pointer_type> const& polyrect_patterns() const = 0;
    /// \return patterns 
    virtual std::vector<rectangle_pointer_type> const& pattern_bboxes() const {return vPatternBbox;}
    /// set color for patterns 
    /// \param pattern_id is the index of vPatternBbox
    virtual void set_color(uint32_t pattern_id, int8_t color) = 0;
    /// mainly used for outputing edges 
    /// \return a point that is on the pattern and close to its center with given pattern id (polygon id for LayoutDBPolygon)
    /// default is to return the center of rectangle in vPatternBbox
    virtual point_type get_point_closest_to_center(uint32_t pattern_id) const = 0; 
    virtual void report_data() const = 0;
	void report_data_kernel() const;

    /// \return the euclidean distance of two patterns 
    virtual coordinate_difference euclidean_distance(rectangle_type const& r1, rectangle_type const& r2) const = 0;

    /// helper functions 
    /// update bounding box of layout 
    void update_bbox(base_type const& bbox);
    void check_layer_and_color(int32_t layer, bool& pattern_layer_flag, int8_t& color) const;
	/// remove overlapping patterns from vTargetPattern
	/// based on scanline approach, horizontal sort and then vertical sort 
	/// O(nlogn)
	void remove_overlap(std::vector<rectangle_pointer_type>& vTargetPattern);

    /// accesser functions
    inline double coloring_distance_nm() const {return parms.coloring_distance_nm;}
    inline int32_t color_num() const {return parms.color_num;}
    inline int32_t simplify_level() const {return parms.simplify_level;}
    inline int32_t thread_num() const {return parms.thread_num;}
    inline bool verbose() const {return parms.verbose;}
    inline uint32_t dbg_comp_id() const {return parms.dbg_comp_id;}
    inline AlgorithmType algo() const {return parms.algo;}
    inline ShapeMode shape_mode() const {return parms.shape_mode;}
    inline std::string const& input_gds() const {return parms.input_gds;}
    inline std::string const& output_gds() const {return parms.output_gds;}
};

SIMPLEMPL_END_NAMESPACE

#endif
