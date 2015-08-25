/*************************************************************************
    > File Name: db.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 06 Nov 2014 09:04:57 AM CST
 ************************************************************************/

#ifndef SIMPLEMPL_DB_H
#define SIMPLEMPL_DB_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <boost/version.hpp>

#if (BOOST_VERSION/100)%1000 > 55
// this is to fix the problem in boost 1.57.0 (1.55.0 works fine)
// it reports problem to find abs 
namespace boost { namespace polygon {
	using std::abs;
}} // namespace boost // namespace polygon
#endif

#include <boost/cstdint.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>
//#include <boost/polygon/polygon.hpp>
#include <boost/geometry.hpp>
// use adapted boost.polygon in boost.geometry, which is compatible to rtree
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "GeometryApi.h"
#include "msg.h"
#include "enums.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;
namespace bg = boost::geometry;
namespace bgi = bg::index;
using boost::int32_t;
using boost::int64_t;
using gtl::rectangle_concept;
using gtl::polygon_90_concept;
using gtl::polygon_90_set_concept;
using gtl::point_data;
using gtl::segment_data;
using gtl::rectangle_data;
using gtl::polygon_90_data;
using gtl::polygon_90_set_data;

using namespace gtl::operators;

class Shape 
{
    public:
		/// default constructor 
		Shape() {this->initialize();}
		/// copy constructor
		Shape(Shape const& rhs) {this->copy(rhs);}
		/// assignment 
		Shape& operator=(Shape const& rhs)
		{
			if (this != &rhs)
				this->copy(rhs);
			return *this;
		}
		virtual ~Shape() {}

#ifdef DEBUG
		long internal_id() {return m_internal_id;}
#endif

		int8_t color() const {return m_color;}
		void color(int8_t c) {m_color = c;}

		int32_t layer() const {return m_layer;}
		void layer(int32_t l) {m_layer = l;}

//		uint32_t comp_id() const {return m_comp_id;}
//		void comp_id(uint32_t c) {m_comp_id = c;}

		uint32_t pattern_id() const {return m_pattern_id;}
		void pattern_id(uint32_t p) {m_pattern_id = p;}

    private:
		void initialize()
		{
			m_color = m_layer = -1;
			//m_comp_id = std::numeric_limits<uint32_t>::max();
			m_pattern_id = std::numeric_limits<uint32_t>::max();
#ifdef DEBUG
#ifdef _OPENMP
#pragma omp critical 
#endif
			m_internal_id = generate_id();
#endif
		}
		void copy(Shape const& rhs)
		{
			this->m_color = rhs.m_color;
			this->m_layer = rhs.m_layer;
			//this->m_comp_id = rhs.m_comp_id;
			this->m_pattern_id = rhs.m_pattern_id;
#ifdef DEBUG
			this->m_internal_id = rhs.m_internal_id;
#endif
		}
		static long generate_id()
		{
			static long cnt = -1;
			mplAssert(cnt < std::numeric_limits<long>::max());
#ifdef _OPENMP
#pragma omp atomic 
#endif
			cnt += 1;
			return cnt;
		}

#ifdef DEBUG
		long m_internal_id; ///< internal id 
#endif

	protected:
		int32_t m_color : 4; ///< color, 4-bit is enough  
		int32_t m_layer : 28; ///< input layer, 28-bit is enough
		//uint32_t m_comp_id; ///< independent component id 
		uint32_t m_pattern_id; ///< index in the pattern array 
};

template <typename T>
class Rectangle : public rectangle_data<T>, public Shape
{
	public:
		typedef T coordinate_type;
		typedef rectangle_data<coordinate_type> base_type;
        typedef Shape shape_base_type;
		typedef rectangle_concept geometry_type; // important 
		typedef point_data<coordinate_type> point_type;
		using typename base_type::interval_type;

		/// default constructor 
		Rectangle() : base_type(), shape_base_type() {}
		Rectangle(interval_type const& hor, interval_type const& ver) : base_type(hor, ver), shape_base_type() {}
		Rectangle(coordinate_type xl, coordinate_type yl, coordinate_type xh, coordinate_type yh) : base_type(xl, yl, xh, yh), shape_base_type() {}
		/// copy constructor
		Rectangle(Rectangle const& rhs) : base_type(rhs), shape_base_type(rhs) {}
		/// assignment 
		Rectangle& operator=(Rectangle const& rhs)
		{
			if (this != &rhs)
			{
				this->base_type::operator=(rhs);
                this->shape_base_type::operator=(rhs);
			}
			return *this;
		}
		virtual ~Rectangle() {}

		friend std::ostream& operator<<(std::ostream& os, Rectangle const& rhs)
		{
			os << "(" << gtl::xl(rhs) << ", " << gtl::yl(rhs) << ", " << gtl::xh(rhs) << ", " << gtl::yh(rhs) << ")";
			return os;
		}
};

template <typename T>
class Polygon : public polygon_90_data<T>, public Shape
{
	public:
		typedef T coordinate_type;
		typedef polygon_90_data<coordinate_type> base_type;
        typedef Shape shape_base_type;
		typedef Rectangle<coordinate_type> rectangle_type;
		typedef point_data<coordinate_type> point_type;
		using typename base_type::geometry_type;
		using typename base_type::compact_iterator_type;
		using typename base_type::iterator_type;
		using typename base_type::area_type;

		/// default constructor 
		Polygon() : base_type(), shape_base_type() {}
		/// copy constructor
		Polygon(Polygon const& rhs) : base_type(rhs), shape_base_type(rhs) {}
		/// assignment 
		Polygon& operator=(Polygon const& rhs)
		{
			if (this != &rhs)
			{
				this->base_type::operator=(rhs);
                this->shape_base_type::operator=(rhs);
			}
			return *this;
		}
		virtual ~Polygon() {}

		friend std::ostream& operator<<(std::ostream& os, Polygon const& rhs)
		{
            os << "(";
            for (iterator_type it = rhs.begin(), ite = rhs.end(); it != ite; ++it)
                os << "(" << it->x() << ", " << it->y() << ")";
            os << ")";
			return os;
		}
};

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

	/// options 
	std::set<int32_t> sUncolorLayer;                ///< layers that represent uncolored patterns 
	std::set<int32_t> sPrecolorLayer;               ///< layers that represent precolored features, they should have the same number of colors 
	std::set<int32_t> sPathLayer;                   ///< path layers that represent conflict edges 
	coordinate_difference coloring_distance;   ///< minimum coloring distance, std::set from coloring_distance_nm and unit
	double coloring_distance_nm;               ///< minimum coloring distance in nanometer, std::set from command line 
	int32_t color_num;                         ///< number of colors available, only support 3 or 4
	int32_t simplify_level;                    ///< simplification level 0|1|2, default is 2
	int32_t thread_num;                        ///< number of maximum threads for parallel computation 
	bool verbose;                              ///< control screen message 

	std::string   input_gds;                        ///< input gdsii filename 
	std::string   output_gds;                       ///< output gdsii filename 

	/// algorithm options 
	AlgorithmType algo;                        ///< control algorithms used to solve coloring problem 

	/// layout information 
    /// for rectangle-only layout, they denote patterns; otherwise, they denote decomposed rectangles from polygon patterns 
	rtree_type tPattern;                       ///< rtree for components that intersects the LayoutDB
	std::vector<rectangle_pointer_type> vPattern;   ///< uncolored and precolored patterns 

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
	virtual void add(int32_t layer, std::vector<point_type> const& vPoint);
	virtual void add_pattern(int32_t layer, std::vector<point_type> const& vPoint) = 0;
	virtual void add_path(int32_t layer, std::vector<point_type> const& vPoint);
	virtual void initialize_data();
    virtual void report_data() const;
	void report_data_kernel() const;

    /// \return true if two patterns belong to the same parent polygon
    virtual bool same_parent(uint32_t id1, uint32_t id2) const = 0;

    /// helper functions 
    /// update bounding box of layout 
    void update_bbox(base_type const& bbox);
    void check_layer_and_color(int32_t layer, bool& pattern_layer_flag, int8_t& color) const;
	/// remove overlapping patterns 
	/// based on scanline approach, horizontal sort and then vertical sort 
	/// O(nlogn)
	void remove_overlap();
};

/// current implementation assume all the input patterns are rectangles 
struct LayoutDBRect : public LayoutDB
{
    typedef LayoutDB base_type;
	typedef base_type::coordinate_type coordinate_type;

	LayoutDBRect();
	LayoutDBRect(coordinate_type xl, coordinate_type yl, coordinate_type xh, coordinate_type yh);
	LayoutDBRect(LayoutDBRect const& rhs);
	virtual ~LayoutDBRect();
	LayoutDBRect& operator=(LayoutDBRect const& rhs);

	void copy(LayoutDBRect const& rhs);
	virtual void add_pattern(int32_t layer, std::vector<point_type> const& vPoint);
	/// call it to initialize rtree 
	/// it should be faster than gradually insertion 
	virtual void initialize_data();
	virtual void report_data() const;

    /// always return false as each rectangle is its own parent 
    virtual bool same_parent(uint32_t, uint32_t) const {return false;}
};

/// polygon based layout 
struct LayoutDBPolygon : public LayoutDB
{
    typedef LayoutDB base_type;
	typedef base_type::coordinate_type coordinate_type;

    /// layout information
    std::vector<uint32_t> vParentPolygonId; ///< no need to actually store polygon patterns as we stored decomposed rectangles 
                                             ///< but we need to know which polygon a rectangle belongs to 
    uint32_t num_polygons; ///< number of polygons 

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
	virtual void report_data() const;
    void report_data_kernel() const;

    void compute_parent_polygons();
    void depth_first_search(uint32_t source, uint32_t polygon_id);

    /// check vParentPolygonId for parent id 
    virtual bool same_parent(uint32_t id1, uint32_t id2) const {return vParentPolygonId[id1] == vParentPolygonId[id2];}
};

SIMPLEMPL_END_NAMESPACE

#endif 
