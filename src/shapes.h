/*************************************************************************
    > File Name: shapes.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 06 Nov 2014 09:04:57 AM CST
 ************************************************************************/

#ifndef SIMPLEMPL_SHAPES_H
#define SIMPLEMPL_SHAPES_H

#include <iostream>
#include <string>
#include <cmath> // std::abs
#include <boost/version.hpp>

#if (BOOST_VERSION/100)%1000 > 55
// this is to fix the problem in boost 1.57.0 (1.55.0 works fine)
// it reports problem to find abs 
namespace boost { namespace polygon {
	using std::abs;
}} // namespace boost // namespace polygon
#endif

#include <boost/cstdint.hpp>
//#include <boost/polygon/polygon.hpp>
#include <boost/geometry.hpp>
// use adapted boost.polygon in boost.geometry, which is compatible to rtree
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "GeometryApi.h"
#include "msg.h"

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

		uint32_t pattern_id() const {return m_pattern_id;}
		void pattern_id(uint32_t p) {m_pattern_id = p;}

    private:
		void initialize()
		{
			m_color = m_layer = -1;
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
        /// convertion to std::string 
        operator std::string() const 
        {
            std::ostringstream oss;
            print(oss);
            return oss.str();
        }
		virtual ~Rectangle() {}

        void print(std::ostream& os) const 
        {
			os << "(" << gtl::xl(*this) << ", " << gtl::yl(*this) << ", " << gtl::xh(*this) << ", " << gtl::yh(*this) << ")";
        }
		friend std::ostream& operator<<(std::ostream& os, Rectangle const& rhs)
		{
            rhs.print(os);
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
        /// convertion to std::string 
        operator std::string() const 
        {
            std::ostringstream oss;
            print(oss);
            return oss.str();
        }
		virtual ~Polygon() {}

        void print(std::ostream& os) const
        {
            os << "(";
            for (iterator_type it = this->begin(), ite = this->end(); it != ite; ++it)
                os << "(" << it->x() << ", " << it->y() << ")";
            os << ")";
        }
		friend std::ostream& operator<<(std::ostream& os, Polygon const& rhs)
		{
            rhs.print(os);
			return os;
		}
};

SIMPLEMPL_END_NAMESPACE

#endif 
