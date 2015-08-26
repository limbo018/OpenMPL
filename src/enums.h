/*************************************************************************
    > File Name: enums.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Mon Jun 15 21:43:53 2015
 ************************************************************************/

#ifndef SIMPLEMPL_ENUMS_H
#define SIMPLEMPL_ENUMS_H

#include <string>
#include <map>
#include <ostream>
#include "namespace.h"

SIMPLEMPL_BEGIN_NAMESPACE

/// base class for enumeration types 
template <typename EnumType>
class EnumExt
{
	public:
        typedef EnumType enum_type;
		EnumExt() {}
		EnumExt& operator=(EnumExt const& rhs)
		{
			if (this != &rhs)
				m_value = rhs.m_value;
			return *this;
		}
		EnumExt& operator=(enum_type const& rhs)
		{
			m_value = rhs;
			return *this;
		}
		EnumExt& operator=(std::string const& rhs)
        {
            m_value = str2Enum(rhs);
            return *this;
        }
		virtual operator std::string() const
        {
            return enum2Str(m_value);
        }
        EnumType const& get() const {return m_value;}

		bool operator==(EnumExt const& rhs) const {return m_value == rhs.m_value;}
		bool operator==(enum_type const& rhs) const {return m_value == rhs;}
		bool operator==(std::string const& rhs) const {return *this == EnumExt(rhs);}
		bool operator!=(EnumExt const& rhs) const {return m_value != rhs.m_value;}
		bool operator!=(enum_type const& rhs) const {return m_value != rhs;}
		bool operator!=(std::string const& rhs) const {return *this != EnumExt(rhs);}

		friend std::ostream& operator<<(std::ostream& os, const EnumExt& rhs)
		{
			rhs.print(os);
			return os;
		}
	protected:
		virtual void print(std::ostream& os) const {os << this->enum2Str(m_value);}

        virtual std::string enum2Str(enum_type const&) const = 0;
        virtual enum_type str2Enum(std::string const&) const = 0;

        EnumType m_value;
};

/// class AlgorithmType denotes type of algorithm for coloring  
struct AlgorithmTypeEnum
{
    enum EnumType {
		BACKTRACK = 0,    // no dependency 
		ILP_GURBOI = 1,   // only valid when gurobi is available
        ILP_CBC = 2,      // only valid when cbc is available
        LP_GUROBI = 3     // only valid when gurobi is available
    };
};
class AlgorithmType : public EnumExt<AlgorithmTypeEnum::EnumType>
{
	public:
        typedef AlgorithmTypeEnum enum_wrap_type;
        typedef enum_wrap_type::EnumType enum_type;
        typedef EnumExt<enum_type> base_type;

		AlgorithmType() : base_type() {m_value = enum_wrap_type::BACKTRACK;}
		AlgorithmType(AlgorithmType const& rhs) : base_type() {m_value = rhs.m_value;}
		AlgorithmType(enum_type const& rhs) : base_type() {m_value = rhs;}
		AlgorithmType(std::string const& rhs) : base_type() {m_value = str2Enum(rhs);}
		AlgorithmType& operator=(AlgorithmType const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}
		AlgorithmType& operator=(enum_type const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}
		AlgorithmType& operator=(std::string const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}

	protected:
        virtual std::string enum2Str(enum_type const& e) const;
        virtual enum_type str2Enum(std::string const& s) const;
};

/// Shape mode for database 
struct ShapeModeEnum
{
    enum EnumType {
        RECTANGLE = 0, // only non-overlapping rectangles 
        POLYGON = 1 // contain polygons or overlapping rectangles 
    };
};
class ShapeMode : public EnumExt<ShapeModeEnum::EnumType>
{
	public:
        typedef ShapeModeEnum enum_wrap_type;
        typedef enum_wrap_type::EnumType enum_type;
        typedef EnumExt<enum_type> base_type;

		ShapeMode() : base_type() {m_value = enum_wrap_type::RECTANGLE;}
		ShapeMode(ShapeMode const& rhs) : base_type() {m_value = rhs.m_value;}
		ShapeMode(enum_type const& rhs) : base_type() {m_value = rhs;}
		ShapeMode(std::string const& rhs) : base_type() {m_value = str2Enum(rhs);}
		ShapeMode& operator=(ShapeMode const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}
		ShapeMode& operator=(enum_type const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}
		ShapeMode& operator=(std::string const& rhs)
		{
            this->base_type::operator=(rhs);
			return *this;
		}

	protected:
        virtual std::string enum2Str(enum_type const& e) const;
        virtual enum_type str2Enum(std::string const& s) const;
};

SIMPLEMPL_END_NAMESPACE

#endif
