/*************************************************************************
    > File Name: enums.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Mon 03 Aug 2015 10:50:21 AM CDT
 ************************************************************************/

#include "enums.h"

SIMPLEMPL_BEGIN_NAMESPACE

#ifndef ENUM2STR
#define ENUM2STR(map, var) \
    map[enum_wrap_type::var] = #var
#endif

#ifndef STR2ENUM
#define STR2ENUM(map, var) \
    map[#var] = enum_wrap_type::var
#endif

std::string AlgorithmType::enum2Str(AlgorithmType::enum_type const& e) const
{
    static std::map<enum_type, std::string> mEnum2Str;
    static bool init = true;

    if (init)
    {
        ENUM2STR(mEnum2Str, BACKTRACK);
        ENUM2STR(mEnum2Str, ILP_GURBOI);
        ENUM2STR(mEnum2Str, ILP_CBC);
        ENUM2STR(mEnum2Str, LP_GUROBI);
        ENUM2STR(mEnum2Str, SDP_CSDP);
        init = false;
    }

    return mEnum2Str.at(e);
}

AlgorithmType::enum_type AlgorithmType::str2Enum(std::string const& s) const
{
    static std::map<std::string, enum_type> mStr2Enum;
    static bool init = true;

    if (init)
    {
        STR2ENUM(mStr2Enum, BACKTRACK);
        STR2ENUM(mStr2Enum, ILP_GURBOI);
        STR2ENUM(mStr2Enum, ILP_CBC);
        STR2ENUM(mStr2Enum, LP_GUROBI);
        STR2ENUM(mStr2Enum, SDP_CSDP);
        init = false;
    }

    return mStr2Enum.at(s);
}

std::string ShapeMode::enum2Str(ShapeMode::enum_type const& e) const
{
    static std::map<enum_type, std::string> mEnum2Str;
    static bool init = true;

    if (init)
    {
        ENUM2STR(mEnum2Str, RECTANGLE);
        ENUM2STR(mEnum2Str, POLYGON);
        init = false;
    }

    return mEnum2Str.at(e);
}

ShapeMode::enum_type ShapeMode::str2Enum(std::string const& s) const
{
    static std::map<std::string, enum_type> mStr2Enum;
    static bool init = true;

    if (init)
    {
        STR2ENUM(mStr2Enum, RECTANGLE);
        STR2ENUM(mStr2Enum, POLYGON);
        init = false;
    }

    return mStr2Enum.at(s);
}
SIMPLEMPL_END_NAMESPACE
