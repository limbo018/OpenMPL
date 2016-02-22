/*************************************************************************
    > File Name: parms.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Sat 29 Aug 2015 11:12:45 AM CDT
 ************************************************************************/

#ifndef SIMPLEMPL_PARMS_H
#define SIMPLEMPL_PARMS_H

#include <set>
#include <boost/cstdint.hpp>
#include <limbo/programoptions/ProgramOptions.h>
#include "Msg.h"
#include "Enums.h"

SIMPLEMPL_BEGIN_NAMESPACE

using boost::uint32_t;

struct ControlParameter
{
	/// options 
	std::set<int32_t> sUncolorLayer;           ///< layers that represent uncolored patterns 
	std::set<int32_t> sPrecolorLayer;          ///< layers that represent precolored features, they should have the same number of colors 
	std::set<int32_t> sPathLayer;              ///< path layers that represent conflict edges 
	double coloring_distance_nm;               ///< minimum coloring distance in nanometer, set from command line 
	int32_t color_num;                         ///< number of colors available, only support 3 or 4
	int32_t simplify_level;                    ///< simplification level 0|1|2|3, default is 3
	int32_t thread_num;                        ///< number of maximum threads for parallel computation 
	bool verbose;                              ///< control screen message 
    uint32_t dbg_comp_id;                      ///< component id for debug, if matched, graphs will be dumped before and after coloring  

	std::string   input_gds;                   ///< input gdsii filename 
	std::string   output_gds;                  ///< output gdsii filename 

	/// algorithm options 
	AlgorithmType algo;                        ///< control algorithms used to solve coloring problem 
    /// database options 
    ShapeMode shape_mode;                      ///< it determins the actual derived layout database type 

    /// default constructor
    inline ControlParameter();

    /// swap with another ControlParameter
    inline void swap(ControlParameter& rhs);
};

inline ControlParameter::ControlParameter()
{
    coloring_distance_nm     = 0;
    color_num                = 3;
    simplify_level           = 3;
    thread_num               = 1;
    verbose                  = false;
    dbg_comp_id              = std::numeric_limits<uint32_t>::max();
    input_gds                = "";
    output_gds               = "";
    algo                     = AlgorithmTypeEnum::BACKTRACK;
    shape_mode               = ShapeModeEnum::RECTANGLE;
}

inline void ControlParameter::swap(ControlParameter& rhs)
{
    std::swap(sUncolorLayer, rhs.sUncolorLayer);
    std::swap(sPrecolorLayer, rhs.sPrecolorLayer);
    std::swap(sPathLayer, rhs.sPathLayer);
    std::swap(coloring_distance_nm, rhs.coloring_distance_nm);
    std::swap(color_num, rhs.color_num);
    std::swap(simplify_level, rhs.simplify_level);
    std::swap(thread_num, rhs.thread_num);
    std::swap(verbose, rhs.verbose);
    std::swap(dbg_comp_id, rhs.dbg_comp_id);
    std::swap(input_gds, rhs.input_gds);
    std::swap(output_gds, rhs.output_gds);
    std::swap(algo, rhs.algo);
    std::swap(shape_mode, rhs.shape_mode);
}

/// parse command line arguments 
struct CmdParser
{
	ControlParameter& parms;

	CmdParser(ControlParameter& p) : parms(p) {}

	bool operator()(int argc, char** argv);
};


SIMPLEMPL_END_NAMESPACE

#endif
