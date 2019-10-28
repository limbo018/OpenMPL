/*************************************************************************
    > File Name: Params.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed 07 Oct 2015 06:58:01 PM CDT
 ************************************************************************/

#include "Params.h"

SIMPLEMPL_BEGIN_NAMESPACE

bool CmdParser::operator()(int argc, char** argv)
{
    ControlParameter defaultParms; // get default value from default constructor 
    bool help = false;
    std::string algo_str;
    std::string shape_str;
    // append options here 
    typedef limbo::programoptions::ProgramOptions po_type;
    using limbo::programoptions::Value;
    po_type desc (std::string("SimpleMPL 1.1 Usage"));
    desc.add_option(Value<bool>("-help", &help, "toggle printing help message").toggle(true).default_value(false).toggle_value(true).help(true))
        .add_option(Value<std::string>("-in", &parms.input_gds, "input gds file name").required(true))
        .add_option(Value<std::string>("-out", &parms.output_gds, "output gds file name").default_value(defaultParms.output_gds))
        .add_option(Value<double>("-coloring_distance", &parms.coloring_distance_nm, "a floating point number indicating number of coloring distance in nanometer").default_value(defaultParms.coloring_distance_nm))
        .add_option(Value<int32_t>("-color_num", &parms.color_num, "an integer indicating number of masks (colors) < 3|4 >").required(true))
        .add_option(Value<int32_t>("-simplify_level", &parms.simplify_level, "an integer indicating graph simplification level < 0|1|2|3 >").default_value(defaultParms.simplify_level))
        .add_option(Value<int32_t>("-thread_num", &parms.thread_num, "an integer indicating maximum thread number").default_value(defaultParms.thread_num))
        .add_option(Value<std::set<int32_t> >("-path_layer", &parms.sPathLayer, "an integer indicating layer for conflict edges"))
        .add_option(Value<std::set<int32_t> >("-precolor_layer", &parms.sPrecolorLayer, "an integer indicating layer for pre-colored patterns"))
        .add_option(Value<std::set<int32_t> >("-uncolor_layer", &parms.sUncolorLayer, "an integer indicating layer for coloring"))
        .add_option(Value<std::string>("-algo", &algo_str, "algorithm type < ILP|BACKTRACK|LP|SDP >").default_value(std::string(defaultParms.algo)))
        .add_option(Value<std::string>("-shape", &shape_str, "shape mode < RECTANGLE|POLYGON >").default_value(std::string(defaultParms.shape_mode)))
        .add_option(Value<bool>("-verbose", &parms.verbose, "toggle controling screen messages").toggle(true).default_value(defaultParms.verbose).toggle_value(true))
        .add_option(Value<bool>("-use_stitch", &parms.use_stitch, "whether use stitch projection").toggle(true).default_value(defaultParms.use_stitch).toggle_value(true))
		.add_option(Value<bool>("-gen_stitch", &parms.gen_stitch, "control whether only generate and output stitches").toggle(true).default_value(defaultParms.gen_stitch).toggle(true))
		.add_option(Value<uint32_t>("-dbg_comp_id", &parms.dbg_comp_id, "debug component id").default_value(defaultParms.dbg_comp_id))
        .add_option(Value<int32_t>("-record", &parms.record, "record level, which controls the degree of recording information").default_value(defaultParms.record))
        .add_option(Value<double>("-weighted_stitch", &parms.weight_stitch, "a floating point number indicating the weight of stitch").default_value(defaultParms.weight_stitch))
        ;
    try
    {
        desc.parse(argc, argv);

        // print help message 
        if (help)
        {
            std::cout << desc << "\n";
            exit(1);
        }

        // post processing algo_str 
        if (limbo::iequals(algo_str, "ILP")) 
        {
#if GUROBI == 1
            parms.algo = AlgorithmTypeEnum::ILP_GURBOI;
#elif LEMONCBC == 1
            parms.algo = AlgorithmTypeEnum::ILP_CBC;
#else 
            mplPrint(kWARN, "ILP is not available without GUROBI or CBC, set to default algorithm\n");
#endif
        }
        if (limbo::iequals(algo_str, "ILP_UPDATED")) 
        {
#if GUROBI == 1
            parms.algo = AlgorithmTypeEnum::ILP_UPDATED_GURBOI;
#else 
            mplPrint(kWARN, "ILP updated is not available without GUROBI, set to default algorithm\n");
#endif
        }
        else if (limbo::iequals(algo_str, "LP"))
        {
#if GUROBI == 1
            parms.algo = AlgorithmTypeEnum::LP_GUROBI;
#else 
            mplPrint(kWARN, "LP is not available without GUROBI, set to default algorithm\n");
#endif
        }
        else if (limbo::iequals(algo_str, "SDP"))
        {
#if CSDP == 1
            parms.algo = AlgorithmTypeEnum::SDP_CSDP;
#else 
            mplPrint(kWARN, "SDP is not available without CSDP, set to default algorithm\n");
#endif
        }
        else if (limbo::iequals(algo_str, "MIS"))
        {
#if GUROBI == 1
            parms.algo = AlgorithmTypeEnum::MIS_GUROBI; 
#else 
            mplPrint(kWARN, "MIS is not available without GUROBI, set to default algorithm\n");
#endif
        }
        else if (limbo::iequals(algo_str, "BACKTRACK"))
            parms.algo = AlgorithmTypeEnum::BACKTRACK;
        else if (limbo::iequals(algo_str, "DL"))
            parms.algo = AlgorithmTypeEnum::DANCING_LINK;
        else mplPrint(kWARN, "Unknown algorithm type %s, set to default algorithm\n", algo_str.c_str());

        // post processing shape_str
        if (limbo::iequals(shape_str, "RECTANGLE"))
            parms.shape_mode = ShapeModeEnum::RECTANGLE;
        else if (limbo::iequals(shape_str, "POLYGON"))
            parms.shape_mode = ShapeModeEnum::POLYGON;
        else mplPrint(kWARN, "Unknown shape mode %s, set to default\n", shape_str.c_str());

        // check condition 
        mplAssertMsg(parms.coloring_distance_nm > 0 || !parms.sPathLayer.empty(), "should set positive coloring_distance_nm or specify path_layer for conflict edges");
    }
    catch (std::exception& e)
    {
        // print help message and error message 
        std::cout << desc << "\n";
        mplPrint(kERROR, "%s\n", e.what());
        return false;
    }
    return true;
}


SIMPLEMPL_END_NAMESPACE
