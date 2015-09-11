/*************************************************************************
    > File Name: io.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed 26 Aug 2015 10:59:58 AM CDT
 ************************************************************************/

#include "io.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;

void GdsWriter::operator() (std::string const& filename, GdsWriter::layoutdb_type const& db, 
        std::vector<std::pair<uint32_t, uint32_t> > const& vConflict, 
        std::vector<std::vector<uint32_t> > const& mAdjVertex, 
        std::string const& strname, double unit) const 
{
    GdsParser::GdsWriter gw (filename.c_str());
    gw.gds_create_lib("POLYGONS", unit /* um per bit */ );
    gw.gds_write_bgnstr();
    gw.gds_write_strname(strname.c_str());

    // if there are precolored patterns, keep the same layer convention 
    int32_t layer_offset = (db.parms.sPrecolorLayer.empty())? 100 : *db.parms.sPrecolorLayer.begin();
    // basic operation
    // will add more 
    write_rectangles(gw, db.polyrect_patterns(), layer_offset);
    write_conflicts(gw, db, vConflict,  layer_offset+db.color_num());   // conflict layer 
    write_edges(gw, db, mAdjVertex, layer_offset+db.color_num()+1); // draw edges 
    //(*this)(gw, db.hPath); // draw edges if there exits 

    gw.gds_write_endstr();
    gw.gds_write_endlib(); 
}
void GdsWriter::write_rectangles(GdsParser::GdsWriter& gw, std::vector<GdsWriter::rectangle_pointer_type> const& vRect, const int32_t layer_offset) const 
{
    for (std::vector<rectangle_pointer_type>::const_iterator it = vRect.begin(); it != vRect.end(); ++it)
    {
        rectangle_type const& rect = **it;
        gw.write_box(layer_offset+rect.color(), 0, 
                gtl::xl(rect), gtl::yl(rect), 
                gtl::xh(rect), gtl::yh(rect));
#ifdef DEBUG
        mplAssert(rect.color() >= 0);
#endif
    }
}
void GdsWriter::write_conflicts(GdsParser::GdsWriter& gw, GdsWriter::layoutdb_type const& db, 
        std::vector<std::pair<uint32_t, uint32_t> > const& vConflict, const int32_t layer) const
{
    for (std::vector<std::pair<uint32_t, uint32_t> >::const_iterator it = vConflict.begin(); it != vConflict.end(); ++it)
    {
        // create a path
        gw.gds_write_path();
        gw.gds_write_layer(layer);
        gw.gds_write_datatype(0);
        gw.gds_write_pathtype(2); // extended square ends
        gw.gds_write_width(5); // 5 nm wide

        point_type vc = db.get_point_closest_to_center(it->first);
        point_type uc = db.get_point_closest_to_center(it->second);
        int32_t x[2] = {gtl::x(vc), gtl::x(uc)};
        int32_t y[2] = {gtl::y(vc), gtl::y(uc)};

        gw.gds_write_xy(x, y, 2);
        gw.gds_write_endel();
    }
}
void GdsWriter::write_paths(GdsParser::GdsWriter& gw, std::map<int32_t, std::vector<GdsWriter::path_type> > const& hPath) const 
{
    for (std::map<int32_t, std::vector<path_type> >::const_iterator it1 = hPath.begin(); it1 != hPath.end(); ++it1)
    {
        const int32_t layer = it1->first;
        std::vector<path_type> const& vPath = it1->second;
        for (std::vector<path_type>::const_iterator it2 = vPath.begin(); it2 != vPath.end(); ++it2)
        {
            path_type const& path = *it2;
            // create a path
            gw.gds_write_path();
            gw.gds_write_layer(layer);
            gw.gds_write_datatype(0);
            gw.gds_write_pathtype(2); // extended square ends
            gw.gds_write_width(5); // 5 nm wide

            int32_t x[2] = {gtl::x(gtl::low(path)), gtl::x(gtl::high(path))};
            int32_t y[2] = {gtl::y(gtl::low(path)), gtl::y(gtl::high(path))};

            gw.gds_write_xy(x, y, 2);
            gw.gds_write_endel();
        }
    }
}
void GdsWriter::write_edges(GdsParser::GdsWriter& gw, GdsWriter::layoutdb_type const& db, std::vector<std::vector<uint32_t> > const& mAdjVertex, const int32_t layer) const 
{
    for (uint32_t i = 0; i != mAdjVertex.size(); ++i)
    {
        for (uint32_t j = 0; j != mAdjVertex[i].size(); ++j)
        {
            uint32_t v = i;
            uint32_t u = mAdjVertex[i][j];
            // create a path from v to u 
            if (v < u) // avoid duplicate 
            {
                // create a path
                gw.gds_write_path();
                gw.gds_write_layer(layer);
                gw.gds_write_datatype(0);
                gw.gds_write_pathtype(2); // extended square ends
                gw.gds_write_width(5); // 5 nm wide

                point_type vc = db.get_point_closest_to_center(v);
                point_type uc = db.get_point_closest_to_center(u);
                int32_t x[2] = {gtl::x(vc), gtl::x(uc)};
                int32_t y[2] = {gtl::y(vc), gtl::y(uc)};

                gw.gds_write_xy(x, y, 2);
                gw.gds_write_endel();
            }
        }
    }
}

bool CmdParser::operator()(int argc, char** argv)
{
    ControlParameter defaultParms; // get default value from default constructor 
    bool help = false;
    std::string algo_str;
    std::string shape_str;
    // append options here 
    typedef limbo::programoptions::ProgramOptions po_type;
    using limbo::programoptions::Value;
    po_type desc (std::string("SimpleMPL 1.X Usage"));
    desc.add_option(Value<bool>("-help", &help, "print help message").toggle(true).default_value(false).toggle_value(true).help(true))
        .add_option(Value<std::string>("-in", &parms.input_gds, "input gds file name").required(true))
        .add_option(Value<std::string>("-out", &parms.output_gds, "output gds file name").default_value(defaultParms.output_gds))
        .add_option(Value<double>("-coloring_distance", &parms.coloring_distance_nm, "a floating point number indicating number of coloring distance in nanometer").default_value(defaultParms.coloring_distance_nm))
        .add_option(Value<int32_t>("-color_num", &parms.color_num, "an integer indicating number of masks (colors)").required(true))
        .add_option(Value<int32_t>("-simplify_level", &parms.simplify_level, "an integer indicating graph simplification level < 0|1|2 >").default_value(defaultParms.simplify_level))
        .add_option(Value<int32_t>("-thread_num", &parms.thread_num, "an integer indicating maximum thread number").default_value(defaultParms.thread_num))
        .add_option(Value<std::set<int32_t> >("-path_layer", &parms.sPathLayer, "an integer indicating layer for conflict edges"))
        .add_option(Value<std::set<int32_t> >("-precolor_layer", &parms.sPrecolorLayer, "an integer indicating layer for pre-colored patterns"))
        .add_option(Value<std::set<int32_t> >("-uncolor_layer", &parms.sUncolorLayer, "an integer indicating layer for coloring"))
        .add_option(Value<std::string>("-algo", &algo_str, "algorithm type < ILP|BACKTRACK >").default_value(std::string(defaultParms.algo)))
        .add_option(Value<std::string>("-shape", &shape_str, "shape mode < RECTANGLE|POLYGON >").default_value(std::string(defaultParms.shape_mode)))
        .add_option(Value<bool>("-verbose", &parms.verbose, "control screen messages").toggle(true).default_value(defaultParms.verbose).toggle_value(true))
        .add_option(Value<uint32_t>("-dbg_comp_id", &parms.dbg_comp_id, "debug component id").default_value(defaultParms.dbg_comp_id))
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
        else if (limbo::iequals(algo_str, "BACKTRACK"))
            parms.algo = AlgorithmTypeEnum::BACKTRACK;
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
