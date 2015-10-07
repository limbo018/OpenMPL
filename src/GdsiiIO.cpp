/*************************************************************************
    > File Name: GdsiiIO.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Wed 26 Aug 2015 10:59:58 AM CDT
 ************************************************************************/

#include "GdsiiIO.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;

bool GdsReader::operator() (std::string const& filename)  
{
    // calculate file size 
    std::ifstream in (filename.c_str());
    if (!in.good()) return false;
    std::streampos begin = in.tellg();
    in.seekg(0, std::ios::end);
    std::streampos end = in.tellg();
    file_size = (end-begin);
    in.close();
    // read gds 
    return GdsParser::read(*this, filename);
}

void GdsReader::bit_array_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vBitArray)
{
    this->integer_cbk(record_type, data_type, vBitArray);
}
void GdsReader::integer_2_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vInteger)
{
    this->integer_cbk(record_type, data_type, vInteger);
}
void GdsReader::integer_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vInteger)
{
    this->integer_cbk(record_type, data_type, vInteger);
}
void GdsReader::real_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vFloat) 
{
    this->float_cbk(record_type, data_type, vFloat);
}
void GdsReader::real_8_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vFloat) 
{
    this->float_cbk(record_type, data_type, vFloat);
}
void GdsReader::string_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::string const& str) 
{
    mplAssert(data_type == GdsParser::GdsData::STRING);
    switch (record_type)
    {
        case GdsParser::GdsRecords::STRNAME:
			db.strname.assign(str);
            break;
        case GdsParser::GdsRecords::LIBNAME:
        case GdsParser::GdsRecords::STRING:
        default: break;
    }
}
void GdsReader::begin_end_cbk(GdsParser::GdsRecords::EnumType record_type)
{
    switch (record_type)
    {
        case GdsParser::GdsRecords::BOX:
        case GdsParser::GdsRecords::BOUNDARY:
        case GdsParser::GdsRecords::PATH:
            vPoint.clear();
            layer = 0;
            status = record_type;
            break;
        case GdsParser::GdsRecords::ENDEL:
            {
                switch (status)
                {
                    case GdsParser::GdsRecords::BOX:
                    case GdsParser::GdsRecords::BOUNDARY:
                    case GdsParser::GdsRecords::PATH:
                        mplAssert(layer != -1);
                        db.add(layer, vPoint);
                        break;
                    default: break;
                }
                status = GdsParser::GdsRecords::UNKNOWN;
            }
            break;
        case GdsParser::GdsRecords::ENDLIB: // notify database on the end of lib 
            db.end_lib();
            break;
        case GdsParser::GdsRecords::ENDSTR: // currently not interested, add stuff here if needed 
            db.end_str();
            break;
        default: // be careful here, you may dump a lot of unnecessary error message for unknown record_type 
            mplPrint(kERROR, "%s() unsupported record_type = %s", __func__, GdsParser::gds_record_ascii(record_type));
            break;
    }
}
/// helper functions 
void GdsReader::integer_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType /*data_type*/, std::vector<int> const& vData)
{
    switch (record_type)
    {
        case GdsParser::GdsRecords::LAYER:
            layer = vData[0];
            break;
        case GdsParser::GdsRecords::XY:
            if (status == GdsParser::GdsRecords::BOX || status == GdsParser::GdsRecords::BOUNDARY || status == GdsParser::GdsRecords::PATH)
            {
				mplAssert((vData.size() & 1) == 0 && vData.size() >= 4);
				vPoint.clear();
				uint32_t end = vData.size();
				// skip last point for BOX and BOUNDARY
				if (status == GdsParser::GdsRecords::BOX || status == GdsParser::GdsRecords::BOUNDARY) end -= 2;
				for (uint32_t i = 0; i < end; i += 2)
					vPoint.push_back(gtl::construct<point_type>(vData[i], vData[i+1]));
            }
            break;
        case GdsParser::GdsRecords::BGNLIB: // notify database on the begin of lib 
            db.begin_lib(); 
            break;
        case GdsParser::GdsRecords::BGNSTR: // just date of creation, not interesting
            db.begin_str();
            break;
        default: // other not interested record_type
            //mplPrint(kERROR, "%s() invalid record_type = %s, data_type = %s", __func__, GdsParser::gds_record_ascii(record_type), GdsParser::gds_data_ascii(data_type));
            break;
    }
}
void GdsReader::float_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vData)
{
    switch (record_type)
    {
        case GdsParser::GdsRecords::UNITS:
			db.unit = vData[1]; 
            break;
        default:
            mplPrint(kERROR, "%s() invalid record_type = %s, data_type = %s", __func__, GdsParser::gds_record_ascii(record_type), GdsParser::gds_data_ascii(data_type));
            break;
    }
}

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
    desc.add_option(Value<bool>("-help", &help, "toggle printing help message").toggle(true).default_value(false).toggle_value(true).help(true))
        .add_option(Value<std::string>("-in", &parms.input_gds, "input gds file name").required(true))
        .add_option(Value<std::string>("-out", &parms.output_gds, "output gds file name").default_value(defaultParms.output_gds))
        .add_option(Value<double>("-coloring_distance", &parms.coloring_distance_nm, "a floating point number indicating number of coloring distance in nanometer").default_value(defaultParms.coloring_distance_nm))
        .add_option(Value<int32_t>("-color_num", &parms.color_num, "an integer indicating number of masks (colors) < 3|4 >").required(true))
        .add_option(Value<int32_t>("-simplify_level", &parms.simplify_level, "an integer indicating graph simplification level < 0|1|2 >").default_value(defaultParms.simplify_level))
        .add_option(Value<int32_t>("-thread_num", &parms.thread_num, "an integer indicating maximum thread number").default_value(defaultParms.thread_num))
        .add_option(Value<std::set<int32_t> >("-path_layer", &parms.sPathLayer, "an integer indicating layer for conflict edges"))
        .add_option(Value<std::set<int32_t> >("-precolor_layer", &parms.sPrecolorLayer, "an integer indicating layer for pre-colored patterns"))
        .add_option(Value<std::set<int32_t> >("-uncolor_layer", &parms.sUncolorLayer, "an integer indicating layer for coloring"))
        .add_option(Value<std::string>("-algo", &algo_str, "algorithm type < ILP|BACKTRACK|LP|SDP >").default_value(std::string(defaultParms.algo)))
        .add_option(Value<std::string>("-shape", &shape_str, "shape mode < RECTANGLE|POLYGON >").default_value(std::string(defaultParms.shape_mode)))
        .add_option(Value<bool>("-verbose", &parms.verbose, "toggle controling screen messages").toggle(true).default_value(defaultParms.verbose).toggle_value(true))
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
