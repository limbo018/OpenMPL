/*************************************************************************
    > File Name: io.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 06 Nov 2014 08:53:46 AM CST
 ************************************************************************/

#ifndef SIMPLEMPL_GDSIIIO_H
#define SIMPLEMPL_GDSIIIO_H

#include <fstream>
#include <limits>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <limbo/parsers/gdsii/stream/GdsWriter.h>

#include "LayoutDB.h"

SIMPLEMPL_BEGIN_NAMESPACE

namespace gtl = boost::polygon;
using boost::int32_t;
using boost::int64_t;
using boost::array;
using gtl::point_concept;
using gtl::rectangle_concept;
using gtl::polygon_90_concept;
using gtl::polygon_90_set_concept;
using gtl::point_data;
using gtl::rectangle_data;
using gtl::polygon_90_data;
using gtl::polygon_90_set_data;

using namespace gtl::operators;

/// read gds file 
struct GdsReader : GdsParser::GdsDataBaseKernel
{
	typedef LayoutDB layoutdb_type;
	typedef layoutdb_type::coordinate_type coordinate_type;
	typedef layoutdb_type::point_type              point_type;
	typedef layoutdb_type::rectangle_type          rectangle_type;
	typedef layoutdb_type::polygon_type            polygon_type;
	typedef layoutdb_type::polygon_pointer_type    polygon_pointer_type;
	typedef layoutdb_type::rectangle_pointer_type  rectangle_pointer_type;
	typedef layoutdb_type::path_type               path_type;

	//string strname; // TOPCELL name, useful for dump out gds files 
	//double unit;
	int32_t layer;
    GdsParser::GdsRecords::EnumType status; 
    std::vector<point_type> vPoint;
	int64_t file_size; // in bytes 

	layoutdb_type& db;  

	GdsReader(layoutdb_type& _db) : db(_db) {}

    /// top api 
	bool operator() (std::string const& filename); 

	// required callbacks in parser 
	virtual void bit_array_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vBitArray);
	virtual void integer_2_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vInteger);
	virtual void integer_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vInteger);
	virtual void real_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vFloat);
	virtual void real_8_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vFloat);
	virtual void string_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::string const& str);
	virtual void begin_end_cbk(GdsParser::GdsRecords::EnumType record_type);

    /// helper functions 
    void integer_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<int> const& vData);
    void float_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, std::vector<double> const& vData);
};

/// write gds file 
struct GdsWriter
{
	typedef LayoutDB layoutdb_type;
    typedef layoutdb_type::coordinate_type coordinate_type;
	typedef layoutdb_type::point_type point_type;
	typedef layoutdb_type::rectangle_type rectangle_type;
	typedef layoutdb_type::polygon_type polygon_type;
	typedef layoutdb_type::polygon_pointer_type polygon_pointer_type;
	typedef layoutdb_type::rectangle_pointer_type rectangle_pointer_type;
	typedef layoutdb_type::path_type path_type;

	void operator() (std::string const& filename, layoutdb_type const& db, 
			std::vector<std::pair<uint32_t, uint32_t> > const& vConflict, 
			std::vector<std::vector<uint32_t> > const& mAdjVertex, 
			std::string const& strname = "TOPCELL", double unit = 0.001) const;
    /// write rectangles 
	void write_rectangles(GdsParser::GdsWriter& gw, std::vector<rectangle_pointer_type> const& vRect, const int32_t layer_offset) const;
    /// write conflicts 
	void write_conflicts(GdsParser::GdsWriter& gw, layoutdb_type const& db, 
			std::vector<std::pair<uint32_t, uint32_t> > const& vConflict, const int32_t layer) const;
    /// write paths before constructing conflict edges 
	void write_paths(GdsParser::GdsWriter& gw, std::map<int32_t, std::vector<path_type> > const& hPath) const;
    /// write conflict edges 
	void write_edges(GdsParser::GdsWriter& gw, layoutdb_type const& db, std::vector<std::vector<uint32_t> > const& mAdjVertex, const int32_t layer) const; 

};

SIMPLEMPL_END_NAMESPACE

#endif 
