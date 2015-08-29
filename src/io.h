/*************************************************************************
    > File Name: io.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Thu 06 Nov 2014 08:53:46 AM CST
 ************************************************************************/

#ifndef SIMPLEMPL_IO_H
#define SIMPLEMPL_IO_H

#include <fstream>
#include <limits>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <limbo/parsers/gdsii/stream/GdsWriter.h>
#include <limbo/programoptions/ProgramOptions.h>

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
struct GdsReader : GdsParser::GdsDataBase
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
	int32_t status; // 0: not in any block, 1 in BOUNDARY or BOX block, 2 in PATH   
    std::vector<point_type> vPoint;
	int64_t file_size; // in bytes 

	layoutdb_type& db;  

	GdsReader(layoutdb_type& _db) : db(_db) {}

	bool operator() (std::string const& filename)  
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

	template <typename ContainerType>
	void general_cbk(std::string const& ascii_record_type, std::string const& ascii_data_type, ContainerType const& vData)
	{
		if (ascii_record_type == "UNITS")
		{
			db.unit = vData[1]; 
		}
		else if (ascii_record_type == "BOUNDARY" || ascii_record_type == "BOX")
		{
			vPoint.clear();
			layer = 0;
			status = 1;
		}
		else if (ascii_record_type == "PATH")
		{
			vPoint.clear();
			layer = 0;
			status = 2;
		}
		else if (ascii_record_type == "LAYER")
		{
			layer = vData[0];
		}
		else if (ascii_record_type == "XY")
		{
			if (status == 1 || status == 2)
			{
				mplAssert((vData.size() & 1) == 0 && vData.size() >= 4);
				vPoint.clear();
				uint32_t end = vData.size();
				// skip last point for BOX and BOUNDARY
				if (status == 1) end -= 2;
				for (uint32_t i = 0; i < end; i += 2)
					vPoint.push_back(gtl::construct<point_type>(vData[i], vData[i+1]));
			}
		}
		else if (ascii_record_type == "ENDEL")
		{
			if (status == 1 || status == 2)
			{
				mplAssert(layer != -1);

				db.add(layer, vPoint);

				status = 0;
			}
		}
		else if (ascii_record_type == "STRNAME")
		{
			mplAssert(ascii_data_type == "STRING");
			mplAssert(!vData.empty());
			db.strname.assign(vData.begin(), vData.end());
		}
	}

	// required callbacks in parser 
	virtual void bit_array_cbk(const char* ascii_record_type, const char* ascii_data_type, std::vector<int> const& vBitArray)
	{this->general_cbk(ascii_record_type, ascii_data_type, vBitArray);}
	virtual void integer_2_cbk(const char* ascii_record_type, const char* ascii_data_type, std::vector<int> const& vInteger)
	{this->general_cbk(ascii_record_type, ascii_data_type, vInteger);}
	virtual void integer_4_cbk(const char* ascii_record_type, const char* ascii_data_type, std::vector<int> const& vInteger)
	{this->general_cbk(ascii_record_type, ascii_data_type, vInteger);}
	virtual void real_4_cbk(const char* ascii_record_type, const char* ascii_data_type, std::vector<double> const& vFloat) 
	{this->general_cbk(ascii_record_type, ascii_data_type, vFloat);}
	virtual void real_8_cbk(const char* ascii_record_type, const char* ascii_data_type, std::vector<double> const& vFloat) 
	{this->general_cbk(ascii_record_type, ascii_data_type, vFloat);}
	virtual void string_cbk(const char* ascii_record_type, const char* ascii_data_type, std::string const& str) 
	{this->general_cbk(ascii_record_type, ascii_data_type, str);}
	virtual void begin_end_cbk(const char* ascii_record_type)
	{this->general_cbk(ascii_record_type, "", std::vector<int>());}

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

/// parse command line arguments 
struct CmdParser
{
	ControlParameter& parms;

	CmdParser(ControlParameter& p) : parms(p) {}

	bool operator()(int argc, char** argv);
};

SIMPLEMPL_END_NAMESPACE

#endif 
