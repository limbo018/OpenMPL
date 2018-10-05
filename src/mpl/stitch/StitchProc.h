/*************************************************************************
> File Name: StitchProc.h
> Author: Qi Sun
> Mail: qsun@cse.cuhk.edu.hk
> Created Time: 27/08/2018
************************************************************************/


#ifndef _STITCHPROC_H_
#define _STITCHPROC_H_

#include <iostream>
#include <stack>
#include <deque>
#include <limbo/geometry/api/GeoBoostPolygonApi.h>
#include <limbo/geometry/Geometry.h>
#include "GdsiiIO.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include <stack>
#include <boost/graph/graphviz.hpp>
#include <boost/timer/timer.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>
#include <boost/polygon/interval_data.hpp>


SIMPLEMPL_BEGIN_NAMESPACE

class StitchProc 
{
public:
	typedef LayoutDB layoutdb_type;
	typedef layoutdb_type::coordinate_type coordinate_type;
	typedef layoutdb_type::coordinate_difference   coordinate_difference;
	typedef layoutdb_type::point_type              point_type;
	typedef layoutdb_type::rectangle_type          rectangle_type;
	typedef layoutdb_type::polygon_type            polygon_type;
	typedef layoutdb_type::polygon_pointer_type    polygon_pointer_type;
	typedef layoutdb_type::rectangle_pointer_type  rectangle_pointer_type;
	typedef layoutdb_type::path_type               path_type;
	typedef layoutdb_type::rtree_type              rtree_type;
	typedef layoutdb_type::graph_type              graph_type;
	typedef layoutdb_type::vertex_descriptor       vertex_descriptor;
	typedef layoutdb_type::edge_descriptor         edge_descriptor;

	// default constructor
	StitchProc();
	StitchProc(double pitch);
	// destructor
	~StitchProc();

	void read_cmd(int argc, char** argv);
	void read_gds();
	void write_gds();
	void relation4rect();
	// void construct_graph();
	void projection();
    

protected:
	// reset data members
	void reset(bool init);
	coordinate_type getWidth(rectangle_pointer_type rect) { return gtl::xh(*rect) - gtl::xl(*rect); }
	coordinate_type getHeight(rectangle_pointer_type rect) { return gtl::yh(*rect) - gtl::yl(*rect); }
	
	bool whetherHor(rectangle_pointer_type tmp);
	void StitchGenerateDPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, const coordinate_type left, const coordinate_type right, bool vddstitch = false);
	void StitchGenerateTPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, const coordinate_type left, const coordinate_type right, bool vddstitch = false);

	// pointer to layout database and user-defined options
	layoutdb_type *m_db;
	// adjacency list data structure for a graph
	std::vector<uint32_t>				m_vVertexOrder;
	// adjacency list
	// std::vector<std::vector<uint32_t> > m_mAdjVertex;
	// independent component id
	std::vector<uint32_t>			m_vCompId;
	// max of connected components
	uint32_t						m_comp_cnt;
	std::vector<rectangle_type>		m_final_rectangle;
	// conflict report 
	mutable std::vector<std::pair<uint32_t, uint32_t> > m_vConflict;
    
    std::vector<std::vector<uint32_t>> m_Touch;   // This vector is for polygon type data.
    std::vector<std::vector<uint32_t>> m_DPL;     // This vector stores the conflict patterns.
												  
	std::vector<std::vector<uint32_t>> m_Safe;    // This vector stores the safe patterns.
	std::vector<std::vector<rectangle_type>> m_interRect;

	double MINDPLDIST			= 0;		// minimual DPL distance
	double PITCH				= 0;		// still not sure whether use pitch in this problem 
	bool dplstitch				= true;		// whether use DPL method to generate stitch candidates, if false, then use TPL method.

	coordinate_type m_layout_left		= INT32_MAX;
	coordinate_type m_layout_right		= INT32_MIN;
	coordinate_type m_layout_top		= INT32_MIN;
	coordinate_type m_layout_down		= INT32_MAX;

};


SIMPLEMPL_END_NAMESPACE


#endif // !_STITCHPROC_H_