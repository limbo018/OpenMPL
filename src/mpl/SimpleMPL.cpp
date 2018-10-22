/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin, Qi Sun
    > Mail: yibolin@utexas.edu, qsun@cse.cuhk.edu.hk
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"

#include <boost/graph/graphviz.hpp>
#include <boost/timer/timer.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>

// only valid when gurobi is available 
#if GUROBI == 1
#include <limbo/algorithms/coloring/ILPColoring.h>
#include <limbo/algorithms/coloring/LPColoring.h>
#include <limbo/algorithms/coloring/MISColoring.h>
#endif
// only valid when lemon cbc api is available 
#if LEMONCBC == 1
#include <limbo/algorithms/coloring/ILPColoringLemonCbc.h>
#endif
#if CSDP == 1
#include <limbo/algorithms/coloring/SDPColoringCsdp.h>
#endif
#include <limbo/algorithms/coloring/BacktrackColoring.h>

#ifdef DEBUG_NONINTEGERS
std::vector<unsigned int> vLP1NonInteger; 
std::vector<unsigned int> vLP1HalfInteger; 
std::vector<unsigned int> vLP2NonInteger; 
std::vector<unsigned int> vLP2HalfInteger; 
std::vector<unsigned int> vLPEndNonInteger; 
std::vector<unsigned int> vLPEndHalfInteger; 
std::vector<unsigned int> vLPNumIter; 
#endif

SIMPLEMPL_BEGIN_NAMESPACE

SimpleMPL::SimpleMPL()
{
    this->reset(true);
}

SimpleMPL::~SimpleMPL()
{
    if (m_db) delete m_db;
}

void SimpleMPL::run(int argc, char** argv)
{
    this->reset(false);
	this->read_cmd(argc, argv);
	this->read_gds();
	this->solve();
	this->report();
	this->write_gds();
}
/*
void SimpleMPL::produce_graph_run(int32_t argc, char **argv)
{
    this->reset(false);
    this->read_cmd(argc, argv);
    this->read_gds();
    std::string vertex_file_name = m_db->output_gds().c_str();
    std::string adjacency_file_name = m_db->output_gds().c_str();
    vertex_file_name.replace(vertex_file_name.end()-4, vertex_file_name.end(),"_vertex.csv");
    adjacency_file_name.replace(adjacency_file_name.end()-4, adjacency_file_name.end(), "_adjacency.csv");
    std::cout<<"vertex_file_name : "<<vertex_file_name<<std::endl;
    std::cout<<"adjacency_file_name : "<<adjacency_file_name<<std::endl;
    this->construct_graph_with_outputs(vertex_file_name, adjacency_file_name);
    this->solve(m_db->output_gds().c_str());

}
*/
void SimpleMPL::reset(bool init)
{
    // release memory and set to initial value 
    if (!init)
    {
        if (m_db) delete m_db;
        std::vector<uint32_t>().swap(m_vVertexOrder);
        std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
        std::vector<uint32_t>().swap(m_vCompId);
        std::vector<uint32_t>().swap(m_vColorDensity);
        std::vector<std::pair<uint32_t, uint32_t> >().swap(m_vConflict);
		std::vector<std::vector<uint32_t> >().swap(SplitMapping);
		std::vector<uint32_t>().swap(new2ori);
    }
    m_db = NULL;
    m_comp_cnt = 0;
}
void SimpleMPL::read_cmd(int argc, char** argv)
{
    // a little bit tricky here 
    // in order to support run-time switch of layoutdb_type
    // we need to construct a dummy ControlParameter for CmdParser
    // then construct actual layoutdb_type according to the option of CmdParser
    // construct a dummy ControlParameter
    ControlParameter tmpParms;
	// read command 
	CmdParser cmd (tmpParms);
	// check options 
    mplAssertMsg(cmd(argc, argv), "failed to parse command");

    // construct actual layout database according to the options 
    if (tmpParms.shape_mode == ShapeModeEnum::RECTANGLE)
        m_db = new LayoutDBRect;
    else 
        m_db = new LayoutDBPolygon;
    // get options from dummy ControlParameter
    tmpParms.swap(m_db->parms);
    std::cout<<"In SimpleMPL::read_cmd() end"<<std::endl;
}

void SimpleMPL::read_gds()
{
    std::cout<<"In SimpleMPL::read_gds() begin"<<std::endl;
    char buf[256];
    mplSPrint(kINFO, buf, "reading input files takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	mplPrint(kINFO, "Reading input file %s\n", m_db->input_gds().c_str());
	// read input gds file 
	GdsReader reader (*m_db);
    mplAssertMsg(reader(m_db->input_gds()), "failed to read %s", m_db->input_gds().c_str());
	// must call initialize after reading 
	m_db->initialize_data();
	// report data 
	m_db->report_data();
    std::cout<<"In SimpleMPL::read_cmd() end"<<std::endl;
}

void SimpleMPL::write_gds()
{
    char buf[256];
    mplSPrint(kINFO, buf, "writing output file takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	if (m_db->output_gds().empty()) 
	{
        mplPrint(kWARN, "Output file not specified, no file generated\n");
		return;
	}
	// write output gds file 
	GdsWriter writer;
	mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
	writer(m_db->output_gds(), *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
}

void SimpleMPL::solve()
{
    // skip if no uncolored layer 
    if (m_db->parms.sUncolorLayer.empty())
        return;
    std::cout<<"In SimpleMPL::solve() begin"<<std::endl;
    char buf[256];
    mplSPrint(kINFO, buf, "coloring takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	if (m_db->vPatternBbox.empty())
	{
        mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}

	this->construct_graph();
    if (m_db->simplify_level() > 0) // only perform connected component when enabled 
        this->connected_component();
    else 
    {
        uint32_t vertex_num = m_vVertexOrder.size();
        uint32_t order_id = 0;
        for (uint32_t v = 0; v != vertex_num; ++v)
        {
            m_vCompId[v] = 0;
            m_vVertexOrder[v] = order_id++;
        }
        m_comp_cnt = 1;
    }

    // create bookmark to index the starting position of each component 
    std::vector<uint32_t> vBookmark (m_comp_cnt);
    for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
    {
        if (i == 0 || m_vCompId[m_vVertexOrder[i-1]] != m_vCompId[m_vVertexOrder[i]])
            vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
    }
#ifdef DEBUG
    std::cout << "============== Before runProjection.. =============="<< std::endl;
    for(uint32_t i = 0; i < vBookmark.size(); i++)
    {
        std::cout << "Component " << i << " starts at : " << vBookmark[i] << std::endl;
    }
#endif

    if(m_db->stitch())
    {
		std::cout << "==== Stitch Insertion ====" << std::endl;
		runProjection(vBookmark);
	}

#ifdef DEBUG
    std::cout << "============== After runProjection.. ==============" << std::endl;
    for(uint32_t i = 0; i < vBookmark.size(); i++)
    {
        std::cout << "Component " << i << " starts at : " << vBookmark[i] << std::endl;
    }
#endif

	mplPrint(kINFO, "Solving %u independent components...\n", m_comp_cnt);
	// thread number controled by user option. Solve the components on parallel.
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
    for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
    {
		std::cout << "solve each component!" << std::endl;
        // construct a component 
        std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+vBookmark[comp_id];
        std::vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+vBookmark[comp_id+1] : m_vVertexOrder.end();
        
        // solve component
        // pass iterators to save memory
        //this->store_component(itBgn, itEnd, comp_id);
        // Here we only want to store the components instead of solving them. If you want to solve them, please uncomment the next line.
        this->solve_component(itBgn, itEnd, comp_id);
	}

#ifdef DEBUG_NONINTEGERS
    mplPrint(kNONE, "vLP1NonInteger vLP1HalfInteger vLP2NonInteger vLP2HalfInteger vLPEndNonInteger vLPEndHalfInteger vLPNumIter\n"); 
    try 
    {
        // I make it simple, so it may go out of range 
        for (uint32_t i = 0; i < vLPNumIter.size(); ++i)
        {
            mplPrint(kNONE, "%u %u %u %u %u %u %u\n", 
                    vLP1NonInteger.at(i), vLP1HalfInteger.at(i), 
                    vLP2NonInteger.at(i), vLP2HalfInteger.at(i), 
                    vLPEndNonInteger.at(i), vLPEndHalfInteger.at(i), 
                    vLPNumIter.at(i)); 
        }
    }
    catch (std::exception const& e)
    {
        mplPrint(kERROR, "%s\n", e.what()); 
    }
    mplPrint(kNONE, "sum of %lu: %u %u %u %u %u %u %u\n", 
            vLPNumIter.size(), 
            limbo::sum(vLP1NonInteger.begin(), vLP1NonInteger.end()), limbo::sum(vLP1HalfInteger.begin(), vLP1HalfInteger.end()), 
            limbo::sum(vLP2NonInteger.begin(), vLP2NonInteger.end()), limbo::sum(vLP2HalfInteger.begin(), vLP2HalfInteger.end()), 
            limbo::sum(vLPEndNonInteger.begin(), vLPEndNonInteger.end()), limbo::sum(vLPEndHalfInteger.begin(), vLPEndHalfInteger.end()), 
            limbo::sum(vLPNumIter.begin(), vLPNumIter.end()));
#endif
    std::cout<<"In SimpleMPL::solve() end"<<std::endl;
}

void SimpleMPL::report() const 
{
    mplPrint(kINFO, "Conflict number = %u\n", conflict_num());
	for (int32_t i = 0, ie = m_db->color_num(); i != ie; ++i)
        mplPrint(kINFO, "Color %d density = %u\n", i, m_vColorDensity[i]);
}

void SimpleMPL::construct_graph_with_outputs(std::string vertex_file_name, std::string adjacency_file_name)
{
    // skip if no uncolored layer
    if (m_db->parms.sUncolorLayer.empty())
        return;

    char buf[256];
    mplSPrint(kINFO, buf, "coloring takes %%t seconds CPU, %%w seconds real\n");
    boost::timer::auto_cpu_timer timer (buf);
    if (m_db->vPatternBbox.empty())
    {
        mplPrint(kWARN, "No patterns found in specified layers\n");
        return;
    }

    mplPrint(kINFO, "Constructing graph for %lu patterns...\n", m_db->vPatternBbox.size());
    // construct vertices
    // assume vertices start from 0 and end with vertex_num - 1
    uint32_t vertex_num = m_db->vPatternBbox.size();
    std::ofstream vertex_out(vertex_file_name.c_str());
    for (uint32_t i=0; i < vertex_num; i++)
    {
        vertex_out << gtl::xl(*(m_db->vPatternBbox[i])) << ","
                   << gtl::yl(*(m_db->vPatternBbox[i])) << ","
                   << gtl::xh(*(m_db->vPatternBbox[i])) << ","
                   << gtl::yh(*(m_db->vPatternBbox[i])) << std::endl;
                   //<< m_db->vPatternBbox[i]->pattern_id() << " "
    }
    vertex_out.close();

    uint32_t edge_num = 0;
    m_vVertexOrder.resize(vertex_num, std::numeric_limits<uint32_t>::max());

    // construct edges
    m_mAdjVertex.resize(vertex_num);
    if(m_db->hPath.empty()) // construct from distance
    {
        std::ofstream adjacency_out(adjacency_file_name.c_str());
        for(uint32_t v=0; v < vertex_num; ++v)
        {
            rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
            std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];

            // find patterns connected with pPattern
            // query tPatternBbox in m_db
            rectangle_type rect(*pPattern);
            // bloat pPattern with minimum coloring distance
            gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
            gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);

            for(rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
                itq != m_db->tPatternBbox.qend(); ++itq)
            {
                rectangle_pointer_type const& pAdjPattern = *itq;
                if(pAdjPattern != pPattern) // skip pPattern itself
                {
                    mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
                    // we consider euclidean distance
                    // use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
                    coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
                    if(distance < m_db->coloring_distance)
                        vAdjVertex.push_back(pAdjPattern->pattern_id());
                }
            }
            vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory
            // parallel with omp reduction here
            uint32_t nei_num = vAdjVertex.size();
            for(uint32_t i = 0; i < nei_num; i++)
            {
                adjacency_out << pPattern->pattern_id() << "," << vAdjVertex[i] << "\n";
            }
            edge_num += vAdjVertex.size();
        }
        adjacency_out.close();
    }
    else
        edge_num = construct_graph_from_paths(vertex_num);
    m_vCompId.resize(vertex_num, std::numeric_limits<uint32_t>::max());
    m_vColorDensity.assign(m_db->color_num(), 0);
    m_vConflict.clear();
    edge_num = edge_num>>1;
    std::cout << "edge_num : " << edge_num << std::endl;
}


void SimpleMPL::construct_graph()
{
    std::cout<<"In SimpleMPL::construct_graph() begin"<<std::endl;
	mplPrint(kINFO, "Constructing graph for %lu patterns...\n", m_db->vPatternBbox.size());
	// construct vertices 
	// assume vertices start from 0 and end with vertex_num-1
	uint32_t vertex_num = m_db->vPatternBbox.size();
#ifdef DEBUG
    // Added by Qi Sun.
    // Used to generate the vertex file.
    std::ofstream vertex_out("/home/qisun/vertex_out.txt");
    for (uint32_t i=0; i < vertex_num; i++)
    {
        vertex_out << m_db->vPatternBbox[i]->pattern_id() << std::endl;
    }
    vertex_out.close();
#endif
	uint32_t edge_num = 0;
	m_vVertexOrder.resize(vertex_num, std::numeric_limits<uint32_t>::max()); 

	// construct edges 
	m_mAdjVertex.resize(vertex_num);
	if (m_db->hPath.empty()) // construct from distance 
        edge_num = construct_graph_from_distance(vertex_num);
	else // construct from conflict edges in hPath
        edge_num = construct_graph_from_paths(vertex_num);
	m_vCompId.resize(vertex_num, std::numeric_limits<uint32_t>::max());
	m_vColorDensity.assign(m_db->color_num(), 0);
	m_vConflict.clear();
    edge_num = edge_num>>1;

	// report statistics 
	mplPrint(kINFO, "%u vertices, %u edges\n", vertex_num, edge_num);
    std::cout<<"In SimpleMPL::construct_graph() end"<<std::endl;
}

uint32_t SimpleMPL::construct_graph_from_distance(uint32_t vertex_num)
{
    std::cout<<"In SimpleMPL::construct_graph_from_distance() begin"<<std::endl;
    uint32_t edge_num = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif
    for (uint32_t v = 0; v < vertex_num; ++v)
    {
        rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];

        // find patterns connected with pPattern 
        // query tPatternBbox in m_db
        rectangle_type rect (*pPattern);
        // bloat pPattern with minimum coloring distance 
        gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
        gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
        for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
                itq != m_db->tPatternBbox.qend(); ++itq)
        {
            rectangle_pointer_type const& pAdjPattern = *itq;
            if (pAdjPattern != pPattern) // skip pPattern itself 
            {
                mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
                // we consider euclidean distance
                // use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
                coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
                if (distance < m_db->coloring_distance)
                    vAdjVertex.push_back(pAdjPattern->pattern_id());
            }
        }
        vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
        // parallel with omp reduction here

        edge_num += vAdjVertex.size();
    }
    std::cout<<"In SimpleMPL::construct_graph_from_distance() end"<<std::endl;
    return edge_num;
}

uint32_t SimpleMPL::construct_graph_from_paths(uint32_t vertex_num)
{
    std::cout<<"In SimpleMPL::construct_graph_from_paths() begin"<<std::endl;
    // at the same time, estimate a coloring distance from conflict edges 
    m_db->coloring_distance = 0;
    for (set<int32_t>::const_iterator itLayer = m_db->parms.sPathLayer.begin(), itLayerE = m_db->parms.sPathLayer.end(); itLayer != itLayerE; ++itLayer)
    {
        if (!m_db->hPath.count(*itLayer)) continue;

        std::vector<path_type> const& vPath = m_db->hPath[*itLayer];

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
        for (uint32_t i = 0; i < vPath.size(); ++i)
        {
            path_type const& path = vPath[i];

            point_type const& p1 = path.low();
            point_type const& p2 = path.high();
            std::list<rectangle_pointer_type> vPattern1;
            std::list<rectangle_pointer_type> vPattern2;
            m_db->tPatternBbox.query(bgi::contains(p1), std::back_inserter(vPattern1));
            m_db->tPatternBbox.query(bgi::contains(p2), std::back_inserter(vPattern2));
#ifdef DEBUG
            if (vPattern1.size() > 1)
            {
                for (std::list<rectangle_pointer_type>::const_iterator it = vPattern1.begin(), ite = vPattern1.end(); it != ite; ++it)
                    mplPrint(kDEBUG, "multiple patterns found for layer %d, (%d, %d, %d, %d)\n", 
                            (*it)->layer(), gtl::xl(**it), gtl::yl(**it), gtl::xh(**it), gtl::yh(**it));
            }
            if (vPattern2.size() > 1)
            {
                for (std::list<rectangle_pointer_type>::const_iterator it = vPattern2.begin(), ite = vPattern2.end(); it != ite; ++it)
                    mplPrint(kDEBUG, "multiple patterns found for layer %d, (%d, %d, %d, %d)\n", 
                            (*it)->layer(), gtl::xl(**it), gtl::yl(**it), gtl::xh(**it), gtl::yh(**it));
            }
#endif
            mplAssert(vPattern1.size() == 1);
            mplAssert(vPattern2.size() == 1);
            rectangle_pointer_type const& pPattern1 = vPattern1.front();
            rectangle_pointer_type const& pPattern2 = vPattern2.front();
            coordinate_difference distance = gtl::length(path);

#ifdef _OPENMP
#pragma omp critical(dataupdate)
#endif
            {
                // if the input gds file contains duplicate edges 
                // we will have duplicate edges here
                m_mAdjVertex[pPattern1->pattern_id()].push_back(pPattern2->pattern_id());
                m_mAdjVertex[pPattern2->pattern_id()].push_back(pPattern1->pattern_id());
                // estimate coloring distance 
                m_db->coloring_distance = std::max(distance, m_db->coloring_distance);
            }
        }
    }
    // eliminate all duplicates in adjacency list 
    uint32_t edge_num = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif 
    for (uint32_t v = 0; v < vertex_num; ++v)
    {
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[v];
        set<uint32_t> sAdjVertex (vAdjVertex.begin(), vAdjVertex.end());
        vAdjVertex.assign(sAdjVertex.begin(), sAdjVertex.end()); 
        vAdjVertex.swap(vAdjVertex); // shrink to fit 
        // parallel with omp reduction
        edge_num += vAdjVertex.size();
    }
    m_db->parms.coloring_distance_nm = m_db->coloring_distance*(m_db->unit*1e+9);
    mplPrint(kINFO, "Estimated coloring distance from conflict edges = %lld (%g nm)\n", m_db->coloring_distance, m_db->coloring_distance_nm());
    std::cout<<"In SimpleMPL::construct_graph_from_paths() end"<<std::endl;
    return edge_num;
}

void SimpleMPL::connected_component()
{
    std::cout<<"In SimpleMPL::connected_component() begin"<<std::endl;
	mplPrint(kINFO, "Computing connected components...\n");
	uint32_t comp_id = 0;
	uint32_t order_id = 0; // position in an ordered pattern array with grouped components 

	uint32_t vertex_num = m_vVertexOrder.size();
	// here we employ the assumption that the graph data structure always has vertices labeled with 0~vertex_num-1
	// m_vVertexOrder only saves an order of it 
    for (uint32_t v = 0; v != vertex_num; ++v)
    {
        if (m_vCompId[v] == std::numeric_limits<uint32_t>::max()) // not visited
        {
            depth_first_search(v, comp_id, order_id);
            comp_id += 1;
        }
    }
	// record maximum number of connected components
    m_comp_cnt = comp_id;

#ifdef DEBUG 
	// check visited 
	for (uint32_t v = 0; v != vertex_num; ++v)
		mplAssert(m_vCompId[v] != std::numeric_limits<uint32_t>::max()); 
#endif

    // reorder with order_id 
    // maybe there is a way to save memory 
    std::vector<uint32_t> vTmpOrder (m_vVertexOrder.size());
    for (uint32_t v = 0; v != vertex_num; ++v)
        vTmpOrder[m_vVertexOrder[v]] = v;
    std::swap(m_vVertexOrder, vTmpOrder);

#ifdef DEBUG
	// check ordered 
	for (uint32_t v = 1; v != vertex_num; ++v)
		mplAssert(m_vCompId[m_vVertexOrder[v-1]] <= m_vCompId[m_vVertexOrder[v]]);
#endif
    std::cout<<"In SimpleMPL::connected_component() end"<<std::endl;
}

void SimpleMPL::depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& order_id)
{
    //std::cout<<"In SimpleMPL::depth_first_search() begin"<<std::endl;
    std::stack<uint32_t> vStack; 
	vStack.push(source);

	while (!vStack.empty())
	{
		uint32_t current = vStack.top();
		vStack.pop();
		if (m_vCompId[current] == std::numeric_limits<uint32_t>::max()) // not visited 
		{
			m_vCompId[current] = comp_id; // set visited 
			m_vVertexOrder[current] = order_id++; // update position 

			for (std::vector<uint32_t>::const_iterator it = m_mAdjVertex[current].begin(); it != m_mAdjVertex[current].end(); ++it)
			{
				uint32_t const& child = *it;
				mplAssert(current != child);
				vStack.push(child);
			}
		}
	}
    //std::cout<<"In SimpleMPL::depth_first_search() end"<<std::endl;
}

///// helper functions for SimpleMPL::solve_component
///// since we may want to call solve_graph_coloring() multiple times 
///// it is better to wrap it as an independent function 

/// create coloring solver pointer according to algorithm type
lac::Coloring<SimpleMPL::graph_type>* SimpleMPL::create_coloring_solver(SimpleMPL::graph_type const& sg) const
{
    typedef lac::Coloring<graph_type> coloring_solver_type;
    coloring_solver_type* pcs = NULL;
    switch (m_db->algo().get())
    {
#if GUROBI == 1
        case AlgorithmTypeEnum::ILP_GURBOI:
            pcs = new lac::ILPColoring<graph_type> (sg); break;
        case AlgorithmTypeEnum::LP_GUROBI:
            pcs = new lac::LPColoring<graph_type> (sg); break;
        case AlgorithmTypeEnum::MIS_GUROBI:
            pcs = new lac::MISColoring<graph_type> (sg); break;
#endif
#if LEMONCBC == 1
        case AlgorithmTypeEnum::ILP_CBC:
            pcs = new lac::ILPColoringLemonCbc<graph_type> (sg); break;
#endif
#if CSDP == 1
        case AlgorithmTypeEnum::SDP_CSDP:
            pcs = new lac::SDPColoringCsdp<graph_type> (sg); break;
#endif
        case AlgorithmTypeEnum::BACKTRACK:
            pcs = new lac::BacktrackColoring<graph_type> (sg);
            break;
        default: mplAssertMsg(0, "unknown algorithm type");
    }
    pcs->stitch_weight(0.1);
    pcs->color_num(m_db->color_num());
    pcs->threads(1); // we use parallel at higher level 

    return pcs;
}

/// given a graph, solve coloring 
/// contain nested call for itself 
uint32_t SimpleMPL::solve_graph_coloring(uint32_t comp_id, SimpleMPL::graph_type const& dg, 
        std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        uint32_t simplify_strategy, std::vector<int8_t>& vColor) const
{
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (dg, m_db->color_num());
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
#ifdef DEBUG
    if (comp_id == 130)
        mplPrint(kDEBUG, "stop\n"); 
#endif
    // set max merge level, actually it only works when MERGE_SUBK4 is on 
    if (m_db->color_num() == 3)
        gs.max_merge_level(3);
    else if (m_db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works better 
        gs.max_merge_level(2);
    gs.simplify(simplify_strategy);
    // gs.write_simplified_graph_by_component(simplified_graph);
    std::cout<<"Simplification Finished."<<std::endl;
	// collect simplified information 
    std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

    // for debug, it does not affect normal run 
	if (comp_id == m_db->dbg_comp_id() && simplify_strategy != graph_simplification_type::MERGE_SUBK4)
		gs.write_simplified_graph_dot("graph_simpl");

	// in order to recover color from articulation points 
	// we have to record all components and mappings 
	// but graph is not necessary 
	std::vector<std::vector<int8_t> > mSubColor (gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig (gs.num_component());
	double acc_obj_value = 0;
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);

		vSubColor.assign(num_vertices(sg), -1);

		// solve coloring 
		typedef lac::Coloring<graph_type> coloring_solver_type;
		coloring_solver_type* pcs = create_coloring_solver(sg);

		// set precolored vertices 
        boost::graph_traits<graph_type>::vertex_iterator vi, vie;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = vColor[vSimpl2Orig[v]];
			if (color >= 0 && color < m_db->color_num())
            {
				pcs->precolor(v, color);
                vSubColor[v] = color; // necessary for 2nd trial
            }
		}
        // 1st trial 
		double obj_value1 = (*pcs)(); // solve coloring 
#ifdef DEBUG
        mplPrint(kDEBUG, "comp_id = %u, %lu vertices, obj_value1 = %g\n", comp_id, num_vertices(sg), obj_value1); 
#endif
        // 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
        double obj_value2 = std::numeric_limits<double>::max();
#ifndef DEBUG_NONINTEGERS
        // very restrict condition to determin whether perform MERGE_SUBK4 or not 
        if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::LP_GUROBI || m_db->algo() == AlgorithmTypeEnum::SDP_CSDP)
                && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
            obj_value2 = solve_graph_coloring(comp_id, sg, itBgn, pattern_cnt, graph_simplification_type::MERGE_SUBK4, vSubColor); // call again
#endif

        // choose smaller objective value 
        if (obj_value1 < obj_value2)
        {
            acc_obj_value += obj_value1;

            // collect coloring results from simplified graph 
            for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
            {
                vertex_descriptor v = *vi;
                int8_t color = pcs->color(v);
                mplAssert(color >= 0 && color < m_db->color_num());
                vSubColor[v] = color;
            }
        }
        else // no need to update vSubColor, as it is already updated by sub call 
            acc_obj_value += obj_value2;
		delete pcs;
	}

	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 
	gs.recover(vColor, mSubColor, mSimpl2Orig);

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
    RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();

    return acc_obj_value;
}

void SimpleMPL::construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt, 
        SimpleMPL::graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor) const
{
	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		vColor[i] = m_db->vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
	}

	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		// conflict edges 
		uint32_t v = *(itBgn+i);
		for (std::vector<uint32_t>::const_iterator it = m_mAdjVertex[v].begin(); it != m_mAdjVertex[v].end(); ++it)
		{
			uint32_t u = *it;
#ifdef DEBUG
			mplAssert(mGlobal2Local.count(u));
#endif
			uint32_t j = mGlobal2Local[u];
			if (i < j) // avoid duplicate 
			{
				std::pair<edge_descriptor, bool> e = edge(i, j, dg);
				if (!e.second) // make sure no duplicate 
				{
					e = add_edge(i, j, dg);
					mplAssert(e.second);
                    boost::put(boost::edge_weight, dg, e.first, 1);
				}
			}
		}
		// stitch edges
		uint32_t parId = new2ori[v];
		std::vector<uint32_t> childVec = SplitMapping[parId];
		for (uint32_t count = 0; count != childVec.size() - 1; count++)
		{
#ifdef DEBUG
			mplAssert(mGlobal2Local.count(childVec[count]));
			mplAssert(mGlobal2Local.count(childVec[count + 1]));
#endif // DEBUG
			uint32_t s = mGlobal2Local[childVec[count]];
			uint32_t t = mGlobal2Local[childVec[count + 1]];
			if (s < t) // avoid duplicate 
			{
				std::pair<edge_descriptor, bool> e = edge(s, t, dg);
				if (!e.second)	// make sure no duplicate
				{
					e = add_edge(s, t, dg);
					mplAssert(e.second);
					// stitch edge has edge_weight -1, different from conflict edges.
					// The coloring solvers will rely on this property.
					boost::put(boost::edge_weight, dg, e.first, -1);
				}
			}
		}
	}
	
}

uint32_t SimpleMPL::solve_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	if (itBgn == itEnd) return 0;
#ifdef DEBUG
    // check order 
	for (std::vector<uint32_t>::const_iterator it = itBgn+1; it != itEnd; ++it)
	{
		uint32_t v1 = *(it-1), v2 = *it;
		mplAssert(m_vCompId[v1] == m_vCompId[v2]);
	}
#endif

    uint32_t acc_obj_value = std::numeric_limits<uint32_t>::max();
    // if current pattern does not contain uncolored patterns, directly calculate conflicts 
    if (check_uncolored(itBgn, itEnd)) 
        acc_obj_value = coloring_component(itBgn, itEnd, comp_id);

	// update global color density map 
	// if parallelization is enabled, there will be uncertainty in the density map 
	// because the density is being updated while being read 
	for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
        int8_t color = m_db->vPatternBbox[*it]->color();
		uint32_t& color_density = m_vColorDensity[color];
#ifdef _OPENMP
#pragma omp atomic
#endif
		color_density += 1;
	}

	uint32_t component_conflict_num = conflict_num(itBgn, itEnd);
    // only valid under no stitch 
    if (acc_obj_value != std::numeric_limits<uint32_t>::max())
        mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, (uint32_t)(itEnd-itBgn), component_conflict_num);

	return component_conflict_num;
}

void SimpleMPL::store_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
    std::cout << "In SimpleMPL::store_component. This is Component "<< comp_id << " !" << std::endl;
    if (itBgn == itEnd)
    {
        std::cout << "This Component " << comp_id << " has only one pattern. We don't need to consider this component!" << std::endl;
        return;
    }
    std::ostringstream ss;
    ss << comp_id;
    std::string id_str = ss.str();
    std::string id_suffix = "_" + id_str;
    std::string file_name = m_db->output_gds().c_str();
    file_name.replace(file_name.end()-4, file_name.end(), id_suffix);
    std::string vertex_file_name = file_name + "_vertex.txt";
    std::string adjacency_file_name = file_name + "_adjacency.txt";

    std::ofstream vertex_out(vertex_file_name.c_str());
    std::ofstream adjacency_out(adjacency_file_name.c_str());
    for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; it++)
    {
        rectangle_pointer_type const& pPattern = m_db->vPatternBbox[*it];
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[*it];
        uint32_t nei_num = vAdjVertex.size();
        adjacency_out<<pPattern->pattern_id()<< " ";
        for (uint32_t i = 0; i < nei_num; i++)
        {
            adjacency_out << vAdjVertex[i] << " ";
        }
        adjacency_out<< std::endl;
        vertex_out << m_db->vPatternBbox[*it]->pattern_id() << " "
                   << gtl::xl(*(m_db->vPatternBbox[*it])) << " "
                   << gtl::yl(*(m_db->vPatternBbox[*it])) << " "
                   << gtl::xh(*(m_db->vPatternBbox[*it])) << " "
                   << gtl::yh(*(m_db->vPatternBbox[*it])) << std::endl;
    }
    vertex_out.close();
    adjacency_out.close();
    std::cout<<" Store Component " << comp_id <<" Successfully!"<<std::endl;
}
/*
void SimpleMPL::store_component_dlx(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	std::cout << "In SimpleMPL::store_component. This is Component "<< comp_id << " !" << std::endl;
    if (itBgn == itEnd)
    {
        std::cout << "This Component " << comp_id << " has only one pattern. We don't need to consider this component!" << std::endl;
        return;
    }
    std::ostringstream ss;
    ss << comp_id;
    std::string id_str = ss.str();
    std::string id_suffix = "_" + id_str;
    std::string file_name = m_db->output_gds().c_str();
    file_name.replace(file_name.end()-4, file_name.end(), id_suffix);
    std::string dlx_file_name = file_name + "_dlx.txt";

    std::ofstream dlx_out(dlx_file_name.c_str());
    int vertex_number = 0;
    int edge_numbers = 0;
    for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; it++)
    {
		vertex_number += 1;
    	//vertex_number += m_db->vPatternBbox[*it].size();
        rectangle_pointer_type const & pPattern = m_db->vPatternBbox[*it];
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[*it];
        edge_numbers += vAdjVertex.size();
    }
    dlx_out << vertex_number << std::endl << edge_numbers << std::endl;
    for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; it++)
    {
    	rectangle_pointer_type const& pPattern = m_db->vPatternBbox[*it];
        std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[*it];
        uint32_t nei_num = vAdjVertex.size();
        for (uint32_t i = 0; i < nei_num; i++)
        {
            dlx_out<<pPattern->pattern_id()<< " " << vAdjVertex[i] << std::endl;
        }
    }
    dlx_out.close();
    std::cout<<" Store Component " << comp_id <<" Successfully!"<<std::endl;
}
*/
uint32_t SimpleMPL::coloring_component(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	// construct a graph for current component 
	uint32_t pattern_cnt = itEnd-itBgn;

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...\n", comp_id, pattern_cnt);

	// decomposition graph 
    // must allocate memory here 
	graph_type dg (pattern_cnt);
	std::vector<int8_t> vColor (pattern_cnt, -1); // coloring results 
	map<uint32_t, uint32_t> mGlobal2Local; // global vertex id to local vertex id 

    // construct decomposition graph for component 
    construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor);

    // for debug, it does not affect normal run 
	if (comp_id == m_db->dbg_comp_id())
        write_graph(dg, "graph_init");

	// graph simplification 
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
    uint32_t simplify_strategy = graph_simplification_type::NONE;
	// keep the order of simplification 
	if (m_db->simplify_level() > 1)
        simplify_strategy |= graph_simplification_type::HIDE_SMALL_DEGREE;
	if (m_db->simplify_level() > 2)
        simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

    // solve graph coloring 
    uint32_t acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor);

#ifdef DEBUG
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		mplAssert(vColor[i] >= 0 && vColor[i] < m_db->color_num());
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
        m_db->set_color(v, vColor[i]);
	}

    return acc_obj_value;
}

uint32_t SimpleMPL::conflict_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd) const
{
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	uint32_t cnt = 0;
    for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
		uint32_t v = *it;
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				int8_t color2 = vPatternBbox[u]->color();
				if (color2 >= 0 && color2 < m_db->color_num())
				{
					if (color1 == color2) ++cnt;
				}
				else ++cnt; // uncolored vertex is counted as conflict 
			}
		}
		else ++cnt; // uncolored vertex is counted as conflict 
	}
	// conflicts will be counted twice 
	return (cnt>>1);
}

uint32_t SimpleMPL::conflict_num() const
{
	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 < m_db->color_num())
		{
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
                if (v < u) // avoid duplicate 
                {
                    int8_t color2 = vPatternBbox[u]->color();
                    if (color2 >= 0 && color2 < m_db->color_num())
                    {
                        if (color1 == color2) 
                            m_vConflict.push_back(std::make_pair(v, u));
                    }
                    else // uncolored vertex is counted as conflict 
                        mplAssertMsg(0, "uncolored vertex %u = %d", u, color2);
                }
			}
		}
        else // uncolored vertex is counted as conflict 
            mplAssertMsg(0, "uncolored vertex %u = %d", v, color1);
	}
	return m_vConflict.size();
}

bool SimpleMPL::check_uncolored(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd) const
{
    for (; itBgn != itEnd; ++itBgn)
        if (m_db->vPatternBbox[*itBgn]->color() < 0)
            return true;
    return false;
}

void SimpleMPL::write_graph(SimpleMPL::graph_type& g, std::string const& filename) const
{
    // in order to make the .gv file readable by boost graphviz reader 
    // I dump it with boost graphviz writer 
    boost::dynamic_properties dp;
    dp.property("id", boost::get(boost::vertex_index, g));
    dp.property("node_id", boost::get(boost::vertex_index, g));
    dp.property("label", boost::get(boost::vertex_index, g));
    // somehow edge properties need mutable graph_type& 
    dp.property("weight", boost::get(boost::edge_weight, g));
    dp.property("label", boost::get(boost::edge_weight, g));
    std::ofstream out ((filename+".gv").c_str());
    boost::write_graphviz_dp(out, g, dp, string("id"));
    out.close();
    la::graphviz2pdf(filename);
}

bool SimpleMPL::whetherHorizontal(rectangle_pointer_type tmp)
{
    double xl = gtl::xl(*tmp);
    double yl = gtl::yl(*tmp);
    double xh = gtl::xh(*tmp);
    double yh = gtl::yh(*tmp);
    return (xh - xl) > (yh - yl);
}

void SimpleMPL::GenerateStitchPosition(const rectangle_pointer_type pRect,
	const std::vector<rectangle_type> vinterRect,
	std::vector <coordinate_type> vstitches, const coordinate_type lower,
	const coordinate_type upper)
{
	// ================================================================================
	// step 1 : generate candidate stitches' positions according to the intersections
	// ================================================================================
	bool isHor = whetherHorizontal(pRect);
	std::set<coordinate_type> tempSet;
	tempSet.insert(lower);
	tempSet.insert(upper);
	for (std::vector<rectangle_type>::const_iterator it = vinterRect.begin(); it != vinterRect.end(); it++)
	{
		if (isHor)
		{
			tempSet.insert(gtl::xl(*it));
			tempSet.insert(gtl::xh(*it));
		}
		else
		{
			tempSet.insert(gtl::yl(*it));
			tempSet.insert(gtl::yh(*it));
		}
	}
	std::vector<coordinate_type> tempVec;
	for(std::set<coordinate_type>::iterator itr = tempSet.begin(); itr != tempSet.end(); itr++)
		tempVec.push_back(*itr);
	// sort all the positions
	sort(tempVec.begin(), tempVec.end());
	
	// ================================================================================
	// step 2 : generate stages according to the stitch positions
	//			that is each pairwise of neighboring positions generates a stage.
	// ================================================================================
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> > vStages;
	for (uint32_t i = 1; i < tempVec.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(tempVec[i - 1], tempVec[i]), 0));
	// calculate the times every stage covered by all the intersections.
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
	{
		coordinate_type left = it->first.first;
		coordinate_type right = it->first.second;
		int overlapping_count = 0;
		for (std::vector<rectangle_type>::const_iterator itInt = vinterRect.begin(); itInt != vinterRect.end(); itInt++)
		{
			if (isHor)
			{
				if (left < gtl::xl(*itInt)) continue;
				if (right > gtl::yh(*itInt)) continue;
				overlapping_count++;
			}
			else
			{
				if (left < gtl::yl(*itInt)) continue;
				if (right > gtl::yh(*itInt)) continue;
				overlapping_count++;
			}
		}
		it->first.second = overlapping_count;
	}

	// ================================================================================
	// step 3 : add default terminal zeros
	//			This will be used in the next step.
	// ================================================================================
	if (vStages.size() <= 0) return;
	if (vStages[0].second != 0)
		vStages.insert(vStages.begin(), std::make_pair(std::make_pair(lower, lower), 0));
	if (vStages[vStages.size() - 1].second != 0)
		vStages.push_back(std::make_pair(std::make_pair(upper, upper), 0));
#ifdef DEBUG
	std::cout << "DEBUG_PROJECTION| vStages = ";
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
		std::cout << it->second;
	std::cout << std::endl;
#endif

	// ================================================================================
	// step 4: find the stages with zero overlapping_count
	//		   The stitches will be chosen from these stages.
	// ================================================================================
	std::vector<uint32_t> vZeroIds;
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		if (vStages[i].second > 0) continue;
		vZeroIds.push_back(i);
	}

	// ================================================================================
	// step 5: choose stitches from vZeroIds
	// ================================================================================
	std::vector<coordinate_type>().swap(vstitches);
	// The operations here are very confusing.
	for (uint32_t i = 0; i < vZeroIds.size() - 1; i++)
	{
		uint32_t pos1 = vZeroIds[i];
		uint32_t pos2 = vZeroIds[i + 1];
		// since ((lower, lower), 0) has been added into vStages, so pos1 must be 0.
		if (i == 0) mplAssertMsg(0 == pos1, "pos1 %d doesn't equal to 0", pos1);
		// remove the useless stitches
		else if (pos1 == 2)
		{
			bool find = false;
			if (vStages.size() < 5) find = true;
			else if (i != 1) find = true;
			else if (1 != vStages[1].second) find = true;
			else if (1 != vStages[3].second) find = true;
			else if (0 != vStages[4].second) find = true;
			if (find == false)
				continue;
			coordinate_type position = (vStages[2].first.first + vStages[2].first.second) / 2;
			vstitches.push_back(position);
		}
		else if (pos1 == vStages.size() - 3)
		{
			mplAssert(i == vZeroIds.size() - 2);
			bool find = false;
			uint32_t zsize = vZeroIds.size();
			if (vStages.size() < 5) find = true;
			else if (i != zsize - 2) find = true;
			else if (1 != vStages[pos1 + 1].second) find = true;
			else if (1 != vStages[pos1 - 1].second) find = true;
			else if (0 != vStages[pos1 - 2].second) find = true;
			if (find == false) continue;
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}
		else
		{
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}

		// search lost stitch in vStages[pos1 --> pos2]
		double maxValue = 0.9;
		uint32_t posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			uint32_t mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
			uint32_t diff = std::abs(static_cast<int>(vStages[i + 1].second - vStages[i - 1].second));
			double value = (double)mind + (double)diff * 0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			uint32_t position = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
			vstitches.push_back(position);
		}
	}
	sort(vstitches.begin(), vstitches.end());
#ifdef DEBUG
	std::cout << "DEBUG| output the vStages: " << std::endl;
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
		std::cout << it->second << "[" << it->first.first << ", " << it->first.second << "]" << std::endl;
	std::cout << std::endl;

	std::cout << "DEBUG| output the vstitches: " << std::endl;
	std::cout << "lower = " << lower << std::endl;
	std::cout << "upper = " << upper << std::endl;
	for (std::vector<coordinate_type>::iterator it = vstitches.begin(); it != vstitches.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
#endif
	return;
}

void SimpleMPL::runProjection(std::vector<uint32_t> & vBookmark)
{
	// In order to execute the projection operations in parallel, I think it's needed to take up more 
	// spaces to store the intermediate results.

	// ===========================================================================================
	// step 1 : run projections in parallel, I think it's needed to take up more 
	// spaces to store the intermediate results.
	//			The newly-generated patterns are stored in Component_Pattern_list, 
	//			with the index of pattern ids.
	// ===========================================================================================
	
	// resize SplitMapping to store all patterns in original graph
	SplitMapping.resize(m_vVertexOrder.size());
	std::vector<std::vector<rectangle_type> > Component_Pattern_list;
    
//#ifdef _OPENMP
//#pragma omp parallel for private(comp_id) num_threads(m_db->thread_num())
//#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		
        std::vector<rectangle_type> new_PatternBox_temp;
		projection(itBgn, itEnd, new_PatternBox_temp);
		Component_Pattern_list.push_back(new_PatternBox_temp);
	}
	// ================================================================================================
	// step 2 : update the graph information, generate ids for new patterns.
	// I don't know how to implement these in parallel. It seems 'omp critical' is not useful.
	// ================================================================================================

	// generate new vBookmark
    std::vector<uint32_t>().swap(vBookmark);
    uint32_t total_pattern_number = 0;

    for(uint32_t comp_id = 0; comp_id < m_comp_cnt; comp_id++)
    {
        vBookmark.push_back(total_pattern_number);
        total_pattern_number += Component_Pattern_list[comp_id].size();
    }

    // push total_pattern_number  into vBookmark to help generate new graph information
    vBookmark.push_back(total_pattern_number);

    std::vector<rectangle_pointer_type>().swap(m_db->vPatternBbox);
    std::vector<uint32_t>().swap(m_vVertexOrder);
    std::vector<uint32_t>().swap(m_vCompId);
    std::vector<uint32_t>().swap(new2ori);

    uint32_t pattern_number = 0;
    // used as pivot 
    uint32_t comp_id = 0;
	for (uint32_t itVec = 0; itVec < SplitMapping.size(); itVec++)
	{
		for (std::vector<uint32_t>::iterator itsplit = SplitMapping[itVec].begin(); itsplit != SplitMapping[itVec].end(); itsplit++)
		{
			Component_Pattern_list[comp_id][*itsplit].pattern_id(pattern_number);
			m_db->vPatternBbox.push_back(&Component_Pattern_list[comp_id][*itsplit]);
			// update the mapping relationship in SplitMapping, which will be used in step 3.
			*itsplit = pattern_number;
			// map new patterns back to original patterns 
			new2ori.push_back(itVec);
			m_vVertexOrder.push_back(pattern_number);
            if(pattern_number >= vBookmark[comp_id + 1]) comp_id++;
            m_vCompId.push_back(comp_id);
			pattern_number++;
		}
	}
    // pop total_pattern_number added just now from vBookmark
    vBookmark.pop_back();

	// ==============================================================================================
	// step 3 : compute the adjacency matrix.
	//			Now only consider the conflicts between 
	//			For the patterns in one original pattern, all the relationships will be measured 
	//			during the component solving process.
	// ==============================================================================================
	std::vector<std::vector<uint32_t> > new_mAdjVertex;
	adj4NewPatterns(new_mAdjVertex);
	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex = new_mAdjVertex;
	return;
}

// Use each component as the input, the output is the of the new component. 
// The relevant information should also be modified.
// itBgn : component begin node in the vector
// itEnd : component end node in the vector
void SimpleMPL::projection(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, std::vector<rectangle_type>& new_PatternVec)
{
	// ===================================================================
	// traverse all the vertices in the layout graph
	// ===================================================================
    // pattern_count : used to regenerate pattern ids for all patterns in this component
    uint32_t pattern_count = 0;
	for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; it++)
	{
		// ===================================================================
		// step 1 : capture the interaction parts with its neighbors.
		// ===================================================================
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[*it];
		bool isHor = whetherHorizontal(pPattern);
		// vInterRect stores the intersection parts of pPattern and its neighbors.
		std::vector<rectangle_type> vInterRect;
		vInterRect.clear();
		std::vector<uint32_t>& vAdjVertex = m_mAdjVertex[m_db->vPatternBbox[*it]->pattern_id()];
        vInterRect.resize(vAdjVertex.size());
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
		for (uint32_t nei = 0; nei < vAdjVertex.size(); nei ++)
		{
			rectangle_type extendPattern(*m_db->vPatternBbox[vAdjVertex[nei]]);
			gtl::bloat(extendPattern, gtl::HORIZONTAL, m_db->coloring_distance);
			gtl::bloat(extendPattern, gtl::VERTICAL, m_db->coloring_distance);
			
#ifdef BOOST_REG_INTERSECTION
			vInterRect[nei] = interSectionRectBoost(extendPattern, *pPattern);
#else
			vInterRect[nei] = interSectionRect(extendPattern, *pPattern);
#endif
		} // for nei, Traverse all the intersections with its neighbors.

		  // ===================================================================
		  // step 2 : generate all the candidate stitches for this pattern.
		  // ===================================================================
		std::vector<coordinate_type> vstitches;
		rectangle_pointer_type tempRect = new rectangle_type(gtl::xl(*pPattern), gtl::yl(*pPattern), gtl::xh(*pPattern), gtl::yh(*pPattern));
		// the lower bound and upper bound of the stitch position.
		// If the pattern is horizontal, the stitches will be horizontal.
		// If the pattern is vertical, the stitches will be vertical.
		coordinate_type lower = 0, upper = 0;
		if (isHor) {
			lower = gtl::xl(*tempRect);
			upper = gtl::xh(*tempRect);
		}
		else
		{
			lower = gtl::yl(*tempRect);
			upper = gtl::yh(*tempRect);
		}
		// Generate stitch points, based on Bei Yu's method.
		GenerateStitchPosition(tempRect, vInterRect, vstitches, lower, upper);

#ifdef DEBUG
		// ===================================================================
		// step 3 : check the positions' legalities
		// ===================================================================
		if (isHor) {
			// If pPattern is horizontal, all the stitches' positions should be in (xl, xh)
			for (uint32_t j = 0; j < vstitches.size(); j++)
			{
				int pos = vstitches[j];
				assert(pos > gtl::xl(*pPattern) && pos < gtl::xh(*pPattern));
			}
		}
		else {
			// If pPattern is vertical, all the stitches' position should be in (yl, yh)
			for (uint32_t j = 0; j < vstitches.size(); j++)
			{
				int pos = vstitches[j];
				assert(pos > gtl::yl(*pPattern) && pos < gtl::yh(*pPattern));
			}
		}
#endif
		// ===============================================================================================
		// step 4 : split the patterns according to the stitches
		//			new_PatternVec	: stores the newly-generated patterns
		//			SplitMapping	: mapping relationships between original patterns and split patterns
		// ===============================================================================================
		std::vector<rectangle_type>().swap(new_PatternVec);
        // If this pattern hasn't been split
		if (vstitches.size() <= 0)
		{
            // shouldn't change pPattern, because vPatternBbox will be used in the following steps.
			rectangle_pointer_type new_Pattern = new rectangle_type(*pPattern);
			new_Pattern->pattern_id(pattern_count);
			new_PatternVec.push_back(*new_Pattern);
			SplitMapping[pPattern->pattern_id()].push_back(pattern_count);
			pattern_count++;
		}
		// This pattern has been splited.
		else
		{
			if (isHor)
			{
				// vstitches : position order, from left to right
				// If horizontal, insert xl and xh into vstitches.
				vstitches.insert(vstitches.begin(), gtl::xl(*pPattern));
				vstitches.push_back(gtl::xh(*pPattern));
				for (uint32_t j = 0; j < vstitches.size() - 1; j++)
				{
					rectangle_pointer_type new_Pattern = new rectangle_type;
					gtl::xl(*new_Pattern, vstitches[j]);
					gtl::yl(*new_Pattern, gtl::yl(*pPattern));
					gtl::xh(*new_Pattern, vstitches[j + 1]);
					gtl::yh(*new_Pattern, gtl::yh(*pPattern));
					new_Pattern->pattern_id(pattern_count);
					// for precolored patterns
					new_Pattern->color(pPattern->color());
					new_PatternVec.push_back(*new_Pattern);
					SplitMapping[pPattern->pattern_id()].push_back(pattern_count);
					pattern_count++;
				}
			}
			else
			{
				// vstitches : position order, from bottom to up
				// If vertical, insert yl and yh into vstitches.
				vstitches.insert(vstitches.begin(), gtl::yl(*pPattern));
				vstitches.push_back(gtl::yh(*pPattern));
				for (uint32_t j = 0; j < vstitches.size() - 1; j++)
				{
					rectangle_pointer_type new_Pattern = new rectangle_type();
					gtl::xl(*new_Pattern, gtl::xl(*pPattern));
					gtl::yl(*new_Pattern, vstitches[j]);
					gtl::xh(*new_Pattern, gtl::xh(*pPattern));
					gtl::yh(*new_Pattern, vstitches[j + 1]);
					new_Pattern->pattern_id(pattern_count);
					// for precolored patterns
					new_Pattern->color(pPattern->color());
					new_PatternVec.push_back(*new_Pattern);
					SplitMapping[pPattern->pattern_id()].push_back(pattern_count);
					pattern_count++;
				}
			}
		}
		mplAssert(new_PatternVec.size() > 0);
	}
	return;
}

void SimpleMPL::adj4NewPatterns(std::vector<std::vector<uint32_t> > & new_mAdjVertex)
{
//#ifdef _OPENMP
//#pragma omp parallel for 
//#endif
	for (uint32_t newPattern = 0; newPattern < new2ori.size(); newPattern++)
	{
		uint32_t parentId = new2ori[newPattern];
		// traverse parent's neighbor list
		for (uint32_t parNei = 0; parNei < m_mAdjVertex[parentId].size(); parNei++)
		{
			uint32_t parNeiId = m_mAdjVertex[parentId][parNei];
			std::vector<uint32_t> possNeiVec = SplitMapping[parNeiId];
			for (std::vector<uint32_t>::iterator it = possNeiVec.begin(); it != possNeiVec.end(); it++)
			{
				coordinate_difference distance = m_db->euclidean_distance(*m_db->vPatternBbox[*it], *m_db->vPatternBbox[newPattern]);
				if (distance < m_db->coloring_distance)
//#ifdef _OPENMP
//#pragma omp critical
//#endif
                {
                    new_mAdjVertex[newPattern].push_back(*it);
				}
			}
		}
	}
	return;
}

LayoutDB::rectangle_type SimpleMPL::interSectionRect(rectangle_type rect1,  rectangle_type rect2)
{
	coordinate_type xl = std::max(gtl::xl(rect1), gtl::xl(rect2));
	coordinate_type yl = std::max(gtl::yl(rect1), gtl::yl(rect2));
	coordinate_type xh = std::min(gtl::xh(rect1), gtl::xh(rect2));
	coordinate_type yh = std::min(gtl::yh(rect1), gtl::yh(rect2));
	rectangle_pointer_type output = new rectangle_type(xl, yl, xh, yh);
	return *output;
}



void SimpleMPL::print_welcome() const
{
  mplPrint(kNONE, "\n\n");
  mplPrint(kNONE, "=======================================================================\n");
  mplPrint(kNONE, "                        OpenMPL - Version 1.1                        \n");
  mplPrint(kNONE, "                                by                                   \n");  
  mplPrint(kNONE, "                Yibo Lin, Bei Yu, Qi Sun and  David Z. Pan           \n");
  mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
  mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
  mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
  mplPrint(kNONE, "            Contact Authors:  {yibolin, dpan}@cerc.utexas.edu        \n");
  mplPrint(kNONE, "                              {byu, qsun}@cse.cuhk.edu.hk            \n");
  mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE
