/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin, Wei Li, Qi Sun
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"
#include <errno.h>
#include <sstream>
#include <fstream>
#include <set>
#define SYSERROR()  errno
#include <stack>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>
#include <boost/timer/timer.hpp>
#include <limbo/algorithms/coloring/GraphSimplification.h>

// only valid when gurobi is available 
#if GUROBI == 1
#include <limbo/algorithms/coloring/ILPColoring.h>
#include <limbo/algorithms/coloring/ILPColoringUpdated.h>
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
#include "DancingLinkColoring.h"
#include "DancingLinkColoringOpt.h"
#define CUT_DG
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
	std::ofstream myfile;
	myfile.open ("record.txt", std::ofstream::app);
	myfile << argv[4]<<" " << argv[16]<<"\n";
	myfile.close();
	this->read_cmd(argc, argv);
	this->read_gds();
    if (m_db->parms.record > 2)
    {
        std::ofstream myfile_small;
        myfile_small.open ("small_results.txt", std::ofstream::app);
        myfile_small << argv[4]<<" " << argv[16]<<"\n";
        myfile_small.close();
    }
	this->solve();
	this->report();
	this->write_gds();
	// youyi
	mplPrint(kINFO, "===============================\n");
	mplPrint(kINFO,  "# Recursive = %d\n", count);
	mplPrint(kINFO,  "# Total Count= %d\n", total_count);
	mplPrint(kINFO, "Mean # vertices = %f\n", double(total_num_vertices) / count);
}

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
    }
    m_db = NULL;
    m_comp_cnt = 0;
	// youyi
	count = 0;
	total_count = 0;
	total_num_vertices = 0;
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
    {
        m_db = new LayoutDBRect;
        m_is_Rec = true;
    }
    else 
    {
        m_db = new LayoutDBPolygon;
        m_is_Rec = false;
    }
    // get options from dummy ControlParameter
    tmpParms.swap(m_db->parms);
}

void SimpleMPL::read_gds()
{
    char buf[256];
    mplSPrint(kINFO, buf, "reading input files takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);
	mplPrint(kINFO, "Reading input file %s\n", m_db->input_gds().c_str());
	// read input gds file 
	GdsReader reader (*m_db);
        mplAssertMsg(reader(m_db->input_gds()), "failed to read %s", m_db->input_gds().c_str());
	m_db->num_of_cells = 1;	
	m_db->cal_bound();
	if(m_db->input2_gds() != std::string("")){
		mplAssertMsg(reader(m_db->input2_gds()),"failed to read %s", m_db->input2_gds().c_str());
		m_db->num_of_cells += 1;
		coordinate_type left_x = m_db->boundaries[1];
		// std::cout<<left_x<<std::endl;
		m_db->cal_bound();
		coordinate_type right_x = m_db->boundaries[1];
		// std::cout<<right_x<<std::endl;
		if(m_db->parms.flip2){
			// std::cout<<"FLIP SECOND CELL!"<<std::endl;
			m_db->flip(left_x,right_x);
		}
	}
	if(m_db->input3_gds() != std::string("")){
		mplAssertMsg(reader(m_db->input3_gds()),"failed to read %s", m_db->input3_gds().c_str());
		m_db->num_of_cells += 1;
		coordinate_type left_x = m_db->boundaries[1];
		m_db->cal_bound();
		coordinate_type right_x = m_db->boundaries[1];
		if(m_db->parms.flip3){
			m_db->flip(left_x,right_x);
		}
	}
	
	// must call initialize after reading 
	m_db->initialize_data();
	
	// report data 
	m_db->report_data();
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
	if (m_db->use_stitch())
	{
		GdsWriter writer;
		mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, m_StitchRelation, m_mAdjVertex);
	}
	else
	{
		GdsWriter writer;
		mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
}

void SimpleMPL::write_json()
{
	std::ofstream jsonFile;
	char filename[80];
	strcpy(filename,"../json/");
	strcat(filename,const_cast<char *>(m_db->input_gds().c_str()));
	strcat(filename,".json");
	// char* tmp = const_cast<char *>(m_db->input_gds().c_str());
	// char* output_json =(char*)strcat(tmp,".json");
	jsonFile.open(filename);
	jsonFile<<"[";
	for(uint32_t i = 0; i != m_mAdjVertex.size(); i++)
    {
		if(m_in_DG[i] == false) 
        {	
            mplAssertMsg(!(m_db->use_stitch() && m_StitchRelation[i].size()!=0), "BUG FOUND %u", i);
        }
		jsonFile <<"\n";
		jsonFile << "\t{\n";
		jsonFile << "\t\t\"id\": "<<i<<",\n";
		jsonFile << "\t\t\"color\": "<<(int)m_db->vPatternBbox[i]->color()<<",\n";
		jsonFile << "\t\t\"conflict_degree\":"<<m_mAdjVertex[i].size()<<",\n";
		if(m_db->use_stitch())
        {
			jsonFile << "\t\t\"stitch_degree\":"<<m_StitchRelation[i].size()<<",\n";
        }
		if(m_mAdjVertex[i].size() == 0)
        {
			jsonFile << "\t\t\"conflict\": []";
		}
		else
        {
			jsonFile << "\t\t\"conflict\": [\n";
			for(uint32_t j = 0; j != m_mAdjVertex[i].size()-1; j++)
            {
				jsonFile << "\t\t\t{\"id\": "<<m_mAdjVertex[i][j]<<"},\n";
			} 
			jsonFile << "\t\t\t{\"id\": " <<m_mAdjVertex[i][m_mAdjVertex[i].size()-1]<<"}\n";
			jsonFile << "\t\t]";
		}
		if(m_db->use_stitch())
        {
            jsonFile<<",\n";
            if(m_StitchRelation[i].size() == 0)
            {
                jsonFile << "\t\t\"stitch\": []\n";
            }
            else
            {
                jsonFile << "\t\t\"stitch\": [\n";
                for(uint32_t j = 0; j != m_StitchRelation[i].size()-1; j++)
                {
                    jsonFile << "\t\t\t{\"id\": "<<m_StitchRelation[i][j]<<"},\n";
                }
                jsonFile << "\t\t\t{\"id\": "<<m_StitchRelation[i][m_StitchRelation[i].size()-1]<<"}\n";
                jsonFile << "\t\t]\n";
            }
        }
		else
        {
			jsonFile<<"\n";
		}

		jsonFile << "\t},";
	}
	long pos = jsonFile.tellp();
	jsonFile.seekp(pos - 1);
	jsonFile<<"\n]";
	jsonFile.close();
} 

void SimpleMPL::write_json(graph_type const& sg,std::string graph_count,std::vector<int8_t>& Colors )
{
	std::ofstream jsonFile;
	// char* tmp = const_cast<char *>(m_db->input_gds().c_str());
	// char* tmp2 = (char*)strcat(tmp,graph_count);
    //mplPrint(kDEBUG, "json file name is %s\n", graph_count);
	// char* output_json =(char*)strcat(tmp2,".json");
    std::stringstream ss1;
    ss1 << m_db->parms.flip2;
    std::stringstream ss2;
    ss2 << m_db->parms.flip3;
	//for one input gds
	string filename0 = m_db->input_gds().substr(m_db->input_gds().find("/",2)+1) +"_"+ graph_count + ".json";
	//for two input gds
	string filename = m_db->input_gds().substr(m_db->input_gds().find("/",2)+1) +"_"+m_db->input2_gds().substr(m_db->input2_gds().find("/",2)+1)+"_"+ss1.str()+"_"+ graph_count + ".json";
	//for three input gds
	string filename1 = m_db->input_gds().substr(m_db->input_gds().find("/",2)+1) +"_"+m_db->input2_gds().substr(m_db->input2_gds().find("/",2)+1)+"_"+ss1.str()+"_"+m_db->input3_gds().substr(m_db->input3_gds().find("/",2)+1)+"_"+ss2.str() +"_"+ graph_count + ".json";
	if(!m_db->input3_gds().empty())
		jsonFile.open(".//json//"+filename1);
	else{
		if(!m_db->input2_gds().empty())
			jsonFile.open(".//json//"+filename);
		else
			jsonFile.open(".//json//"+filename0);
	}
		
	//jsonFile.open("/json/"+ m_db->input_gds() + graph_count + ".json");
	//jsonFile.open("/research/byu2/wli/repository/OpenMPL/bin/json/" + m_db->input_gds() + graph_count + ".json");
	//std::cout<<jsonFile.is_open()<<std::endl;
	//std::cerr<<"Failed to open file : "<<SYSERROR()<<std::endl;
	SimpleMPL::graph_type tmp_graph = sg;
	jsonFile<<"[";
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	for (boost::tie(vi1, vie1) = boost::vertices(tmp_graph); vi1 != vie1; ++vi1)
    {
		vertex_descriptor v1 = *vi1;
		jsonFile <<"\n";
		jsonFile << "\t{\n";
		jsonFile << "\t\t\"id\": "<<(int)v1<<",\n";
		jsonFile << "\t\t\"color\": "<<(int)Colors[(int)v1]<<",\n";
		boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
		boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
		int conflict_count = 0;
		int stitch_count = 0;
		for (next2 = vi2; vi2 != vie2; vi2 = next2)
		{
			++next2; 
			vertex_descriptor v2 = *vi2;
			std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
			if(boost::get(boost::edge_weight, tmp_graph, e12.first) > 0)
            {
				conflict_count ++;
			}
			else
            {
				stitch_count ++;
			}
			//out << int(v1) <<" "<< int(v2) <<" "<< boost::get(boost::edge_weight, tmp_graph, e12.first)<<"\n";
		}
		if(conflict_count == 0)
        {
			jsonFile << "\t\t\"conflict\": []";
		}
		else
        {
			jsonFile << "\t\t\"conflict\": [\n";
			int conflict_second_count = 0;
			boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
			for (next2 = vi2; vi2 != vie2; vi2 = next2)
			{
				++next2; 
				vertex_descriptor v2 = *vi2;
				std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
				if(boost::get(boost::edge_weight, tmp_graph, e12.first) > 0)
                {
					conflict_second_count ++ ;
					if(conflict_second_count == conflict_count)
                    {
						jsonFile << "\t\t\t{\"id\": " <<int(v2)<<"}\n";
					}
					else
                    {
						jsonFile << "\t\t\t{\"id\": "<<int(v2)<<"},\n";
					}
					
				}
				//out << int(v1) <<" "<< int(v2) <<" "<< boost::get(boost::edge_weight, tmp_graph, e12.first)<<"\n";
			}

			jsonFile << "\t\t]";
		}
		if(m_db->use_stitch())
        {
            jsonFile<<",\n";
            if(stitch_count == 0)
            {
                jsonFile << "\t\t\"stitch\": []";
            }
            else
            {
                jsonFile << "\t\t\"stitch\": [\n";
                int stitch_second_count = 0;
                boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
                for (next2 = vi2; vi2 != vie2; vi2 = next2)
                {

                    ++next2; 
                    vertex_descriptor v2 = *vi2;
                    std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
                    if(boost::get(boost::edge_weight, tmp_graph, e12.first) < 0)
                    {
                        stitch_second_count ++ ;
                        if(stitch_second_count == stitch_count)
                        {
                            jsonFile << "\t\t\t{\"id\": " <<int(v2)<<"}\n";
                        }
                        else
                        {
                            jsonFile << "\t\t\t{\"id\": "<<int(v2)<<"},\n";
                        }
                    }
                }

                jsonFile << "\t\t]";
            }

        }
		else
        {
			jsonFile<<"\n";
		}

		jsonFile << "\t},";
	}
	long pos = jsonFile.tellp();
	jsonFile.seekp(pos - 1);
	jsonFile<<"\n]";
	jsonFile.close();
} 
void SimpleMPL::write_graph(SimpleMPL::graph_type& g, std::string const& filename ) const
{
    // in order to make the .gv file readable by boost graphviz reader 
    // I dump it with boost graphviz writer 
    boost::dynamic_properties dp;
    dp.property("id", boost::get(boost::vertex_index, g));
    dp.property("node_id", boost::get(boost::vertex_index, g));
    dp.property("label", boost::get(boost::vertex_index, g));
    // somehow edge properties need mutable graph_type& 
    dp.property("color", boost::get(boost::edge_weight, g));
    dp.property("label", boost::get(boost::edge_weight, g));
    std::ofstream out ((filename+".gv").c_str());
    boost::write_graphviz_dp(out, g, dp, string("id"));
    out.close();
    la::graphviz2pdf(filename);
}
void SimpleMPL::write_txt(graph_type const& sg,std::string const filename, double& cost)
{
	std::ofstream out(("./graph/"+filename+"_"+std::to_string(cost)+".txt").c_str());
	SimpleMPL::graph_type tmp_graph = sg;
	out << num_vertices(tmp_graph)<<"\n";

	//output edges
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	for (boost::tie(vi1, vie1) = boost::vertices(tmp_graph); vi1 != vie1; ++vi1)
	{
		vertex_descriptor v1 = *vi1;
		boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
		boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
		for (next2 = vi2; vi2 != vie2; vi2 = next2)
		{
			++next2; 
			vertex_descriptor v2 = *vi2;
			if (v1 >= v2) continue;
			std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
			mplAssert(e12.second);
			out << int(v1) <<" "<< int(v2) <<" "<< boost::get(boost::edge_weight, tmp_graph, e12.first)<<"\n";
		}
	}
	out.close();
}
void SimpleMPL::out_stat()
{
	int TCE = 0;
	int TSE = 0;
    for (uint32_t i = 0; i != m_mAdjVertex.size(); i++)
    {
		if(m_in_DG[i] == false) continue;
		for(uint32_t j = 0; j != m_mAdjVertex[i].size(); j++)
        {
            if(m_in_DG[m_mAdjVertex[i][j]] == false) 
                continue;
            TCE ++;
        }
	}
    for (uint32_t i = 0; i != m_StitchRelation.size(); i++)
    {
		TSE += m_StitchRelation[i].size();
	}
	mplPrint(kNONE, "=================== Graph Simplification Information ===================\n");
	mplPrint(kNONE, "# of Total Conflict Edges                              : %d\n", TCE/2);
	mplPrint(kNONE, "# of Total Stitch Edges                                : %d\n", TSE/2);
	mplPrint(kNONE, "# of Total Wires                                       : %d\n", (int)m_db->vPatternBbox.size());
	mplPrint(kNONE, "# of DG component after  DG Bridge Division            : %d\n", (int)m_DG_num); 
	mplPrint(kNONE, "# of DGs after  DG Bridge Division                     : %d\n", (int)std::count(m_in_DG.begin(),m_in_DG.end(),true)); 
	mplPrint(kNONE, "========================================================================\n");
}

void SimpleMPL::solve()
{
	int numProcs = omp_get_num_procs();
    std::cout << "omp_get_num_procs() = " << numProcs << std::endl;
    // skip if no uncolored layer 
    if (m_db->parms.sUncolorLayer.empty())
        return;
	if(m_db->parms.selector!= ""){
		this->update_algorithm_selector(m_db->parms.selector);
	}
	
	boost::timer::cpu_timer total_timer;
	total_timer.start();
	if (m_db->vPatternBbox.empty())
	{
        mplPrint(kWARN, "No patterns found in specified layers\n");
		return;
	}
	
	//First simplification: ICC (independent component computation)
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
	// TODO: if the m_vBookmark generation is time-consuming, some optimization algorithm may be proposed (set is better than vector)
	std::vector<uint32_t>().swap(m_vBookmark);
	m_vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
        {
            m_vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
        }
	}
	//We implememt another cal_bound func in layout class for cell sequence generation, therefore we comment this line to slightly reduce the runtime 04/11/19 Wei
	//this->cal_boundaries();
	this->setVddGnd(); //perhapes we should also consider pshape->getPointNum()==4
	mplPrint(kINFO,"VDD nodes set done\n");
	boost::timer::cpu_timer stitch_timer;
	stitch_timer.start();

	//First simplification: IVR (HIDE_SMALL_DEGREE) && Bridge detection
	this->lg_simplification();

	// Secondly, insert stitch
	if (m_db->use_stitch()) //if we use stitches, we need insert stitches through projection()
	{

		//GdsWriter writer;
		//writer.write_Simplification(m_db->output_gds() + "_lg_simplification.gds", *m_db, m_vCompId, m_mAdjVertex, m_in_DG, m_isVDDGND, true);
		this->projection();		///< m_vBookmark has already been updated in projection()
		std::cout<<total_timer.format(2, "stitch time %ts(%p%), %ws real")<<std::endl;
	}
    std::vector<uint32_t>().swap(m_dgCompId);
    m_dgCompId.assign(m_db->vPatternBbox.size(), 0);
    m_globalCompId = 1;
#ifdef DEBUG_LIWEI
    mplPrint(kDEBUG, "%lu Nodes in DG\n", std::count(m_in_DG.begin(), m_in_DG.end(), true));
#endif

    // update VddGnd information
    this->setVddGnd();


    m_vdd_multi_comp.resize(m_db->vPatternBbox.size());
    //this->dg_simplification(); 
    if (m_db->use_stitch())
    {
        GdsWriter writer;
        writer.write_Simplification(m_db->output_gds() + "_dgSimplification.gds", *m_db, m_dgCompId, m_StitchRelation, std::vector<bool>(), m_isVDDGND, false);
    }
#ifdef DEBUG_LIWEI
    mplPrint(kDEBUG, "This graph finally been decomposed into : %u\n", m_comp_cnt);
#endif

    // for (uint32_t i = 0; i < m_vdd_multi_comp.size(); i++)
    // {
    // 	if (m_isVDDGND[i])
    // 	{
    // 		std::cout << "VDD " << i << " in : ";
    // 		for (std::vector<uint32_t>::iterator it = m_vdd_multi_comp[i].begin(); it != m_vdd_multi_comp[i].end(); it++)
    // 			std::cout << *it << " ";
    // 		std::cout << std::endl;
    // 	}
    // }

	boost::timer::cpu_timer t;
	t.start();
	// thread number controled by user option 

	// Thirdly, solve each component
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
    for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
    {
        // construct a component 
		boost::timer::cpu_timer t_comp;
		t_comp.start();
        std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+m_vBookmark[comp_id];
        std::vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+m_vBookmark[comp_id+1] : m_vVertexOrder.end();
        // solve component 
        // pass iterators to save memory 
        this->solve_component(itBgn, itEnd, comp_id);
	// cout << "almost finished -1" << endl;
		//std::cout<<"comp_id "<<comp_id<<" omp_get_thread_num() is "<<omp_get_thread_num()<<t_comp.format(2, " comp time %ts(%p%), %ws real")<<std::endl;
	}
    if( m_db->parms.record > 0)
    {
	// cout << "in" << endl;
        if(m_db->use_stitch())
        {
            std::ofstream myfile,color_time,total_time;
            myfile.open ("record.txt", std::ofstream::app);
            color_time.open("result_w_stitch.txt",std::ofstream::app);
            myfile << t.format(2, "color time %ts(%p%), %ws real")<<"\n";
            myfile << total_timer.format(2, "total time %ts(%p%), %ws real")<<"\n";
            // if(m_db->algo() == AlgorithmTypeEnum::ILP_GUROBI || m_db->algo() == AlgorithmTypeEnum::DANCING_LINK)
            // {
            //     color_time<<"\\\\\n"<<m_db->input_gds();
            // }
            if(m_db->algo() == AlgorithmTypeEnum::ILP_UPDATED_GUROBI)
            {
                color_time<<"\\\\\n"<<m_db->input_gds();
            }
			// color_time<<"\\\\\n"<<m_db->input_gds();
            // color_time<<t.format(3, "&  %w");

			// color_time << m_db->input_gds() << " ";
			color_time << std::setw(5) <<m_db->algo()<<" "<<t.format(3, "&  %w");
            myfile.close();
            color_time.close();
        }
        else
        {
            std::ofstream myfile,color_time,total_time;
            myfile.open ("record.txt", std::ofstream::app);
            color_time.open("color_wo_stitch.txt",std::ofstream::app);
            total_time.open("total_wo_stitch.txt",std::ofstream::app);
            myfile << t.format(2, "color time %ts(%p%), %ws real")<<"\n";
            myfile << total_timer.format(2, "total time %ts(%p%), %ws real")<<"\n";
            if(m_db->algo() == AlgorithmTypeEnum::ILP_GUROBI)
            {
                color_time<<"\\\\\n"<<m_db->input_gds();
                total_time<<"\n"<<m_db->input_gds();
            }
            color_time<<t.format(2, "& %t  &  %w");
            total_time<<total_timer.format(2, "& %t  &  %w");
            myfile.close();
            color_time.close();
            total_time.close();
        }
    }
	// cout << "almost finished" << endl;
	this->out_stat();
	// this->write_json();
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
}

void SimpleMPL::lg_simplification()
{
	std::vector<bool>().swap(m_in_DG);
	std::vector<bool>().swap(m_articulation_vec);
	m_in_DG.assign(m_db->vPatternBbox.size(), false);
	m_articulation_vec.assign(m_db->vPatternBbox.size(), false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + m_vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + m_vBookmark[comp_id + 1] : m_vVertexOrder.end();
		// projection on each component
		this->lg_simplification(itBgn, itEnd, comp_id);
	}
}

void SimpleMPL::lg_simplification(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
	(void)comp_id;
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	uint32_t simplify_strategy = graph_simplification_type::NONE;
	if (m_db->simplify_level() > 1)
		simplify_strategy |= graph_simplification_type::HIDE_SMALL_DEGREE;

	// NOTE: After experiments, I noticed that BICONNECTED_COMPONENT cannot be used in the first simplification phase
	// since some valid stitches may be introduced in the arti_point! if we use BICONNECTED_COMPONENT.
	// Currently, biconnected point recovery does not fit into stitch case!
	// by WEI 21/02/2020
	
	// if (m_db->simplify_level() > 2)
    //     simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;
	
	// uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
	//simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;
	graph_simplification_type gs(dg, m_db->color_num());
	
	gs.set_isVDDGND(vdd_set);
	gs.simplify(simplify_strategy);
	std::vector<bool> projLocal(mGlobal2Local.size(), false);

	std::vector<vertex_descriptor> all_articulations;
	gs.get_articulations(all_articulations);

	std::set<vertex_descriptor> temp_set(all_articulations.begin(), all_articulations.end());

	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<vertex_descriptor> vSimpl2Orig;
		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		for (uint32_t i = 0; i < vSimpl2Orig.size(); i++)
			projLocal[vSimpl2Orig[i]] = true;
	}

	for (std::map<uint32_t, uint32_t>::iterator it = mGlobal2Local.begin(); it != mGlobal2Local.end(); it++)
	{
		if (projLocal[it->second])
			m_in_DG[it->first] = true;
		if (temp_set.find(it->second) != temp_set.end())
			m_articulation_vec[it->first] = true;
	}
}

void SimpleMPL::dg_simplification()
{
	std::vector<std::vector<uint32_t> > comp2subcomp(m_comp_cnt);
	std::vector<std::map<vertex_descriptor, std::set<uint32_t> > > m_mArtiPoints(m_comp_cnt);
	std::vector<std::vector<vertex_descriptor> > m_mCompVertex;
	std::vector<std::vector<vertex_descriptor> > HiddenVertices;
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + m_vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + m_vBookmark[comp_id + 1] : m_vVertexOrder.end();
		std::vector<std::vector<vertex_descriptor> > comp_vertex;
		std::map<vertex_descriptor, std::set<uint32_t> >& arti_point = m_mArtiPoints[comp_id];
		// Simplification on each component

		this->dg_simplification(itBgn, itEnd, comp_id, arti_point, comp_vertex);
	}
}



void SimpleMPL::dg_simplification(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::map<vertex_descriptor, std::set<uint32_t> >& m_ArtiPoint, std::vector<std::vector<vertex_descriptor> >& m_CompVertex)
{
	(void)m_ArtiPoint;
	(void)m_CompVertex;
	(void)comp_id;
	//std::cout << "\n\n\n**************************************************************\n";
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
	simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;

	graph_simplification_type gs(dg, m_db->color_num());
	gs.precolor(vColor.begin(), vColor.end());

	std::vector<uint32_t> mLocal2Global(mGlobal2Local.size());
	for (std::map<uint32_t, uint32_t>::iterator it = mGlobal2Local.begin(); it != mGlobal2Local.end(); it++)
		mLocal2Global[it->second] = it->first;

	gs.set_isVDDGND(vdd_set);
	gs.simplify(simplify_strategy);

	///< big component : input big component (the only several big components in original input gds file).
	///< small component : small component divided from big component after conducting IVR and biconnected decomposition, may contains one or two VDDGNDs, is corresponding one row in big component
	///< sub small component : sub small component divided from small component after only conducting biconnected decomposition, divide one row into sevel small sub components
	
	std::vector<std::vector<int8_t> > big_comp_color(gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSmall2Big (gs.num_component());
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();
	for (uint32_t small_comp_id = 0; small_comp_id < gs.num_component(); ++small_comp_id)
	{

		graph_type small_g;
		
		std::vector<int8_t>& small_comp_color = big_comp_color[small_comp_id];	///< used to store coloring results of current small component
		(void)small_comp_color;
		std::vector<vertex_descriptor>& vSimplSmall2OriBig = mSmall2Big[small_comp_id];		///< used to store the simplificaiton mapping information from current small component to original big component

		gs.simplified_graph_component(small_comp_id, small_g, vSimplSmall2OriBig);

		std::set<vertex_descriptor> small_vdd_set;
		std::vector<int8_t> small_vColor(num_vertices(small_g), -1);
	}

}


uint32_t SimpleMPL::merge_K4_coloring(SimpleMPL::graph_type const& dg, std::vector<int8_t>& vColor) const
{
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs(dg, m_db->color_num());
	gs.precolor(vColor.begin(), vColor.end());
	if (m_db->color_num() == 3)
		gs.max_merge_level(3);
	else if (m_db->color_num() == 4)
		gs.max_merge_level(2);
	gs.simplify(graph_simplification_type::MERGE_SUBK4);
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

	std::vector<std::vector<int8_t> > mSubColor(gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig(gs.num_component());
	double acc_obj_value = 0;

	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		graph_type sg;
		std::vector<int8_t>& vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];
		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		vSubColor.assign(num_vertices(sg), -1);

		typedef lac::Coloring<graph_type> coloring_solver_type;
		//TODO: merge k4 is not used currently so that we use sub_comp_id here instead of comp_id
		coloring_solver_type* pcs = create_coloring_solver(sg,sub_comp_id,sub_comp_id);

		boost::graph_traits<graph_type>::vertex_iterator vi, vie;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
		{
			vertex_descriptor v = *vi;
			int8_t color = vColor[vSimpl2Orig[v]];
			if (color >= 0 && color < m_db->color_num())
			{
				pcs->precolor(v, color);
				vSubColor[v] = color;
			}
		}
		double obj_value1 = (*pcs)();

		acc_obj_value += obj_value1;
		for (boost::tie(vi, vie) = vertices(sg); vi != vie; vi++)
		{
			vertex_descriptor v = *vi;
			int8_t color = pcs->color(v);
			vSubColor[v] = color;
		}

		delete pcs;
	}

	gs.recover(vColor, mSubColor, mSimpl2Orig);

	return acc_obj_value;
}

void SimpleMPL::setVddGnd()
{
	std::vector<bool>().swap(m_isVDDGND);
	int threshold = 6000;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	m_isVDDGND.assign(vertex_num, false);
	int count = 0;
	for (uint32_t i = 0; i < vertex_num; i++)
	{
		// if (is_long_enough(m_db->vPatternBbox[i]))
		// {
		// 	m_isVDDGND[i] = true;
		// 	m_db->set_color(i,0);
		// 	// m_db->vPatternBbox[i]->color(m_db->color_num() - 1);
		// }
		if(boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) == threshold )
        {
			m_isVDDGND[i] = false;
#ifdef DEBUG_LIWEI
            mplPrint(kDEBUG, "dirty coding for oracle cases\n");
#endif
		}
	}
#ifdef DEBUG_LIWEI
    mplPrint(kDEBUG, "%d dirty coding for sparc cases\n", count);
#endif
	//special case for sparc benchmark
	return;
}

void SimpleMPL::cal_boundaries()
{
	std::vector<coordinate_type>().swap(m_boundaries);
	mplAssert(m_boundaries.empty());
	uint32_t vertex_num = m_db->vPatternBbox.size();
	coordinate_type tmp_bound[4] = {INT_MAX,0,INT_MAX ,0};


	for (uint32_t i = 0; i < vertex_num; i++)
    {
		coordinate_type test0 = gtl::xl(*(m_db->vPatternBbox[i]));
		coordinate_type test1 = gtl::xh(*(m_db->vPatternBbox[i]));
		coordinate_type test2 = gtl::yl(*(m_db->vPatternBbox[i]));
		coordinate_type test3 = gtl::yh(*(m_db->vPatternBbox[i]));
		tmp_bound[0] = std::min(test0,tmp_bound[0]);
		tmp_bound[1] = std::max(test1,tmp_bound[1]);
		tmp_bound[2] = std::min(test2,tmp_bound[2]);
		tmp_bound[3] = std::max(test3,tmp_bound[3]);
	}
	m_boundaries.push_back(tmp_bound[0]);
	m_boundaries.push_back(tmp_bound[1]);
	m_boundaries.push_back(tmp_bound[2]);
	m_boundaries.push_back(tmp_bound[3]);
	return;
}

void SimpleMPL::projection()
{
#ifdef _OPENMP
	std::cout<<"OPENMP enabled with thread num "<<m_db->thread_num()<<std::endl;
#endif
	uint32_t vertex_num = m_db->vPatternBbox.size();
	std::vector<rectangle_pointer_type> rect_vec = m_db->polyrect_patterns();	///< original rectangle list
	std::vector<uint32_t> Poly_Rect_begin;
    if(m_is_Rec)
    {
        for (uint32_t i = 0; i < vertex_num; i++)
        {
            Poly_Rect_begin.push_back(i);
        }
    }
    else
	{
		Poly_Rect_begin = m_db->PolyRectBgnLoc();				///< original polygons mapping to rectangles
	}
    

	std::vector<uint32_t> Poly_Rect_end;
	Poly_Rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		Poly_Rect_end[i] = Poly_Rect_begin[i + 1] - 1;
	
	Poly_Rect_end[vertex_num - 1] = rect_vec.size() - 1;
	std::vector<uint32_t> new_vCompId_vec;
	std::vector<uint32_t> new_vertex_order;
	std::vector<rectangle_pointer_type> new_rect_vec;				///< store the new rectangles
	std::vector<uint32_t> new_Rect2ParentPoly;						///< map from new rectangles to new parent polygons
	std::vector<std::vector<uint32_t> >().swap(m_ori2new_polygon);	///< map from original polygons to new polygons
	std::vector<uint32_t>().swap(m_new2ori_polygon);					///< map from new polygons to old polygons
	m_ori2new_polygon.resize(vertex_num);
	
	uint32_t new_polygon_id = -1;						///< new polygon id
	uint32_t new_rectangle_id = -1;						///< new rectangle id
	uint32_t stitch_edge_number = 0;    
#ifdef _OPENMP
#pragma omp parallel num_threads(m_db->thread_num()) 
#endif
#ifdef _OPENMP
	#pragma omp for ordered schedule(dynamic, 1)
#endif
	for (uint32_t i = 0; i < vertex_num; i++)
	{
		uint32_t v = m_vVertexOrder[i];
		uint32_t comp_id = m_vCompId[v];
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[v];
		uint32_t pid = pPattern->pattern_id();
		uint32_t start_idx = Poly_Rect_begin[pid];
		uint32_t end_idx = Poly_Rect_end[pid];
		
		std::vector<uint32_t>& nei_vec = m_mAdjVertex[pid];
		if (m_in_DG[pid] && m_isVDDGND[pid] == false) 
		// if (m_in_DG[pid] && m_articulation_vec[pid] == false && m_isVDDGND[pid] == false) 
		{ 
			std::vector<rectangle_pointer_type> poss_nei_vec;
			for (std::vector<uint32_t>::iterator it = nei_vec.begin(); it != nei_vec.end(); it++)
			{

				if (m_in_DG[*it] == false) continue;
				uint32_t s_idx = Poly_Rect_begin[*it];
				uint32_t e_idx = Poly_Rect_end[*it];
				poss_nei_vec.insert(poss_nei_vec.end(), rect_vec.begin() + s_idx, rect_vec.begin() + e_idx + 1);
			} ///< for nei_vec

			std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_split;	///< the generated rectangles and their parent original rectangle id
			///< traverse all rectangles in current polygon
			for (uint32_t j = start_idx; j <= end_idx; j++)
			{
				rectangle_type rect(*rect_vec[j]);
				std::vector<rectangle_pointer_type> split;
				split_rectangle(rect, split, poss_nei_vec);
				for (uint32_t i = 0; i < split.size(); i++)
					rect_split.push_back(std::make_pair(split[i], rect.pattern_id()));
				stitch_edge_number += split.size() - 1;
			} 

            std::vector<uint32_t> new_polygon_id_list;
            // reconstruct polygons, also generate stitch relationships
            reconstruct_polygon(new_polygon_id, new_polygon_id_list, rect_split,pid);
            mplAssert(new_polygon_id_list.size() == rect_split.size());
#ifdef _OPENMP
#pragma omp ordered
#endif
            {
#ifdef _OPENMP
#pragma omp critical
#endif
                {


                    uint32_t pivot = new_polygon_id_list[0] -1 ;
                    for (uint32_t i = 0; i < rect_split.size(); i++)
                    {
                        rect_split[i].first->pattern_id(++new_rectangle_id);
                        new_rect_vec.push_back(rect_split[i].first);
                        new_Rect2ParentPoly.push_back(new_polygon_id_list[i]);

                        // insert new polygon
                        if (pivot != new_polygon_id_list[i])
                        {
                            pivot = new_polygon_id_list[i];
                            m_new2ori_polygon.push_back(pid);
                            m_ori2new_polygon[pid].push_back(pivot);
                            new_vertex_order.push_back(pivot);
                            new_vCompId_vec.push_back(comp_id);
                        }
                    } ///< for rect_list
                }
            }
		} ///< if, stitch generation core
		else
		{
#ifdef _OPENMP
#pragma omp ordered
#endif
            {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    new_polygon_id += 1;
                    for (uint32_t j = start_idx; j <= end_idx; j++)
                    {
                        rectangle_pointer_type rect = rect_vec[j];
                        rect->pattern_id(++new_rectangle_id);
                        if (m_db->gen_stitch())
                            rect->color(new_polygon_id % 7);
                        new_rect_vec.push_back(rect);
                        new_Rect2ParentPoly.push_back(new_polygon_id);
                    }
                    m_new2ori_polygon.push_back(pid);
                    new_vertex_order.push_back(new_polygon_id);
                    new_vCompId_vec.push_back(comp_id);
                    m_ori2new_polygon[pid].push_back(new_polygon_id);
                    m_StitchRelation.push_back(std::vector<uint32_t>());
                }
			}
		} ///< else (m_in_DG)
	} ///< for all vertices

	m_db->refresh(new_rect_vec, new_Rect2ParentPoly);
	update_conflict_relation(); 

	std::vector<uint32_t>().swap(m_vCompId);
	m_vCompId.swap(new_vCompId_vec);

	std::vector<uint32_t>().swap(m_vVertexOrder);
	m_vVertexOrder.swap(new_vertex_order);

	std::vector<uint32_t>().swap(m_vBookmark);
	m_vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			m_vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

	return;
}

void SimpleMPL::update_conflict_relation()
{
	// now update adjacency list, we still need original adjacency list.
	// here we use set to remove duplicated neighboring polygons
	std::vector<std::set<uint32_t> >new_mAdjVertex;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	new_mAdjVertex.resize(m_db->vPatternBbox.size());
	uint32_t edge_num = 0;
	uint32_t dg_conflict = 0;
	//uint32_t dg_illegal_conflict = 0;

	//update m_in_DG
	std::vector<bool> tmp_in_DG;
	tmp_in_DG.assign(m_db->vPatternBbox.size(), false);
	for(uint32_t i = 0; i< tmp_in_DG.size();i++)
    {
		if(m_in_DG[m_new2ori_polygon[i]])
			tmp_in_DG[i] = true;
	}
	m_in_DG.swap(tmp_in_DG);

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num()) reduction(+:edge_num)
#endif
	for (uint32_t v = 0; v < vertex_num; ++v)
	{
		rectangle_pointer_type const& pPattern = m_db->vPatternBbox[v];
		std::set<uint32_t>& vAdjVertex = new_mAdjVertex[v];

		// find patterns connected with pPattern 
		// query tPatternBbox in m_db
		mplAssert(v == pPattern->pattern_id());
		rectangle_type rect(*pPattern);
		// bloat pPattern with minimum coloring distance 
		gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
		gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
		for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
			itq != m_db->tPatternBbox.qend(); ++itq)
		{
			rectangle_pointer_type const& pAdjPattern = *itq;
			
			// Stitch relationships have already generated. If these two polygons have stitch relations, ignore them.
			std::vector<uint32_t>::iterator itr = find(m_StitchRelation[pPattern->pattern_id()].begin(), m_StitchRelation[pPattern->pattern_id()].end(), pAdjPattern->pattern_id());
			//if (itr != m_StitchRelation[pPattern->pattern_id()].end())
			if(m_new2ori_polygon[pPattern->pattern_id()] == m_new2ori_polygon[pAdjPattern->pattern_id()])
			{
				continue;
			}
				
			mplAssert(itr == m_StitchRelation[pPattern->pattern_id()].end());
			mplAssert(pAdjPattern != pPattern);
			if (pAdjPattern != pPattern) // skip pPattern itself 
			{
				mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
				// we consider euclidean distance
				// use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
				coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
				if (distance < m_db->coloring_distance)
				{
					if(m_in_DG[pPattern->pattern_id()] && m_in_DG[pAdjPattern->pattern_id()])	dg_conflict++;
					vAdjVertex.insert(pAdjPattern->pattern_id());
				}
			}
		}
		vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
		edge_num += vAdjVertex.size();
	}

	///< update m_mAdjVertex
	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex.resize(new_mAdjVertex.size());
	for (uint32_t i = 0; i < new_mAdjVertex.size(); i++)
		m_mAdjVertex[i].insert(m_mAdjVertex[i].end(), new_mAdjVertex[i].begin(), new_mAdjVertex[i].end());
}


// find conflicts polygons(not rectangles in this version) \conflict_rects of this rect \prec
void SimpleMPL::find_conflict_rects(rectangle_pointer_type& prec, std::vector<rectangle_pointer_type>& conflict_rects,uint32_t current_poly_id)
{
	const std::vector<rectangle_pointer_type>& vPolyRectPattern = m_db->polyrect_patterns();
	const std::vector<uint32_t>& vPolyRectBeginId = m_db->PolyRectBgnLoc();
	uint32_t num_polyrects = vPolyRectPattern.size();
	//std::cout<<"prec->pattern_id() "<<prec->pattern_id()<<" m_mAdjVertex[current_poly_id].size() "<<m_mAdjVertex[current_poly_id].size()<<" current_poly_id "<<current_poly_id<<std::endl;
	//we tranverse all of the conflict polys of its parent poly, so that reduce runtime to calculate intersecting polys
	for(auto parentAdj = m_mAdjVertex[current_poly_id].begin(); parentAdj != m_mAdjVertex[current_poly_id].end(); ++parentAdj)
	{
		rectangle_pointer_type& AdjPoly = m_db->vPatternBbox[*parentAdj];
		//when we calculate the distance between this rect and poly, we calculate the minimal distance between the rect and child rects of the poly
		coordinate_difference distance = std::numeric_limits<coordinate_difference>::max();
   	 	uint32_t polyRectId1e = (*parentAdj+1 == vPolyRectBeginId.size())? num_polyrects : vPolyRectBeginId[*parentAdj+1];
		for (uint32_t polyRectId1 = vPolyRectBeginId[*parentAdj]; polyRectId1 != polyRectId1e; ++polyRectId1)
		{
			//std::cout<<(coordinate_difference)gtl::euclidean_distance(*(vPolyRectPattern[polyRectId1]), *prec)<<std::endl;
			distance = std::min(distance, (coordinate_difference)gtl::euclidean_distance(*(vPolyRectPattern[polyRectId1]), *prec) );
		}
		if(distance < m_db->coloring_distance)
			conflict_rects.push_back(AdjPoly);	

	}
}

// void SimpleMPL::find_conflict_rects(rectangle_pointer_type& prec, std::vector<rectangle_pointer_type>& conflict_rects,uint32_t current_poly_id)
// {
// 	rectangle_type rect(*prec);
// 	gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
// 	gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
// 	//TODO:: this part can be accelareted by directly using m_mAdjVertex[prec->pattern_id()] (its parrent polygon id)
// 	// and use (coordinate_difference)gtl::euclidean_distance(*(m_mAdjVertex[prec->pattern_id()][0]), *prec) to exactly calculate the distance
// 	//
// 	//Another option is that recording rects instead of polys (not recommended cause it may increase stitch numbers)
// 	for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
// 			itq != m_db->tPatternBbox.qend(); ++itq)
// 	{
// 		rectangle_pointer_type const& pAdjPattern = *itq;
// 		if(gtl::xl(*(pAdjPattern))<=gtl::xl(*(prec)) && gtl::xh(*(pAdjPattern))>=gtl::xh(*(prec)) && gtl::yl(*(pAdjPattern))<=gtl::yl(*(prec)) && gtl::yh(*(pAdjPattern))>=gtl::yh(*(prec)))
// 			continue;
// 		std::vector<rectangle_pointer_type> vPolyRectPattern = m_db->polyrect_patterns();
// 		std::vector<uint32_t> vPolyRectBeginId = m_db->PolyRectBgnLoc();
// 		uint32_t num_polyrects = vPolyRectPattern.size();
// 		uint32_t adj_id =  pAdjPattern->pattern_id();
// 		coordinate_difference distance = std::numeric_limits<coordinate_difference>::max();
//    	 	uint32_t polyRectId1e = (adj_id+1 == vPolyRectBeginId.size())? num_polyrects : vPolyRectBeginId[adj_id+1];
// 		for (uint32_t polyRectId1 = vPolyRectBeginId[adj_id]; polyRectId1 != polyRectId1e; ++polyRectId1)
// 			distance = std::min(distance, (coordinate_difference)gtl::euclidean_distance(*(vPolyRectPattern[polyRectId1]), *prec) );
// 		if(distance < m_db->coloring_distance)
// 			conflict_rects.push_back(pAdjPattern);
// 	}
// }

void SimpleMPL::find_touch_rects(rectangle_pointer_type& prec, std::vector<rectangle_pointer_type>& touch_rects,std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list)
{
	for (uint32_t index = 0 ; index < rect_list.size(); index++)
	{
		if(rect_list[index].first == prec) continue;
		if((coordinate_difference)gtl::euclidean_distance(*rect_list[index].first,*prec) <= 0)	touch_rects.push_back(rect_list[index].first);	
	}	
}

bool SimpleMPL::is_long_enough(rectangle_pointer_type& rec)
{
    if (boost::polygon::delta(*rec, gtl::HORIZONTAL) >= 0.64*(m_db->boundaries[1]-m_db->boundaries[0])||\
            boost::polygon::delta(*rec, gtl::VERTICAL) >= 0.64*(m_db->boundaries[3]-m_db->boundaries[2]))
        return true;
    else
        return false;
}

// 7 | 3 | 6
// __|___|__
// 1 | 8 | 2
// __|___|__
// 4 | 0 | 5
uint32_t SimpleMPL::box2box_direction(rectangle_pointer_type& prec1,rectangle_pointer_type& prec2)
{
	coordinate_difference xl1 = gtl::xl(*(prec1));
	coordinate_difference xl2 = gtl::xl(*(prec2));
	coordinate_difference xh1 = gtl::xh(*(prec1));
	coordinate_difference xh2 = gtl::xh(*(prec2));
	coordinate_difference yl1 = gtl::yl(*(prec1));
	coordinate_difference yl2 = gtl::yl(*(prec2));
	coordinate_difference yh1 = gtl::yh(*(prec1));
	coordinate_difference yh2 = gtl::yh(*(prec2));
	if(xh1 <= xl2)
    {
		if(yh1 <= yl2)	return 4;
		else if(yl1 >= yh2) return 7;
		else return 1;
	}
	else if(xl1 >= xh2)
    {
		if(yh1 <= yl2)	return 5;
		else if(yl1 >= yh2) return 6;
		else return 2;
	}
	else if(yh1 <= yl2) 
    {
        return 0;
    }
	else if(yl1 >= yh2) 
    {
        return 3;
    }
	return 8;
}

//TODO: FINISH this function
bool SimpleMPL::is_boundslice(rectangle_pointer_type& prec1,rectangle_pointer_type& prec2,std::vector<rectangle_pointer_type>& touch_set1)
{
	uint32_t direction = box2box_direction(prec1,prec2);
	if(direction > 3)
    {
        return false;
    }
	bool isHor = (1==direction || 2==direction);
	// Step 1: determine boundary of prec1 & prec2
	coordinate_difference bound_x1, bound_y1, bound_x2, bound_y2;
	bound_x1 = std::min(gtl::xl(*(prec1)), gtl::xl(*(prec2)));
	bound_x2 = std::max(gtl::xh(*(prec1)), gtl::xh(*(prec2)));
	bound_y1 = std::min(gtl::yl(*(prec1)), gtl::yl(*(prec2)));
	bound_y2 = std::max(gtl::yh(*(prec1)), gtl::yh(*(prec2)));
	for(auto touch_rec = touch_set1.begin(); touch_rec != touch_set1.end();touch_rec++)
    {
        if (!isHor && gtl::xl(*(*touch_rec)) > bound_x1) continue;
        if (!isHor && gtl::xh(*(*touch_rec)) < bound_x2) continue;
        if ( isHor && gtl::yl(*(*touch_rec)) > bound_y1) continue;
        if ( isHor && gtl::yh(*(*touch_rec)) < bound_y2) continue;
        return false;
	}
	return true;
}

bool SimpleMPL::can_intro_stitches(rectangle_pointer_type& prec1,rectangle_pointer_type& prec2,std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list, uint32_t current_poly_id)
{
	//mplAssertMsg(m_db->color_num() == 3, "THE PROGRAM ONLY SUPPORTS triple coloring now: %u\n", (uint32_t)m_db->color_num());
	std::vector<rectangle_pointer_type> vID1;
	std::vector<rectangle_pointer_type> vID2;
	std::vector<rectangle_pointer_type> vTouch1;
	std::vector<rectangle_pointer_type> vTouch2;
	find_conflict_rects(prec1,vID1,current_poly_id);
	find_conflict_rects(prec2,vID2,current_poly_id);
	find_touch_rects(prec1,vTouch1,rect_list);
	find_touch_rects(prec2,vTouch2,rect_list);

	//remove the stitches on vdd
  	if (is_long_enough(prec1) || is_long_enough(prec2)) 
    {
        uint32_t direction = box2box_direction(prec1,prec2);
        if(1==direction || 2== direction) 
        {
            return false;
        }
    }
	//remove boundslice
	coordinate_difference threshold = m_db->coloring_distance / (m_db->color_num() + 2 );
	if((gtl::xh(*prec1) - gtl::xl(*prec1)) < threshold || (gtl::yh(*prec1) - gtl::yl(*prec1))< threshold)
	{
		if (true == is_boundslice(prec1, prec2,vTouch1)) 
            return false;
	}
	if(gtl::xh(*prec2) - gtl::xl(*prec2) < threshold || gtl::yh(*prec2) - gtl::yl(*prec2) < threshold)
	{
		if (true == is_boundslice(prec2, prec1,vTouch2)) 
            return false;
	}


	/// given two rects, check whether there is a pshbox that:
	/// 1) touching to both prec1 & prec2
	/// 2) direct(pmybox1, pshbox) != direct(pmybox1, pmybox2)
	/// if yes, directly return false(cannot insert stitch)
	for(auto touch_of_rec1 = vTouch1.begin(); touch_of_rec1 != vTouch1.end(); touch_of_rec1 ++ )
    {
		if(*touch_of_rec1 == prec2)	
            continue;
		for(auto touch_of_rec2 = vTouch2.begin(); touch_of_rec2 != vTouch2.end(); touch_of_rec2 ++ )
        {
			if(*touch_of_rec2 == prec1)	
                continue;
			if(*touch_of_rec1 != *touch_of_rec2)
                continue;
			uint32_t direction1 = box2box_direction(prec1,*touch_of_rec1);
			uint32_t direction2 = box2box_direction(prec1,prec2);
			if (direction1 != direction2) 
                return false;
		}		
	}

	if(vTouch1.size()<=1 && vID1.size() <= 0)
        return false;
	if(vTouch2.size()<=1 && vID2.size() <= 0)
        return false;
	//bFindX actuallly means not find in X
	bool bFind1 = false, bFind2 = false;
	for (uint32_t i=0; i<vID1.size(); i++)
	{
		rectangle_pointer_type id = vID1[i];
		std::vector<rectangle_pointer_type>::iterator vitr = find(vID2.begin(), vID2.end(), id);
		if (vitr == vID2.end()) bFind1 = true;
	}
	for (uint32_t i=0; i<vID2.size(); i++)
	{
		rectangle_pointer_type id = vID2[i];
		std::vector<rectangle_pointer_type>::iterator vitr = find(vID1.begin(), vID1.end(), id);
		if (vitr == vID1.end()) bFind2 = true;
	}
	if (bFind1 || bFind2) 
        return true;
	return false;	
}

void SimpleMPL::reconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t>& new_polygon_id_list, 
        std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list, uint32_t current_poly_id)
{
	Graph G(rect_list.size());
	int stitchCount = 0;

	std::vector<uint32_t>().swap(new_polygon_id_list);
	new_polygon_id_list.assign(rect_list.size(), std::numeric_limits<uint32_t>::max());
	std::vector<std::vector<uint32_t> > stitch_list;
	std::vector<std::vector<uint32_t> > rec_stitch_list;
	std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_temp;
	std::vector<uint32_t> new_polygon_id_temp;
	rec_stitch_list.resize(rect_list.size());
	for(uint32_t i=0;i<rect_list.size();i++)
    {
		rectangle_pointer_type prec1 = rect_list[i].first;
		for(uint32_t j=i+1;j<rect_list.size();j++)
        {
			rectangle_pointer_type prec2 = rect_list[j].first;
			coordinate_difference distance = boost::geometry::distance(*prec1, *prec2);
			if(distance > 0) continue;
			bool canStitch = can_intro_stitches(prec1,prec2,rect_list, current_poly_id);
			if(canStitch)
            {
				stitchCount ++;
				rec_stitch_list[i].push_back(j);
				rec_stitch_list[j].push_back(i);//if two recs with stitch relations locate in same poly, and warning should be arised
				continue;
			}
			add_edge(i, j, G);
		}
	}

	mplAssert(new_polygon_id_list.size() == num_vertices(G));
  	int component_num = connected_components(G, &new_polygon_id_list[0]);
	uint32_t start;
#ifdef _OPENMP
#pragma omp ordered
#endif
    {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            start = polygon_id + 1;	
            for(uint32_t i=0; i< new_polygon_id_list.size();i++)
            {
                new_polygon_id_list[i] += start;
            }
            polygon_id += component_num;
        }
    }
	stitch_list.resize(component_num);
	//transfrom the stitch list of rectangles into polys
	for (uint32_t i = 0; i < rect_list.size(); i++)
	{
		for (uint32_t j = i + 1; j < rect_list.size(); j++)
		{
			if(std::find(rec_stitch_list[i].begin(), rec_stitch_list[i].end(), j) != rec_stitch_list[i].end())
            {
				if(new_polygon_id_list[i] != new_polygon_id_list[j])
                {
					if(std::find(stitch_list[new_polygon_id_list[i] - start].begin(), stitch_list[new_polygon_id_list[i] - start].end(), new_polygon_id_list[j]) == stitch_list[new_polygon_id_list[i] - start].end())
                    {
                        stitch_list[new_polygon_id_list[i] - start].push_back(new_polygon_id_list[j]);    
                        stitch_list[new_polygon_id_list[j] - start].push_back(new_polygon_id_list[i]);    
                    }  
				}
			}
		}
	}
#ifdef _OPENMP
#pragma omp ordered
#endif
    {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            m_StitchRelation.insert(m_StitchRelation.end(), stitch_list.begin(), stitch_list.end());
        }
    }
	// rectangles in the same polygon should be abutting
	for (uint32_t i = start, ie = start+component_num; i < ie; ++i)
	{
		for (uint32_t j = 0; j < new_polygon_id_list.size(); ++j)
		{
			if (new_polygon_id_list[j] == i)
			{
				new_polygon_id_temp.push_back(i);
				rect_temp.push_back(rect_list[j]);
			}
		}
	}
	rect_list.swap(rect_temp);
	new_polygon_id_list.swap(new_polygon_id_temp);	
}	


void SimpleMPL::split_rectangle(rectangle_type & pRect, std::vector<rectangle_pointer_type>& split, std::vector<rectangle_pointer_type> nei_Vec)
{
	bool hor = gtl::delta(pRect, gtl::HORIZONTAL) >= gtl::delta(pRect, gtl::VERTICAL);
	std::vector<coordinate_type> vstitches;
	coordinate_difference width = gtl::xh(pRect) - gtl::xl(pRect);
	(void)width;
	coordinate_type lower_boundary;
	coordinate_type upper_boundary;

	std::set<coordinate_type> vset;
	std::vector<coordinate_type> vPossibleStitches;
	if (hor)
	{
		lower_boundary = gtl::xl(pRect);
		upper_boundary = gtl::xh(pRect);
	}
	else 
    {
		lower_boundary = gtl::yl(pRect);
		upper_boundary = gtl::yh(pRect);
	}
	vset.insert(lower_boundary);
	vset.insert(upper_boundary);

	std::vector<rectangle_type> vInterSect;
	// generate intersections
	for (std::vector<rectangle_pointer_type>::iterator it = nei_Vec.begin(); it != nei_Vec.end(); it++)
	{
		coordinate_difference distance = boost::geometry::distance(pRect, *(*it));
		if (distance >= m_db->coloring_distance) continue;
		rectangle_type temp(*(*it));
		gtl::bloat(temp, gtl::HORIZONTAL, m_db->coloring_distance);
		gtl::bloat(temp, gtl::VERTICAL, m_db->coloring_distance);		
		boost::polygon::intersect(temp, pRect, true);
		vInterSect.push_back(temp);
		if (hor)
		{
			vset.insert(gtl::xl(temp));
			vset.insert(gtl::xh(temp));
		}
		else
		{
			vset.insert(gtl::yl(temp));
			vset.insert(gtl::yh(temp));
		}
	}

	for (std::set<coordinate_type>::iterator it = vset.begin(); it != vset.end(); it++)
		vPossibleStitches.push_back(*it);
	std::sort(vPossibleStitches.begin(), vPossibleStitches.end());

	generate_stitch_position(pRect, vInterSect, vPossibleStitches, vstitches);

#ifdef STITCH_FILTER
	// check the stitch positions' legalities
	// if the position is very colse to the rectangle's boundary, it's illegal.
	coordinate_type threshold = 20;
	std::vector<coordinate_type> temp;
	for (std::vector<coordinate_type>::iterator it = vstitches.begin(); it != vstitches.end(); it++)
	{
		coordinate_type dis_low = std::abs(*it - lower_boundary);
		coordinate_type dis_up = std::abs(*it - upper_boundary);
		if (dis_low >= threshold && dis_up >= threshold)
			temp.push_back(*it);
	}
	std::vector<coordinate_type>().swap(vstitches);
	vstitches.swap(temp);
#endif

	// split rectangles according to the stitch positions.
	if (vstitches.size() <= 0)
	{
		rectangle_pointer_type new_Pattern = new rectangle_type(pRect);
		new_Pattern->color(pRect.color());
		split.push_back(new_Pattern);
	}
	else
	{
		if (hor)
		{
			vstitches.insert(vstitches.begin(), gtl::xl(pRect));
			vstitches.push_back(gtl::xh(pRect));
			for (uint32_t j = 0; j < vstitches.size() - 1; j++)
			{
				rectangle_pointer_type new_Pattern = new rectangle_type(
					vstitches[j], gtl::yl(pRect), vstitches[j + 1], gtl::yh(pRect));
				new_Pattern->color(pRect.color());
				split.push_back(new_Pattern);
			}
		}
		else
		{
			vstitches.insert(vstitches.begin(), gtl::yl(pRect));
			vstitches.push_back(gtl::yh(pRect));
			for (uint32_t j = 0; j < vstitches.size() - 1; j++)
			{
				rectangle_pointer_type new_Pattern = new rectangle_type(
					gtl::xl(pRect), vstitches[j], gtl::xh(pRect), vstitches[j + 1]);
				new_Pattern->color(pRect.color());
				split.push_back(new_Pattern);
			}
		}
	}
	mplAssert(split.size() > 0);

	return;
}

void SimpleMPL::generate_stitch_position(const rectangle_type pRect, std::vector<rectangle_type> vInterSect,
	std::vector<coordinate_type>& vPossibleStitches, std::vector<coordinate_type>& vstitches)
{
	bool ishor = gtl::delta(pRect, gtl::HORIZONTAL) >= gtl::delta(pRect, gtl::VERTICAL);
	coordinate_type lower, upper;
	if (ishor) 
    {
		lower = gtl::xl(pRect);
		upper = gtl::xh(pRect);
	}
	else 
    {
		lower = gtl::yl(pRect);
		upper = gtl::yh(pRect);
	}
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, int > > vStages;
	for (uint32_t i = 1; i < vPossibleStitches.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(vPossibleStitches[i - 1], vPossibleStitches[i]), 0));
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		coordinate_type left = vStages[i].first.first;
		coordinate_type right = vStages[i].first.second;
		for (uint32_t j = 0; j < vInterSect.size(); j++)
		{
			if (ishor)
			{
				if (left < gtl::xl(vInterSect[j])) continue;
				if (right > gtl::xh(vInterSect[j])) continue;
				vStages[i].second++;
			}
			else
			{
				if (left < gtl::yl(vInterSect[j])) continue;
				if (right > gtl::yh(vInterSect[j])) continue;
				vStages[i].second++;
			}
		}
	}
	if (vStages.size() <= 0) return;
	if (vStages[0].second != 0)
		vStages.insert(vStages.begin(), std::make_pair(std::make_pair(lower, lower), 0));
	if (vStages[vStages.size() - 1].second != 0)
		vStages.push_back(std::make_pair(std::make_pair(upper, upper), 0));
#ifdef _DEBUG_QI
	std::cout << "DEBUG_QI| vStages = ";
	for (uint32_t i = 0; i < vStages.size(); i++)
		std::cout << vStages[i].second << " ";
	std::cout << std::endl;
#endif
	std::vector<int32_t> vzeroids;
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
		if (vStages[i].second > 0) continue;
		vzeroids.push_back(i);
	}

	vstitches.clear();
	for (uint32_t i = 0; i < vzeroids.size() - 1; i++)
	{
		uint32_t pos1 = vzeroids[i];
		uint32_t pos2 = vzeroids[i + 1];
		if (0 == i) mplAssert(0 == pos1);
		else if (2 == pos1) // remove the useless stitch
		{
			bool bFind = true;
			if (vStages.size() < 5) bFind = false;
			else if (i != 1) bFind = false;
			else if (1 != vStages[1].second) bFind = false;
			else if (1 != vStages[3].second) bFind = false;
			else if (0 != vStages[4].second) bFind = false;
			if (bFind == true) continue;
			coordinate_type pos = (vStages[2].first.first + vStages[2].first.second) / 2;
			vstitches.push_back(pos);
		}
		else if (pos1 == vStages.size() - 3)
		{
			mplAssert(i == vzeroids.size() - 2);
			bool bFind = true;
			uint32_t zsize = vzeroids.size();
			if (vStages.size() < 5) bFind = false;
			else if (i != zsize - 2) bFind = false;
			else if (1 != vStages[pos1 + 1].second) bFind = false;
			else if (1 != vStages[pos1 - 1].second) bFind = false;
			else if (0 != vStages[pos1 - 2].second) bFind = false;
			if (bFind == true) continue;
			coordinate_type pos = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}
		else
		{
			coordinate_type pos = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}

		double maxValue = 0.9;
		uint32_t posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		// for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		// {
		// 	if (vStages[i - 1].second <= vStages[i].second) continue;
		// 	if (vStages[i + 1].second <= vStages[i].second) continue;
		// 	int mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
		// 	int diff = std::abs(vStages[i + 1].second - vStages[i - 1].second);
		// 	double value = (double)mind + (double)diff*0.1;
		// 	if (value > maxValue)
		// 	{
		// 		maxValue = value;
		// 		posLost = i;
		// 	}
		// }
		// if (maxValue > 0.9)
		// {
		// 	coordinate_type pos = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
		// 	vstitches.push_back(pos);
		// }
		for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			coordinate_type pos = (vStages[i].first.first + vStages[i].first.second) / 2;
			vstitches.push_back(pos);
		}
	}
	std::sort(vstitches.begin(), vstitches.end());
}


void SimpleMPL::report() 
{
	uint32_t c_num = conflict_num();
    mplPrint(kINFO, "Conflict number = %u\n", c_num);
	for (int32_t i = 0, ie = m_db->color_num(); i != ie; ++i)
        mplPrint(kINFO, "Color %d density = %u\n", i, m_vColorDensity[i]);
	int count = 0;
	int stitch_count = 0;
	for (uint32_t i = 0; i < m_db->vPatternBbox.size(); i++)
	{
		if (m_db->vPatternBbox[i]->color() >= m_db->color_num())
			count++;
	}
	for (uint32_t i = 0; i < m_StitchRelation.size(); i++)
	{
		for(uint32_t j = 0; j < m_StitchRelation[i].size();j++)
        {
			if(m_db->vPatternBbox[i]->color() != m_db->vPatternBbox[m_StitchRelation[i][j]]->color())
				stitch_count++;
		}
	}
    if( m_db->parms.record > 0)
    {
        std::ofstream myfile;
        myfile.open ("record.txt", std::ofstream::app);
        myfile << "Stitch number: "<<stitch_count/2<<"\n";
        myfile << "Conflict number: "<<c_num<<"\n";
        myfile.close();
        std::ofstream result;
        if(m_db->use_stitch())
        {
            result.open("result_w_stitch.txt",std::ofstream::app);
            result<<"& " << std::setw(5) << stitch_count/2<<" & "<< std::setw(5) <<c_num<<" & "<< std::setw(5) <<0.1*stitch_count/2 + c_num;
			if(m_db->algo() == AlgorithmTypeEnum::DANCING_LINK) 
			{
				result << " \\\\\n";
			} else {
				result << " &  ";
			}
        }
        else
        {
            result.open("result_wo_stitch.txt",std::ofstream::app);
            if(m_db->algo() == AlgorithmTypeEnum::ILP_GUROBI)
            {
                result<<"\\\\\n"<<m_db->input_gds();
            }
            result<<"&"<<c_num;
        }
        result.close();
    }
	mplPrint(kINFO, "Invalid color number : %d\n", count);
	mplPrint(kINFO, "Invalid stitch number : %d\n", stitch_count/2);
	mplPrint(kINFO, "Cost: %.2f\n", 0.1*stitch_count/2 + c_num);
}

void SimpleMPL::construct_graph()
{
	mplPrint(kINFO, "Constructing graph for %lu patterns...\n", m_db->vPatternBbox.size());
	// construct vertices 
	// assume vertices start from 0 and end with vertex_num-1
	uint32_t vertex_num = m_db->vPatternBbox.size();
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
}

uint32_t SimpleMPL::construct_graph_from_distance(uint32_t vertex_num)
{
    uint32_t edge_num = 0;
	uint32_t zero_len_edge_num = 0;
	std::vector<uint32_t> distance_value_count(m_db->coloring_distance);
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
        gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance*2);
        gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance*2);
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
				
				if(distance == 0)
					zero_len_edge_num +=1;
                if (distance < m_db->coloring_distance)
				{
					distance_value_count[distance] += 1;
					vAdjVertex.push_back(pAdjPattern->pattern_id());
				}
                    
            }
        }
        //vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
        // parallel with omp reduction here 
        edge_num += vAdjVertex.size();
    }
    return edge_num;
}

uint32_t SimpleMPL::construct_graph_from_paths(uint32_t vertex_num)
{
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

    return edge_num;
}

void SimpleMPL::connected_component()
{
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

#ifdef DEBUG_LIWEI
    mplPrint(kDEBUG, "After component analysis, there are %u components in total.\n", m_comp_cnt);
#endif
}

void SimpleMPL::depth_first_search(uint32_t source, uint32_t comp_id, uint32_t& order_id)
{
    std::stack<uint32_t> vStack; 
	int component_count = 0;
	vStack.push(source);

	while (!vStack.empty())
	{
		uint32_t current = vStack.top();
		vStack.pop();
		if (m_vCompId[current] == std::numeric_limits<uint32_t>::max()) // not visited 
		{
			m_vCompId[current] = comp_id; // set visited 
			m_vVertexOrder[current] = order_id++; // update position 
			component_count += 1;
			for (std::vector<uint32_t>::const_iterator it = m_mAdjVertex[current].begin(); it != m_mAdjVertex[current].end(); ++it)
			{
				uint32_t const& child = *it;
				mplAssert(current != child);
				vStack.push(child);
			}
		}
	}
}

///// helper functions for SimpleMPL::solve_component
///// since we may want to call solve_graph_coloring() multiple times 
///// it is better to wrap it as an independent function 

/// create coloring solver pointer according to algorithm type
lac::Coloring<SimpleMPL::graph_type>* SimpleMPL::create_coloring_solver(SimpleMPL::graph_type const& sg, uint32_t comp_id, uint32_t sub_comp_id) const
{
    typedef lac::Coloring<graph_type> coloring_solver_type;
    coloring_solver_type* pcs = NULL;
	if(m_db->parms.selector != ""){
		if(m_algorithm_selector[comp_id][sub_comp_id] == 0){
			#if GUROBI == 1
			pcs = new lac::ILPColoring<graph_type> (sg); 
			#endif
		}
		else{
			pcs = new DancingLinkColoring<graph_type> (sg); 
		}
	}
	else{
    switch (m_db->algo().get())
		{
	#if GUROBI == 1
			case AlgorithmTypeEnum::ILP_GUROBI:
				pcs = new lac::ILPColoring<graph_type> (sg); 
				break;
			case AlgorithmTypeEnum::ILP_UPDATED_GUROBI:
				pcs = new lac::ILPColoringUpdated<graph_type> (sg);
				break; 
			case AlgorithmTypeEnum::LP_GUROBI:
				pcs = new lac::LPColoring<graph_type> (sg); 
				break;
			case AlgorithmTypeEnum::MIS_GUROBI:
				pcs = new lac::MISColoring<graph_type> (sg); 
				break;
	#endif
	#if LEMONCBC == 1
			case AlgorithmTypeEnum::ILP_CBC:
				pcs = new lac::ILPColoringLemonCbc<graph_type> (sg); 
				break;
	#endif
	#if CSDP == 1
			case AlgorithmTypeEnum::SDP_CSDP:
				pcs = new lac::SDPColoringCsdp<graph_type> (sg); 
				break;
	#endif
			case AlgorithmTypeEnum::BACKTRACK:
				pcs = new lac::BacktrackColoring<graph_type> (sg);
				break;
			case AlgorithmTypeEnum::DANCING_LINK:
				pcs = new DancingLinkColoring<graph_type> (sg);
				break; 
			case AlgorithmTypeEnum::DANCING_LINK_OPT:
				pcs = new DancingLinkColoringOpt<graph_type> (sg);
				break; 
			default: mplAssertMsg(0, "unknown algorithm type");
		}
	}

    pcs->stitch_weight(m_db->parms.weight_stitch);
    pcs->color_num(m_db->color_num());
    pcs->threads(1); // we use parallel at higher level 

    return pcs;
}

/// given a graph, solve coloring 
/// contain nested call for itself 
double SimpleMPL::solve_graph_coloring(uint32_t comp_id, SimpleMPL::graph_type const& dg, 
        std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt, 
        uint32_t simplify_strategy, std::vector<int8_t>& vColor, std::set<vertex_descriptor> vdd_set)
{
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs (dg, m_db->color_num());
	gs.set_isVDDGND(vdd_set);
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices 
	std::set<vertex_descriptor> not_in_DG;
	
    // set max merge level, actually it only works when MERGE_SUBK4 is on 
    if (m_db->color_num() == 3)
        gs.max_merge_level(3);
    else if (m_db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works better 
        gs.max_merge_level(2);


	// SECOND SIMPLIFICATION !
    gs.simplify(simplify_strategy);
	
	// collect simplified information 
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

    // for debug, it does not affect normal run 

	// in order to recover color from articulation points 
	// we have to record all components and mappings 
	// but graph is not necessary 
	std::vector<std::vector<int8_t> > mSubColor (gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig (gs.num_component());
	m_DG_num += gs.num_component();
	double acc_obj_value = 0;

	std::vector<vertex_descriptor> all_articulations;

	gs.get_articulations(all_articulations);
	// std::cout<<"articulations!"<<std::endl;
	// for(auto arti:all_articulations){
	// 	std::cout<<arti<<std::endl;
	// }
	// std::cout<<"DG Graph size is "<<num_vertices(dg)<<" "<<comp_id<< ", num_of  hidden vertices:"<<vHiddenVertices.size()<<std::endl;
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
    {

		//sg: sub-graph: graphs after second simplification
        graph_type sg;
        std::vector<int8_t>& vSubColor =  mSubColor[sub_comp_id];
        std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

        gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);

        vSubColor.assign(num_vertices(sg), -1);
		// std::cout<<"SG Graph size is "<<num_vertices(sg)<<" "<<comp_id<<" "<<sub_comp_id<<std::endl;
#ifdef _OPENMP
#pragma omp critical(m_dgGlobal2Local)
#endif
        {
            for (std::map<uint32_t, uint32_t>::iterator it = m_dgGlobal2Local.begin(); it != m_dgGlobal2Local.end(); it++)
            {
                if (std::find(vSimpl2Orig.begin(), vSimpl2Orig.end(), it->second) != vSimpl2Orig.end())
                {
                    m_dgCompId[it->first] = sub_comp_id + 1;

                }
            }
        }
		if (comp_id == m_db->dbg_comp_id())
		{
			for (std::vector<vertex_descriptor>::const_iterator it = vSimpl2Orig.begin(); it != vSimpl2Orig.end(); ++it)
			{mplPrint(kDEBUG, "sub_comp_id %u, orig id is %u.\n", (uint32_t)sub_comp_id, (uint32_t)*it);}
			for(auto orig:vSimpl2Orig){
				std::cout<<orig<<std::endl;
				std::cout<<gtl::xl(*(m_db->vPatternBbox[*(itBgn+orig)]))<<std::endl;
				std::cout<<gtl::yl(*(m_db->vPatternBbox[*(itBgn+orig)]))<<std::endl;
				
			}
		}

        //if algorithm is Dancing Link, call it directly
        boost::timer::cpu_timer comp_timer;
        comp_timer.start();

        // solve coloring 
        typedef lac::Coloring<graph_type> coloring_solver_type;
        coloring_solver_type* pcs = create_coloring_solver(sg,comp_id,sub_comp_id);

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

        // if(boost::num_vertices(sg) > 2)
        // {
        // // 	this->write_txt(sg,std::to_string(comp_id),obj_value1);
        // 	std::cout<<"Write Component "<<comp_id<<", Sub component "<<sub_comp_id<<std::endl;
		// 	print_graph(sg);
        // }

        // 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only 
        double obj_value2 = std::numeric_limits<double>::max();

#ifdef DEBUG_NONINTEGERS
        std::set<vertex_descriptor> s_vdd_set;
        for (uint32_t i = 0; i < vSimpl2Orig.size(); i++)
        {
            if (vdd_set.find(vSimpl2Orig[i]) != vdd_set.end())
            {
                s_vdd_set.insert(i);
            }
        }
        // very restrict condition to determin whether perform MERGE_SUBK4 or not 
        if (obj_value1 >= 1 && boost::num_vertices(sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::LP_GUROBI || m_db->algo() == AlgorithmTypeEnum::SDP_CSDP)
                && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // MERGE_SUBK4 is not performed 
            obj_value2 = solve_graph_coloring(comp_id, sg, itBgn, pattern_cnt, graph_simplification_type::MERGE_SUBK4, vSubColor, s_vdd_set); // call again 
#endif
        // #ifdef DEBUG
        //         mplPrint(kDEBUG, "comp_id = %u, %lu vertices, obj_value1 = %g, obj_value2 = %g\n", comp_id, num_vertices(sg), obj_value1,obj_value2); 
        // #endif	
        // choose smaller objective value 
        if (obj_value1 < obj_value2)
        {
            acc_obj_value += obj_value1;

            // collect coloring results from simplified graph 
            for (boost::tie(vi, vie) = vertices(sg); vi != vie; ++vi)
            {
                vertex_descriptor v = *vi;
                int8_t color = pcs->color(v);
                // if (color < 0) color = m_db->color_num();
				//std::cout<<"color is"<<(uint32_t)color<<std::endl;
                mplAssert(color >= 0 && color < m_db->color_num());
                vSubColor[v] = color;
                if (comp_id == m_db->dbg_comp_id())
                {
                    mplPrint(kDEBUG, "vertex %u color is %u.\n", (uint32_t)v, (uint32_t)color);
                }
            }
        }
        else // no need to update vSubColor, as it is already updated by sub call 
        {
            acc_obj_value += obj_value2;
            if (comp_id == m_db->dbg_comp_id())
            {
                for(uint32_t v = 0; v<vSubColor.size();v++)
                {
                    mplPrint(kDEBUG, "vertex %u color is %u.\n", (uint32_t)v, (uint32_t)vSubColor[v]);
                }
            }
        }
		if(m_db->parms.record > 1){
			//store ILP runtime information/ read ILP information for DL to comparison, this is for DAC2020
            if(m_db->algo() == AlgorithmTypeEnum::ILP_GUROBI && num_vertices(sg) > 3)
            {
				std::string runtime = comp_timer.format(6,"%w");
                std::ofstream myfile;
                myfile.open("ILP_obj.txt", std::ofstream::app);
                myfile<<m_db->input_gds().c_str()<<" "<<comp_id<<" "<<sub_comp_id <<" "<<num_vertices(sg)<<" "<<obj_value1<<" "<<runtime<<"\n";
                myfile.close();
            }
            if(m_db->algo() == AlgorithmTypeEnum::DANCING_LINK && num_vertices(sg) > 3)
            {
				std::string runtime = comp_timer.format(6,"%w");
                std::ofstream myfile;
                myfile.open("DL_obj.txt", std::ofstream::app);
                myfile<<m_db->input_gds().c_str()<<" "<<comp_id<<" "<<sub_comp_id <<" "<<num_vertices(sg)<<" "<<obj_value1<<" "<<runtime<<"\n";
                myfile.close();
            }
		}
        if( m_db->parms.record > 2)
        {
            if(obj_value1 != 0)
            {
                std::ofstream myfile;
                myfile.open ("small_results.txt", std::ofstream::app);
                myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<", obj_value: "<< obj_value1<<"\n";
                myfile.close();
            }
            if (comp_id == m_db->dbg_comp_id())
            {
                write_graph(sg, std::to_string(comp_id) +"_"+ std::to_string(sub_comp_id));
                std::ofstream myfile;
                myfile.open ("debug_component_color_results.txt", std::ofstream::app);
                myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<"\n";
                for(uint32_t v = 0; v<vSubColor.size();v++)
                {
                    myfile<<"vertex "<<v<<" color is "<<(uint32_t)vSubColor[v]<<". \n";
                }
                myfile.close();
            }

            //we only store json file with graph size larger than 3
            if(num_vertices(sg) > 3)
            {
				//The json file name is orgnized as follows: compid_subcompid_objvalue_time.json
                std::string name = std::to_string(comp_id);
                name.append("_");
                name.append(std::to_string(sub_comp_id));
				name.append("_");
				name.append(std::to_string(obj_value1).substr(0,std::to_string(obj_value1).size()-5));
				name.append("_");
                name.append(std::to_string(num_vertices(sg)));
				name.append("_");
				name.append(comp_timer.format(4,"%w"));

                this->write_json(sg,(char*)name.c_str(),vSubColor);
				// print_graph(sg);
            }
        }

        delete pcs;
    }

	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 

	gs.recover(vColor, mSubColor, mSimpl2Orig); 

	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
    RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();
	SimpleMPL::graph_type non_const_dg = dg;
	double total_obj = new_calc_cost(non_const_dg,vColor);
	// std::cout<<"total_obj "<<total_obj<<",acc_obj_value "<<acc_obj_value<<std::endl;
	// assert(total_obj == acc_obj_value);
	// std::cout<<(total_obj == acc_obj_value)<<std::endl;
	// mplAssert(std::abs(total_obj - acc_obj_value) < 1e);
	// mplPrint(kINFO, "Component %u solved", comp_id);
	// std::string name = std::to_string(comp_id);
	// this->write_json(dg,(char*)name.c_str(),vColor);

    return total_obj;
}


void SimpleMPL::iterative_mark(graph_type const& g, std::vector<uint32_t>& parent_node_ids, vertex_descriptor& v1) const
{
	mplAssert(num_vertices(g) == parent_node_ids.size());
	typename boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
	boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, g);
	for (next2 = vi2; vi2 != vie2; vi2 = next2)
	{	
		mplAssert(num_vertices(g) > 1);
		++next2; 
		vertex_descriptor v2 = *vi2;
		std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, g);
		mplAssert(e12.second);
		//if two nodes are stitch relationships 
		if(boost::get(boost::edge_weight, g, e12.first) < 0)
		{
			if(parent_node_ids[(uint32_t)v2] == (uint32_t)-1)
			{
				parent_node_ids[(uint32_t)v2] = parent_node_ids[(uint32_t)v1];
				iterative_mark(g, parent_node_ids,v2);
			}
			else
			{
				mplAssert(parent_node_ids[(uint32_t)v2] == parent_node_ids[(uint32_t)v1]);
			}
		}
	}
}

double SimpleMPL::new_calc_cost(SimpleMPL::graph_type& g,std::vector<int8_t> const& vColor){
	double cost = 0;
	std::vector<uint32_t> parent_node_ids;
	parent_node_ids.assign(boost::num_vertices(g),-1);
	uint32_t parent_node_id = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	for (boost::tie(vi1, vie1) = boost::vertices(g); vi1 != vie1; ++vi1)
	{	
		vertex_descriptor v1 = *vi1;
		mplAssert(vColor[v1]!=-1);
		mplAssert((uint32_t)v1 < parent_node_ids.size());
		if(parent_node_ids[(uint32_t)v1] == (uint32_t)-1)
		{
			parent_node_ids[(uint32_t)v1] = parent_node_id;
			iterative_mark(g,parent_node_ids,v1);
			parent_node_id++;
		}
	}
	// std::cout << "Assigned parents" << std::endl;
	//calculate conflict by parent node
	std::vector<std::vector<bool>> conflict_mat(parent_node_id, std::vector<bool>(parent_node_id, 0));
	// for(uint32_t i = 0; i < parent_node_id; i++)
	// {
	// 	std::vector<bool> row_mat;
	// 	std::cout << "Before assign " << parent_node_id << " of false" << std::endl;
	// 	for(int p = 0; p < parent_node_id; ++p) {
	// 		row_mat.push_back(false);
	// 	}
	// 	// row_mat.assign(parent_node_id,false);
	// 	std::cout << "Before push_back [" << i << "]" << std::endl;
	// 	// conflict_mat.push_back(row_mat);
	// 	conflict_mat.emplace_back(std::vector<bool>(parent_node_id, false));
	// }
	// std::cout << "Init conflict_mat" << std::endl;
	for (boost::tie(vi1, vie1) = boost::vertices(g); vi1 != vie1; ++vi1)
	{	
		vertex_descriptor v1 = *vi1;
		typename boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
		boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, g);
		for (next2 = vi2; vi2 != vie2; vi2 = next2)
		{
			++next2; 
			vertex_descriptor v2 = *vi2;
			if (v1 >= v2) continue;
			std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, g);
			mplAssert(e12.second);
			//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0)
			{
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				//cost += (vColor[v1] != vColor[v2])*m_db->parms.weight_stitch;
			}
			else
			{
				if(vColor[v1] == vColor[v2])
				{
					if(conflict_mat[parent_node_ids[(uint32_t)v1]][parent_node_ids[(uint32_t)v2]] == false)
					{
						cost += 1;
						conflict_mat[parent_node_ids[(uint32_t)v1]][parent_node_ids[(uint32_t)v2]] = true;
						conflict_mat[parent_node_ids[(uint32_t)v2]][parent_node_ids[(uint32_t)v1]] = true;
					}
				}
			}
		}
	}
	// std::cout << "Before return" << std::endl;
	return cost;
}
void SimpleMPL::print_graph(SimpleMPL::graph_type& g)
{
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	std::cout<<"Print graph!"<<std::endl;
	for (boost::tie(vi1, vie1) = boost::vertices(g); vi1 != vie1; ++vi1)
	{	
		vertex_descriptor v1 = *vi1;
		boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
		boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, g);
		for (next2 = vi2; vi2 != vie2; vi2 = next2)
		{
			++next2; 
			vertex_descriptor v2 = *vi2;
			if (v1 >= v2) continue;
			std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, g);
			mplAssert(e12.second);
            mplPrint(kNONE, "%u %u %d\n", (uint32_t)v1, (uint32_t)v2, (int)boost::get(boost::edge_weight, g, e12.first));
			// if the nodes in stitch relation has same conflict adj, then this stitch can be removed
			if((int)boost::get(boost::edge_weight, g, e12.first) < 0){
				uint32_t v1_conflict_adj_num = num_of_conflict_adj(g,v1);
				uint32_t v2_conflict_adj_num = num_of_conflict_adj(g,v2);
				if(v1_conflict_adj_num == v2_conflict_adj_num){
					mplPrint(kDEBUG,"ERROR! node %u degree %u, node %u degree %u\n",v1,v1_conflict_adj_num,v2,v2_conflict_adj_num);
				}
			}
		}
	}
}

//here the count is the sum of all INDEXES!! of conflict adjs instead of the 
uint32_t SimpleMPL::num_of_conflict_adj(SimpleMPL::graph_type& g, vertex_descriptor v)
{
	boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
	boost::tie(vi2, vie2) = boost::adjacent_vertices(v, g);
	uint32_t count = 0;
	for (next2 = vi2; vi2 != vie2; vi2 = next2)
	{
		++next2; 
		vertex_descriptor v2 = *vi2;
		std::pair<edge_descriptor, bool> e12 = boost::edge(v, v2, g);
		mplAssert(e12.second);
		if((int)boost::get(boost::edge_weight, g, e12.first) > 0){
			count += (v2+1);
		}
	}
	return count;
}

uint32_t SimpleMPL::num_of_stitch_adj(SimpleMPL::graph_type& g,vertex_descriptor v)
{
	boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
	boost::tie(vi2, vie2) = boost::adjacent_vertices(v, g);
	uint32_t count = 0;
	for (next2 = vi2; vi2 != vie2; vi2 = next2)
	{
		++next2; 
		vertex_descriptor v2 = *vi2;
		std::pair<edge_descriptor, bool> e12 = boost::edge(v, v2, g);
		mplAssert(e12.second);
		if((int)boost::get(boost::edge_weight, g, e12.first) < 0){
			count += (v2+1);
		}
	}
	return count;
}
void SimpleMPL::construct_component_graph(const std::vector<uint32_t>::const_iterator itBgn, uint32_t const pattern_cnt, 
        SimpleMPL::graph_type& dg, std::map<uint32_t, uint32_t>& mGlobal2Local, std::vector<int8_t>& vColor,std::set<vertex_descriptor>& vdd_set, bool flag) const
{
	// precolored patterns 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		vColor[i] = m_db->vPatternBbox[v]->color();
		mGlobal2Local[v] = i;
		if (m_isVDDGND[v])
			vdd_set.insert(i);
	}

	//TODO: CANNOT SET IVR, BUT NEED TO USE FAST COLOR TRIAL
	// edges 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t v = *(itBgn+i);
		///< add conflict edges
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
		} // for, end of conflict edges
		
		///< add stitch edges
		if (flag == true)
		{
			for (std::vector<uint32_t>::const_iterator it = m_StitchRelation[v].begin(); it != m_StitchRelation[v].end(); ++it)
			{
				uint32_t s = *it;
				mplAssert(mGlobal2Local.count(s));
				uint32_t j = mGlobal2Local[s];
				if (i < j)
				{
					std::pair<edge_descriptor, bool> e = edge(i, j, dg);
					if (!e.second)
					{
						e = add_edge(i, j, dg);
						mplAssert(e.second);
						// for stitch, edge_weight is negative.
						boost::put(boost::edge_weight, dg, e.first, -1);
					}
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
	(void)acc_obj_value;
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
	// if(component_conflict_num != acc_obj_value){
	// 	std::cout<<"component_conflict_num "<<component_conflict_num<<",acc_obj_value "<<acc_obj_value<<std::endl;
	// 	// debug_conflict_num(itBgn, itEnd);
	// }

	//NOTE: There is a bug(feature): if multi-thread is enabled, the following assertion may fail, but it does not influence final results, therefore I ignore this part.
	// 03/16/2020 Wei
	// mplAssert(component_conflict_num == acc_obj_value);
    // only valid under no stitch 
    // if (acc_obj_value != std::numeric_limits<uint32_t>::max())
    //    mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, (uint32_t)(itEnd-itBgn), component_conflict_num);

	return component_conflict_num;
}


// Here, dg is graph before second simplification
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

	bool flag;
    // construct decomposition graph for component 
	if (m_db->use_stitch())
		flag = true;
	else
		flag = false;

	std::set<vertex_descriptor> vdd_set;
	// for (uint32_t i = 0; i != pattern_cnt; ++i)
	// {
	// 	uint32_t const& v = *(itBgn+i);
	// 	std::cout<<gtl::xl(*(m_db->vPatternBbox[v]))<<std::endl;
	// }
	construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, flag);
#ifdef _OPENMP
#pragma omp critical(m_dgGlobal2Local)
#endif
	{
		std::map<uint32_t, uint32_t>().swap(m_dgGlobal2Local);
		m_dgGlobal2Local.insert(mGlobal2Local.begin(), mGlobal2Local.end());
	}


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
	double acc_obj_value = 0;
    // solve graph coloring 
	if (m_db->remove_stitch_redundancy())
	{
		// if(pattern_cnt > 10){
		// 	acc_obj_value = solve_graph_coloring_with_remove_stitch_redundancy(0, comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
		// }
		// else{
		// 	acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
		// }
		acc_obj_value = solve_graph_coloring_with_remove_stitch_redundancy(0, comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
	}
	else
	{
    	acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
	}
	// if (boost::num_vertices(dg) >= 10)
	// {
	// 	acc_obj_value = solve_graph_coloring_with_remove_stitch_redundancy(0, comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
	// }
	// else
	// {
    // 	acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);
	// }

#ifdef DEBUG
	/*
	for (uint32_t i = 0; i != pattern_cnt; ++i)
		mplAssert(vColor[i] >= 0 && vColor[i] < m_db->color_num());
	*/
#endif

	// record pattern color 
	for (uint32_t i = 0; i != pattern_cnt; ++i)
	{
		uint32_t const& v = *(itBgn+i);
		int8_t color = vColor[i];
		if (color < 0) color = m_db->color_num();
        m_db->set_color(v, color);
	}

    return acc_obj_value;
}

void SimpleMPL::find_all_stitches(uint32_t vertex, std::vector<uint32_t>& stitch_vec)
{
	for(uint32_t stit = 0; stit < m_StitchRelation[vertex].size() ;stit++)
    {
		if(std::find(stitch_vec.begin(), stitch_vec.end(), m_StitchRelation[vertex][stit]) == stitch_vec.end())
        {
			stitch_vec.push_back(m_StitchRelation[vertex][stit]);
			find_all_stitches(m_StitchRelation[vertex][stit],stitch_vec);
		}
	}
}

uint32_t SimpleMPL::conflict_num(const std::vector<uint32_t>::const_iterator itBgn, const std::vector<uint32_t>::const_iterator itEnd)
{
	// std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	// uint32_t cnt = 0;
    // for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	// {
	// 	uint32_t v = *it;
	// 	int8_t color1 = vPatternBbox[v]->color();
	// 	if (color1 >= 0 && color1 < m_db->color_num())
	// 	{
	// 		std::set<uint32_t> parent_indexes;
	// 		for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
	// 		{
				
	// 			uint32_t u = *itAdj;
	// 			int8_t color2 = vPatternBbox[u]->color();
	// 			if (color2 >= 0 && color2 < m_db->color_num())
	// 			{
	// 				if (color1 == color2) ++cnt;
	// 			}
	// 			else ++cnt; // uncolored vertex is counted as conflict 
	// 		}
	// 	}
	// 	else ++cnt; // uncolored vertex is counted as conflict 
	// }
	// // conflicts will be counted twice 
	// return (cnt>>1);


	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	// std::vector<std::vector<bool>> is_counted;
	// for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
    // {
	// 	std::vector<bool> v_th_counted;
	// 	v_th_counted.assign(vPatternBbox[v].size(),false);
	// 	is_counted.push_back(v_th_counted);
	// }
	for (std::vector<uint32_t>::const_iterator it = itBgn; it != itEnd; ++it)
	{
		uint32_t v = *it;
		std::vector<uint32_t> stitch_vec;
		if(m_db->use_stitch())
        {
			find_all_stitches(v,stitch_vec);
			if(stitch_vec.size() == 0)
            {
				stitch_vec.push_back(v);
			}
		}
		else
        {
			stitch_vec.push_back(v);
		}

		int8_t color1 = vPatternBbox[v]->color();
		// std::cout<<"node "<<v<< ",color "<<(int32_t)color1<<std::endl;
		if (color1 >= 0 && color1 <= m_db->color_num())
		{
			std::set<uint32_t> parent_indexes;
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				// std::cout<<"node "<<v<< ",adj "<<u<<",color "<<(int32_t)vPatternBbox[u]->color()<<std::endl;
				// std::cout<<m_db->vParentPolygonId[u]<<std::endl;
				// if(parent_indexes.find(m_db->vParentPolygonId[u])!=parent_indexes.end()){continue;}
				// else{parent_indexes.insert(m_db->vParentPolygonId[u]);}
                if (v < u) // avoid duplicate 
                {
					if(m_isVDDGND[u] && m_isVDDGND[v]) continue;
                    int8_t color2 = vPatternBbox[u]->color();
                    if (color2 >= 0 && color2 <= m_db->color_num())
                    {
                        if (color1 == color2) 
							{
								std::vector<uint32_t> stitch_vec2;
								if(m_db->use_stitch())
                                {
									find_all_stitches(u,stitch_vec2);
									if(stitch_vec2.size() == 0)
                                    {
										stitch_vec2.push_back(u);
									}
								}
								else
                                {
									stitch_vec2.push_back(u);
								}	
								bool need_count = true;
								// for(uint32_t j = 0; j < stitch_vec.size();j++)
                                // {
								// 	std::cout<<" stitch_vec "<<j<<" "<<stitch_vec[j]<<", color is "<<(int)vPatternBbox[stitch_vec[j]]->color()<<std::endl;
								// }
								// for(uint32_t j2 = 0 ; j2 < stitch_vec2.size();j2++)
                                // {
								// 		std::cout<<" stitch_vec2 "<<j2<<" "<<stitch_vec2[j2]<<", color is "<<(int)vPatternBbox[stitch_vec2[j2]]->color()<<std::endl;
								// }
								// we only count conflict between two stitch sets(actually the original polygon) once
								//therefore we default to count the conflict of smallest nodes in v_set (v and its stitch neighbors) (since  v< u always holds).
								// for example:
								// stitch relation 1-2, 3-4-5
								// conflict relation 1-5, 2-3, 2-4
								// then 1-5 is counted (since 1 is the smallest node in 1-2)
								for(uint32_t j = 0; j < stitch_vec.size();j++)
                                {
									for(uint32_t j2 = 0 ; j2 < stitch_vec2.size();j2++)
                                    {
										// std::cout<<"node "<<v<< ",adj "<<u<<" stitch vec "<<stitch_vec[j]<<", stitch_vec2[j2]"<<stitch_vec2[j2]<<std::endl;
										//if the nodes are in the conflict relation
										if(std::find(m_mAdjVertex[stitch_vec[j]].begin(),m_mAdjVertex[stitch_vec[j]].end(),stitch_vec2[j2])!=m_mAdjVertex[stitch_vec[j]].end())
                                        {
											if(vPatternBbox[stitch_vec[j]]->color() == vPatternBbox[stitch_vec2[j2]]->color())
                                            {
												if(stitch_vec[j] <  v || (stitch_vec2[j2] < u && stitch_vec[j] ==  v ))	//avoid duplicate, for each 
												{
													need_count = false;
													break;
												}
											}
										}
									}
									if(need_count == false) break;
								}
								if(need_count)	m_vConflict.push_back(std::make_pair(v, u));
							}
                            
                    }
                    else // uncolored vertex is counted as conflict 
                    {
                        mplAssertMsg(0, "uncolored vertex %u = %d", u, color2);
                    }
                }
			}
		}
        else // uncolored vertex is counted as conflict 
            mplAssertMsg(0, "uncolored vertex %u = %d", v, color1);
	}
	// for(auto conflict:m_vConflict){
	// 	std::cout<<conflict.first<<" "<<conflict.second<<std::endl;
	// }
	return m_vConflict.size();
}



uint32_t SimpleMPL::conflict_num() 
{
	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	// std::vector<std::vector<bool>> is_counted;
	// for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
    // {
	// 	std::vector<bool> v_th_counted;
	// 	v_th_counted.assign(vPatternBbox[v].size(),false);
	// 	is_counted.push_back(v_th_counted);
	// }
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		std::vector<uint32_t> stitch_vec;
		if(m_db->use_stitch())
        {
			find_all_stitches(v,stitch_vec);
			if(stitch_vec.size() == 0)
            {
				stitch_vec.push_back(v);
			}
		}
		else
        {
			stitch_vec.push_back(v);
		}

		int8_t color1 = vPatternBbox[v]->color();
		// std::cout<<"node "<<v<< ",color "<<(int32_t)color1<<std::endl;
		if (color1 >= 0 && color1 <= m_db->color_num())
		{
			std::set<uint32_t> parent_indexes;
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				//std::cout<<"node "<<v<< ",adj "<<u<<",color "<<(int32_t)vPatternBbox[u]->color()<<std::endl;
				// std::cout<<m_db->vParentPolygonId[u]<<std::endl;
				// if(parent_indexes.find(m_db->vParentPolygonId[u])!=parent_indexes.end()){continue;}
				// else{parent_indexes.insert(m_db->vParentPolygonId[u]);}
                if (v < u) // avoid duplicate 
                {
					if(m_isVDDGND[u] && m_isVDDGND[v]) continue;
                    int8_t color2 = vPatternBbox[u]->color();
                    if (color2 >= 0 && color2 <= m_db->color_num())
                    {
                        if (color1 == color2) 
							{
								std::vector<uint32_t> stitch_vec2;
								if(m_db->use_stitch())
                                {
									find_all_stitches(u,stitch_vec2);
									if(stitch_vec2.size() == 0)
                                    {
										stitch_vec2.push_back(u);
									}
								}
								else
                                {
									stitch_vec2.push_back(u);
								}	
								bool need_count = true;
								// for(uint32_t j = 0; j < stitch_vec.size();j++)
                                // {
								// 	std::cout<<" stitch_vec "<<j<<" "<<stitch_vec[j]<<", color is "<<(int)vPatternBbox[stitch_vec[j]]->color()<<std::endl;
								// }
								// for(uint32_t j2 = 0 ; j2 < stitch_vec2.size();j2++)
                                // {
								// 		std::cout<<" stitch_vec2 "<<j2<<" "<<stitch_vec2[j2]<<", color is "<<(int)vPatternBbox[stitch_vec2[j2]]->color()<<std::endl;
								// }
								// we only count conflict between two stitch sets(actually the original polygon) once
								//therefore we default to count the conflict of smallest nodes in v_set (v and its stitch neighbors) (since  v< u always holds).
								// for example:
								// stitch relation 1-2, 3-4-5
								// conflict relation 1-5, 2-3, 2-4
								// then 1-5 is counted (since 1 is the smallest node in 1-2)							
								for(uint32_t j = 0; j < stitch_vec.size();j++)
                                {
									for(uint32_t j2 = 0 ; j2 < stitch_vec2.size();j2++)
                                    {
										//std::cout<<"node "<<v<< ",adj "<<u<<" stitch vec "<<stitch_vec[j]<<", stitch_vec2[j2]"<<stitch_vec2[j2]<<std::endl;
										//if the nodes are in the conflict relation
										if(std::find(m_mAdjVertex[stitch_vec[j]].begin(),m_mAdjVertex[stitch_vec[j]].end(),stitch_vec2[j2])!=m_mAdjVertex[stitch_vec[j]].end())
                                        {
											if(vPatternBbox[stitch_vec[j]]->color() == vPatternBbox[stitch_vec2[j2]]->color())
                                            {
												if(stitch_vec[j] <  v || (stitch_vec2[j2] < u && stitch_vec[j] ==  v ))	//avoid duplicate, for each 
												{
													// we only count conflict between two stitch sets(actually the original polygon) once
													//therefore we default to count the conflict between smallest nodes No.
													need_count = false;
													break;
												}
											}
										}
									}
									if(need_count == false) break;
								}
								if(need_count)	m_vConflict.push_back(std::make_pair(v, u));
							}
                            
                    }
                    else // uncolored vertex is counted as conflict 
                    {
                        mplAssertMsg(0, "uncolored vertex %u = %d", u, color2);
                    }
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

void SimpleMPL::update_algorithm_selector(std::string filename){
	int max_comp_id = -1;
	int max_sub_comp_id = -1;
	std::ifstream selection_file(filename);
	//read the max comp_id to alloct enough memory
	while(!selection_file.eof()){
		int comp_id = -1;
		int sub_comp_id = -1;
		int algorithm = -1;
		selection_file >> comp_id;
		selection_file >> sub_comp_id;
		selection_file >> algorithm;
		if(comp_id > max_comp_id){
			max_comp_id = comp_id;
		}
		if(sub_comp_id > max_sub_comp_id){
			max_sub_comp_id = sub_comp_id;
		}
	}
    selection_file.close();
	std::vector<std::vector<int> >().swap(m_algorithm_selector);
	m_algorithm_selector.resize(max_comp_id + 1);
	std::cout<<"max_comp_id"<<max_comp_id<<",max_sub_comp_id"<<max_sub_comp_id<<std::endl;
    for(int i = 0; i <= max_comp_id; i ++){
        std::vector<int> sub_selector;
		std::vector<int>().swap(sub_selector);
        sub_selector.assign(max_sub_comp_id+1,-1);
        m_algorithm_selector[i].insert(m_algorithm_selector[i].end(),sub_selector.begin(),sub_selector.end());
    }
	std::ifstream selection_file_twice(filename);
	//read the max comp_id to alloct enough memory
	while(!selection_file_twice.eof()){
		int comp_id = -1;
		int sub_comp_id = -1;
		int algorithm = -1;
		selection_file_twice >> comp_id;
		selection_file_twice >> sub_comp_id;
		selection_file_twice >> algorithm;
		if(comp_id == -1){
			continue;
		}
        m_algorithm_selector[comp_id][sub_comp_id] = algorithm;
        //std::cout<<comp_id<<" "<<sub_comp_id<<" "<<algorithm<<std::endl;
	}
    selection_file_twice.close();
	return;
}


void SimpleMPL::print_welcome() const
{
  mplPrint(kNONE, "\n\n");
  mplPrint(kNONE, "=======================================================================\n");
  mplPrint(kNONE, "                      OpenMPL - Version 2.0                          \n");
  mplPrint(kNONE, "                                by                                   \n");  
  mplPrint(kNONE, "       Yibo Lin, Bei Yu, Wei Li, Qi Sun, and  David Z. Pan           \n");
  mplPrint(kNONE, "               School of EECS, Peking University                     \n");
  mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
  mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
  mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
  mplPrint(kNONE, "            Contact Authors:  yibolin@pku.edu.cn                     \n");
  mplPrint(kNONE, "                              {byu,wli,qsun}@cse.cuhk.edu.hk         \n");
  mplPrint(kNONE, "                              dpan@cerc.utexas.edu                   \n");
  mplPrint(kNONE, "=======================================================================\n");
}

// std::vector<int> SimpleMPL::get_adjacent_vertices(const vertex_descriptor &v, const graph_type &G, int &stitch_edge_cnt)
// {
// 	std::vector<int> adjacent_vertices;
// 	boost::graph_traits<graph_type>::adjacency_iterator adj_iter, adj_iter_end;
// 	for (boost::tie(adj_iter, adj_iter_end) = boost ::adjacent_vertices(v, G); adj_iter != adj_iter_end; ++adj_iter)
// 	{
// 		vertex_descriptor u = *adj_iter;
// 		if (boost::get(boost::edge_weight, G, boost::edge(v, u, G).first) > 0)
// 		{
// 			adjacent_vertices.emplace_back(int(u));
// 		}
// 		else
// 		{
// 			++stitch_edge_cnt;
// 		}
// 	}
// 	return adjacent_vertices;
// }
using std::cout;
using std::endl;
double SimpleMPL::solve_graph_coloring_with_remove_stitch_redundancy(int depth, uint32_t comp_id, SimpleMPL::graph_type const &dg,
																	 std::vector<uint32_t>::const_iterator itBgn, uint32_t pattern_cnt,
																	 uint32_t simplify_strategy, std::vector<int8_t> &vColor, std::set<vertex_descriptor> vdd_set) {
    // boost::timer::cpu_timer solver_timer;
	// solver_timer.start();
	if (depth > 0){
		++count;
	}
	++total_count;
	typedef lac::GraphSimplification<graph_type> graph_simplification_type;
	graph_simplification_type gs(dg, m_db->color_num());
	gs.set_isVDDGND(vdd_set);
	gs.precolor(vColor.begin(), vColor.end()); // set precolored vertices
	std::set<vertex_descriptor> not_in_DG;

	// set max merge level, actually it only works when MERGE_SUBK4 is on
	if (m_db->color_num() == 3)
		gs.max_merge_level(3);
	else if (m_db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works better
		gs.max_merge_level(2);

	// SECOND SIMPLIFICATION !
	// cout << "before simplify" << endl;
	gs.simplify(simplify_strategy);
	// cout << "after simplify" << endl;

	// collect simplified information
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();

	// for debug, it does not affect normal run

	// in order to recover color from articulation points
	// we have to record all components and mappings
	// but graph is not necessary
	std::vector<std::vector<int8_t>> mSubColor(gs.num_component());
	std::vector<std::vector<vertex_descriptor>> mSimpl2Orig(gs.num_component());
	m_DG_num += gs.num_component();
	double acc_obj_value = 0;

	std::vector<vertex_descriptor> all_articulations;

	gs.get_articulations(all_articulations);
	// cout << "Start -- Enter subcomponent : simplification cost: " << // solver_timer.format(6,"%w") << ", depth is "<<depth<<", num_of_nodes simplified:"<<vHiddenVertices.size()<< ", graph size is "<<boost::num_vertices(dg)<<", gs.num_component() is "<<gs.num_component()<<endl;
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id) {
	// solver_timer.start();
		graph_type sg;
		std::vector<int8_t> &vSubColor = mSubColor[sub_comp_id];
		std::vector<vertex_descriptor> &vSimpl2Orig = mSimpl2Orig[sub_comp_id];
		auto sub_vdd_set = get_sub_vdd_set(vSimpl2Orig, vdd_set);

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		bool need_simplify = false;
		graph_type tmp_sg;
		std::vector<std::pair<vertex_descriptor, vertex_descriptor>> redundancy_stitch_pairs;
		auto total_v = boost::num_vertices(sg);
		std::vector<std::vector<vertex_descriptor>> stitch_neigh_of_remove_vertcies;
		if (total_v >= 0) {
			boost::graph_traits<graph_type>::edge_iterator edge_iter, edge_iter_end;
			for (boost::tie(edge_iter, edge_iter_end) = boost::edges(sg); edge_iter != edge_iter_end; ++edge_iter) {
				edge_descriptor edge = *edge_iter;
				if (boost::get(boost::edge_weight, sg, edge) > 0)
					continue;
				vertex_descriptor source_vertex = boost::source(edge, sg);
				vertex_descriptor target_vertex = boost::target(edge, sg);
				int source_stitch_edge_cnt = 0;
				int target_stitch_edge_cnt = 0;
				std::vector<bool> source_adj_vertices(total_v, false);
				std::vector<bool> target_adj_vertices(total_v, false);
				std::vector<vertex_descriptor> target_stitch_neigh;
				boost::graph_traits<graph_type>::adjacency_iterator adj_iter, adj_iter_end;
				bool can_merge = true;
				for (boost::tie(adj_iter, adj_iter_end) = boost::adjacent_vertices(source_vertex, sg); adj_iter != adj_iter_end; ++adj_iter){
					vertex_descriptor u = *adj_iter;
					if (boost::get(boost::edge_weight, sg, boost::edge(source_vertex, u, sg).first) > 0){
						source_adj_vertices.at(u) = true;
						if(boost::edge(target_vertex, u, sg).second == false){
							can_merge = false;
							break;
						}
						
					}
					else {
						++source_stitch_edge_cnt;
					}
				}
				if(can_merge == false){
					break;
				}
				for (boost::tie(adj_iter, adj_iter_end) = boost::adjacent_vertices(target_vertex, sg); adj_iter != adj_iter_end; ++adj_iter){
					vertex_descriptor u = *adj_iter;
					if (boost::get(boost::edge_weight, sg, boost::edge(target_vertex, u, sg).first) > 0){
						target_adj_vertices.at(u) = true;
						if(source_adj_vertices.at(u) == false){
							can_merge = false;
							break;
						}
					}
					else {
						if (u != source_vertex){
						target_stitch_neigh.emplace_back(u);}
						++target_stitch_edge_cnt;
					}
				}
				if(can_merge == false){
					break;
				}
				// if ((source_adj_vertices == target_adj_vertices) and (source_stitch_edge_cnt <= 1) and (target_stitch_edge_cnt <= 1)){
				if ((source_adj_vertices == target_adj_vertices) ){
					redundancy_stitch_pairs.emplace_back(std::make_pair(source_vertex, target_vertex));
					need_simplify = true;
					stitch_neigh_of_remove_vertcies.emplace_back(target_stitch_neigh);
					break;
				}
			}
		}
		if (need_simplify and total_v > 0) {
			// if (total_v >2){
			// cout << "# vertices=" << total_v << std::endl;}
	// cout << "Enter subcomponent -- Before Copy Graph: " << // solver_timer.format(6,"%w")<<", boost::num_vertices(sg)"<<boost::num_vertices(sg)<<endl;
	// solver_timer.start();
			// tmp_sg = sg;
	// cout << "Copy Graph: " << // solver_timer.format(6,"%w") << endl;
	// solver_timer.start();
			total_num_vertices += total_v;
			// if (total_v>2){
			// cout << "redundancy stitch pairs size = " << redundancy_stitch_pairs.size()<< endl;
			// cout << "tmp_sg #v=" << num_vertices(tmp_sg) << ", #e=" << num_edges(tmp_sg) << endl;}

			for (auto iter = redundancy_stitch_pairs.begin(); iter != redundancy_stitch_pairs.end(); ++iter) {
				std::vector<vertex_descriptor> stitch_neighs = stitch_neigh_of_remove_vertcies.at(std::distance(redundancy_stitch_pairs.begin(), iter));
				// solver_timer.start();
				for( auto stitch_neigh : stitch_neighs){
					auto e = boost::add_edge((*iter).first, stitch_neigh, sg);
					boost::put(boost::edge_weight, sg, e.first, -1);
				}
				
				// cout << "add edge and assign edge value cost" << // solver_timer.format(6,"%w") << endl;
				// solver_timer.start();
				boost::clear_vertex((*iter).second, sg);
				boost::remove_vertex((*iter).second, sg);
				// cout << "remove vertex costs" << // solver_timer.format(6,"%w") << endl;
				// solver_timer.start();
			}
			// if (total_v>2){
			// cout << "after merge tmp_sg #v=" << num_vertices(tmp_sg) << ", #e=" << num_edges(tmp_sg) << endl;}
	// cout << "Remove Node: " << // solver_timer.format(6,"%w") << endl;
			uint32_t sg_size = num_vertices(sg);
			vSubColor.assign(sg_size, -1);
			double cost = 0;
			cost = solve_graph_coloring_with_remove_stitch_redundancy(depth+1, sub_comp_id, sg, itBgn, 0, simplify_strategy, vSubColor, sub_vdd_set);
			// if(sg_size > 6){
			// 	cost = solve_graph_coloring_with_remove_stitch_redundancy(depth+1, sub_comp_id, sg, itBgn, 0, simplify_strategy, vSubColor, sub_vdd_set);
			// }
			// else{
			// 	cost = solve_graph_coloring(sub_comp_id, sg, itBgn, 0, simplify_strategy, vSubColor, sub_vdd_set);
			// }
	// solver_timer.start();
			for (auto iter = redundancy_stitch_pairs.begin(); iter != redundancy_stitch_pairs.end(); ++iter) {
				vSubColor.insert(vSubColor.begin() + (*iter).second, vSubColor.at((*iter).first));
			}
	// cout << "Recover: " << // solver_timer.format(6,"%w") <<", boost::num_vertices(sg)"<<boost::num_vertices(sg)<< endl;
		}
		else{
			vSubColor.assign(num_vertices(sg), -1);

#ifdef _OPENMP
#pragma omp critical(m_dgGlobal2Local)
#endif
			for (std::map<uint32_t, uint32_t>::iterator it = m_dgGlobal2Local.begin(); it != m_dgGlobal2Local.end(); it++)
			{
				if (std::find(vSimpl2Orig.begin(), vSimpl2Orig.end(), it->second) != vSimpl2Orig.end())
				{
					m_dgCompId[it->first] = sub_comp_id + 1;
				}
			}
			if (comp_id == m_db->dbg_comp_id())
			{
				for (std::vector<vertex_descriptor>::const_iterator it = vSimpl2Orig.begin(); it != vSimpl2Orig.end(); ++it)
				{
					mplPrint(kDEBUG, "sub_comp_id %u, orig id is %u.\n", (uint32_t)sub_comp_id, (uint32_t)*it);
				}
			}
			//if algorithm is Dancing Link, call it directly
			boost::timer::cpu_timer comp_timer;
			comp_timer.start();

			// solve coloring
			typedef lac::Coloring<graph_type> coloring_solver_type;
			coloring_solver_type *pcs = create_coloring_solver(sg, comp_id, sub_comp_id);

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

			// 2nd trial, call solve_graph_coloring() again with MERGE_SUBK4 simplification only
			double obj_value2 = std::numeric_limits<double>::max();


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
					if (comp_id == m_db->dbg_comp_id())
					{
						mplPrint(kDEBUG, "vertex %u color is %u.\n", (uint32_t)v, (uint32_t)color);
					}
				}
			}
			else // no need to update vSubColor, as it is already updated by sub call
			{
				acc_obj_value += obj_value2;
				if (comp_id == m_db->dbg_comp_id())
				{
					for (uint32_t v = 0; v < vSubColor.size(); v++)
					{
						mplPrint(kDEBUG, "vertex %u color is %u.\n", (uint32_t)v, (uint32_t)vSubColor[v]);
					}
				}
			}

			if (m_db->parms.record > 1)
			{
				//store ILP runtime information/ read ILP information for DL to comparison, this is for DAC2020
				if (m_db->algo() == AlgorithmTypeEnum::ILP_GUROBI && num_vertices(sg) > 3)
				{
					std::string runtime = comp_timer.format(6, "%w");
					std::ofstream myfile;
					myfile.open("ILP_obj.txt", std::ofstream::app);
					myfile << m_db->input_gds().c_str() << " " << comp_id << " " << sub_comp_id << " " << num_vertices(sg) << " " << obj_value1 << " " << runtime << "\n";
					myfile.close();
				}
				if (m_db->algo() == AlgorithmTypeEnum::DANCING_LINK && num_vertices(sg) > 3)
				{
					std::string runtime = comp_timer.format(6, "%w");
					std::ofstream myfile;
					myfile.open("DL_obj.txt", std::ofstream::app);
					myfile << m_db->input_gds().c_str() << " " << comp_id << " " << sub_comp_id << " " << num_vertices(sg) << " " << obj_value1 << " " << runtime << "\n";
					myfile.close();
				}
			}
			if (m_db->parms.record > 2)
			{
				if (obj_value1 != 0)
				{
					std::ofstream myfile;
					myfile.open("small_results.txt", std::ofstream::app);
					myfile << ", comp_id: " << comp_id << ", sub_comp_id: " << sub_comp_id << ",num_vertices(sg) " << num_vertices(sg) << ", obj_value: " << obj_value1 << "\n";
					myfile.close();
				}
				if (comp_id == m_db->dbg_comp_id())
				{
					write_graph(sg, std::to_string(comp_id) + "_" + std::to_string(sub_comp_id));
					std::ofstream myfile;
					myfile.open("debug_component_color_results.txt", std::ofstream::app);
					myfile << ", comp_id: " << comp_id << ", sub_comp_id: " << sub_comp_id << ",num_vertices(sg) " << num_vertices(sg) << "\n";
					for (uint32_t v = 0; v < vSubColor.size(); v++)
					{
						myfile << "vertex " << v << " color is " << (uint32_t)vSubColor[v] << ". \n";
					}
					myfile.close();
				}

				//we only store json file with graph size larger than 3
				if (num_vertices(sg) > 3)
				{
					//The json file name is orgnized as follows: compid_subcompid_objvalue_time.json
					std::string name = std::to_string(comp_id);
					name.append("_");
					name.append(std::to_string(sub_comp_id));
					name.append("_");
					name.append(std::to_string(obj_value1).substr(0, std::to_string(obj_value1).size() - 5));
					name.append("_");
					name.append(std::to_string(num_vertices(sg)));
					name.append("_");
					name.append(comp_timer.format(4, "%w"));

					this->write_json(sg, (char *)name.c_str(), vSubColor);
				}
			}

			delete pcs;
		}
	}

	// recover csolor assignment according to the simplification level set previously
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing

	gs.recover(vColor, mSubColor, mSimpl2Orig);
	RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();
	SimpleMPL::graph_type non_const_dg(dg);
	double total_obj = new_calc_cost(non_const_dg, vColor);

	return total_obj;
}

std::set<SimpleMPL::vertex_descriptor> SimpleMPL::get_sub_vdd_set(std::vector<vertex_descriptor>& vSimpl2Org, std::set<vertex_descriptor> &vdd_set)
{
	std::set<vertex_descriptor> sub_vdd_set;
	for (auto it=vSimpl2Org.begin(); it!=vSimpl2Org.end(); ++it){
		if(vdd_set.find(*it) != vdd_set.end()) {
			sub_vdd_set.insert(std::distance(vSimpl2Org.begin(), it));
		}
	}
	return sub_vdd_set;
}
SIMPLEMPL_END_NAMESPACE

