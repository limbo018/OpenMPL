/*************************************************************************
    > File Name: SimpleMPL.cpp
    > Author: Yibo Lin, Qi Sun
    > Mail: yibolin@utexas.edu
    > Created Time: Wed May 20 22:38:50 2015
 ************************************************************************/

/*******************************************
in conflict_num() : uncomment color mplAssertmessage
*******************************************/
#include "SimpleMPL.h"
#include "LayoutDBRect.h"
#include "LayoutDBPolygon.h"
#include "RecoverHiddenVertex.h"
#include "DL_MPL.h"
#include "MatrixCover.h"
#include <stack>
#include <time.h>
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
#define CUT_DG
// RECORD is a value to record some information for debug by Wei
#define RECORD 2
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
#if RECORD > 0
	std::ofstream myfile;
	myfile.open ("record.txt", std::ofstream::app);
	myfile << argv[4]<<" " << argv[16]<<"\n";
	myfile.close();
#endif

#if RECORD > 1
	std::ofstream myfile_small;
	myfile_small.open ("small_results.txt", std::ofstream::app);
	myfile_small << argv[4]<<" " << argv[16]<<"\n";
	myfile_small.close();
#endif
	this->read_cmd(argc, argv);
	this->read_gds();
	// graph_type test_graph (2);
	// std::pair<edge_descriptor, bool> e = add_edge(0,1,test_graph);
	// mplAssert(e.second);
	// boost::put(boost::edge_weight, test_graph, e.first, -1);
	// typedef lac::Coloring<graph_type> coloring_solver_type;
	// coloring_solver_type* pcs = create_coloring_solver(test_graph);
	// double obj_value1 = (*pcs)();
	// boost::graph_traits<graph_type>::vertex_iterator vi, vie;
	// for (boost::tie(vi, vie) = vertices(test_graph); vi != vie; ++vi)
	// {
	// 	vertex_descriptor v = *vi;
	// 	int8_t color = pcs->color(v);
	// 	std::cout<<"LIWEI========== color is "<<(int)color<<std::endl;
	// }
	// std::cout<<"LIWEI ======, obj_value1 is "<<obj_value1<<std::endl;
	this->solve();
#ifdef _DGOUT
	return;
#endif
	this->report();
	this->write_gds();
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
        {m_db = new LayoutDBRect;
		is_Rec = true;
		}
    else 
        {m_db = new LayoutDBPolygon;
		is_Rec = false;}
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
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, StitchRelation, m_mAdjVertex);
	}
	else
	{
		GdsWriter writer;
		mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
		writer(m_db->output_gds() + ".gds", *m_db, m_vConflict, m_mAdjVertex, m_db->strname, m_db->unit*1e+6);
	}
}

void SimpleMPL::writeJson(){
	std::ofstream jsonFile;
	char* tmp = const_cast<char *>(m_db->input_gds().c_str());
	char* output_json =(char*)strcat(tmp,".json");
	jsonFile.open(output_json);
	jsonFile<<"[";
	for(uint32_t i = 0; i != m_mAdjVertex.size(); i++){
		if(in_DG[i] == false) 
			{	if(m_db->use_stitch() && StitchRelation[i].size()!=0)
					std::cout<<"BUG FOUND "<<i<<std::endl;
				}
		jsonFile <<"\n";
		jsonFile << "\t{\n";
		jsonFile << "\t\t\"id\": "<<i<<",\n";
		jsonFile << "\t\t\"color\": "<<(int)m_db->vPatternBbox[i]->color()<<",\n";
		jsonFile << "\t\t\"conflict_degree\":"<<m_mAdjVertex[i].size()<<",\n";
		if(m_db->use_stitch())
			jsonFile << "\t\t\"stitch_degree\":"<<StitchRelation[i].size()<<",\n";
		if(m_mAdjVertex[i].size() == 0){
			jsonFile << "\t\t\"conflict\": []";
		}
		else{
			jsonFile << "\t\t\"conflict\": [\n";
			for(uint32_t j = 0; j != m_mAdjVertex[i].size()-1; j++){
				jsonFile << "\t\t\t{\"id\": "<<m_mAdjVertex[i][j]<<"},\n";
			} 
			jsonFile << "\t\t\t{\"id\": " <<m_mAdjVertex[i][m_mAdjVertex[i].size()-1]<<"}\n";
			jsonFile << "\t\t]";
		}
		if(m_db->use_stitch())
			{
				jsonFile<<",\n";
				if(StitchRelation[i].size() == 0){
					jsonFile << "\t\t\"stitch\": []\n";
				}
				else{
					jsonFile << "\t\t\"stitch\": [\n";
					for(uint32_t j = 0; j != StitchRelation[i].size()-1; j++){
						jsonFile << "\t\t\t{\"id\": "<<StitchRelation[i][j]<<"},\n";
					}
					jsonFile << "\t\t\t{\"id\": "<<StitchRelation[i][StitchRelation[i].size()-1]<<"}\n";
					jsonFile << "\t\t]\n";
				}
			}
		else{
			jsonFile<<"\n";
		}

		jsonFile << "\t},";
		//std::cout<<"FIle OUTPUT FINISH"<<i<<std::endl;
	}
	long pos = jsonFile.tellp();
	jsonFile.seekp(pos - 1);
	jsonFile<<"\n]";
	jsonFile.close();
} 


void SimpleMPL::writeJson(graph_type const& sg,std::string graph_count,std::vector<int8_t>& Colors ){
	std::ofstream jsonFile;
	// char* tmp = const_cast<char *>(m_db->input_gds().c_str());
	// char* tmp2 = (char*)strcat(tmp,graph_count);
	 std::cout<<"json file name is"<<graph_count<<std::endl;
	// char* output_json =(char*)strcat(tmp2,".json");
	
	jsonFile.open(m_db->input_gds() + graph_count + ".json");

	SimpleMPL::graph_type tmp_graph = sg;
	jsonFile<<"[";
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	for (boost::tie(vi1, vie1) = boost::vertices(tmp_graph); vi1 != vie1; ++vi1){
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
			if(boost::get(boost::edge_weight, tmp_graph, e12.first) > 0){
				conflict_count ++;
			}
			else{
				stitch_count ++;
			}
			//out << int(v1) <<" "<< int(v2) <<" "<< boost::get(boost::edge_weight, tmp_graph, e12.first)<<"\n";
		}
		if(conflict_count == 0){
			jsonFile << "\t\t\"conflict\": []";
		}
		else{
			jsonFile << "\t\t\"conflict\": [\n";
			int conflict_second_count = 0;
			boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
			for (next2 = vi2; vi2 != vie2; vi2 = next2)
			{
				
				++next2; 
				vertex_descriptor v2 = *vi2;
				std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
				if(boost::get(boost::edge_weight, tmp_graph, e12.first) > 0){
					conflict_second_count ++ ;
					if(conflict_second_count == conflict_count){
						jsonFile << "\t\t\t{\"id\": " <<int(v2)<<"}\n";
					}
					else{
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
				if(stitch_count == 0){
					jsonFile << "\t\t\"stitch\": []";
				}
				else{
					jsonFile << "\t\t\"stitch\": [\n";
					int stitch_second_count = 0;
					boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, tmp_graph);
					for (next2 = vi2; vi2 != vie2; vi2 = next2)
					{
						
						++next2; 
						vertex_descriptor v2 = *vi2;
						std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, tmp_graph);
						if(boost::get(boost::edge_weight, tmp_graph, e12.first) < 0){
							stitch_second_count ++ ;
							if(stitch_second_count == stitch_count){
								jsonFile << "\t\t\t{\"id\": " <<int(v2)<<"}\n";
							}
							else{
								jsonFile << "\t\t\t{\"id\": "<<int(v2)<<"},\n";
							}
							
						}
						//out << int(v1) <<" "<< int(v2) <<" "<< boost::get(boost::edge_weight, tmp_graph, e12.first)<<"\n";
					}

					jsonFile << "\t\t]";
				}
				
			}
		else{
			jsonFile<<"\n";
		}

		jsonFile << "\t},";
		//std::cout<<"FIle OUTPUT FINISH"<<i<<std::endl;
	}
	long pos = jsonFile.tellp();
	jsonFile.seekp(pos - 1);
	jsonFile<<"\n]";
	jsonFile.close();
} 

void SimpleMPL::writeGraph(graph_type const& sg,std::string const filename, double& cost){
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
void SimpleMPL::outStat(){
	int TCE = 0;
	int TSE = 0;
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"Now node size is "<<m_mAdjVertex.size()<<", dg size is "<<in_DG.size()<<std::endl;
	#endif
    for (uint32_t i = 0; i != m_mAdjVertex.size(); i++)
    {
		if(in_DG[i] == false) continue;
		#ifdef DEBUG_LIWEI_DONE
			std::cout<<"LIWEI: for DG node "<<i<<", there are "<<m_mAdjVertex[i].size()<<std::endl;
		#endif
		for(uint32_t j = 0; j != m_mAdjVertex[i].size(); j++)
			{if(in_DG[m_mAdjVertex[i][j]] == false) continue;
			TCE ++;}
	}
    for (uint32_t i = 0; i != StitchRelation.size(); i++)
    {
		TSE += StitchRelation[i].size();
	}
	printf("\n");
	printf("=================== Graph Simplification Information ===================\n");
	printf("# of Total Conflict Edges                              : %d\n", TCE/2);
	printf("# of Total Stitch Edges                                : %d\n", TSE/2);
	printf("# of Total Wires                                       : %d\n", (int)m_db->vPatternBbox.size());
	printf("# of DG component after  DG Bridge Division            : %d\n", (int)DG_num); 
	printf("# of DGs after  DG Bridge Division                     : %d\n", (int)std::count(in_DG.begin(),in_DG.end(),true)); 
	printf("========================================================================\n");



}

void SimpleMPL::solve()
{
    // skip if no uncolored layer 
    if (m_db->parms.sUncolorLayer.empty())
        return;
	boost::timer::cpu_timer total_timer;
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
	// TODO: if the vBookmark generation is time-consuming, some optimization algorithm may be proposed (set is better than vector)
	std::vector<uint32_t>().swap(vBookmark);
	vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}
	this->calBoundaries();
	this->setVddGnd(); //perhapes we should also consider pshape->getPointNum()==4
	
	clock_t begin = clock();
	this->lgSimplification();

	
	if (m_db->use_stitch()) //returns whether use stitch
	{

		//GdsWriter writer;
		//writer.write_Simplification(m_db->output_gds() + "_lgSimplification.gds", *m_db, m_vCompId, m_mAdjVertex, in_DG, isVDDGND, true);
		#ifdef DEBUG_LIWEI_DONE
			std::cout<<"LIWEI: projection starts"<<std::endl;
		#endif
		this->projection();		///< vBookmark has already been updated in projection()
		
		clock_t end = clock();
		mplPrint(kINFO, "Projection takes  %f.\n", (double)(end - begin) / CLOCKS_PER_SEC);
	}
		std::vector<uint32_t>().swap(dgCompId);
		dgCompId.assign(m_db->vPatternBbox.size(), 0);
		globalCompId = 1;
		std::cout << "======= Nodes in DG : " << std::count(in_DG.begin(), in_DG.end(), true) << std::endl;
		
		// update VddGnd information
		this->setVddGnd();

		
		vdd_multi_comp.resize(m_db->vPatternBbox.size());
		//this->dgSimplColoring(); 
		if (m_db->use_stitch())
		{
			GdsWriter writer;
			writer.write_Simplification(m_db->output_gds() + "_dgSimplification.gds", *m_db, dgCompId, StitchRelation, std::vector<bool>(), isVDDGND, false);
		}
		std::cout << "This graph finally been decomposed into : " << globalCompId << std::endl;

		for (uint32_t i = 0; i < vdd_multi_comp.size(); i++)
		{
			if (isVDDGND[i])
			{
				std::cout << "VDD " << i << " in : ";
				for (std::vector<uint32_t>::iterator it = vdd_multi_comp[i].begin(); it != vdd_multi_comp[i].end(); it++)
					std::cout << *it << " ";
				std::cout << std::endl;
			}
		}
#ifdef _DGOUT
		return;
#endif

	boost::timer::cpu_timer t;

	// thread number controled by user option 
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif 
    for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
    {
        // construct a component 
        std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin()+vBookmark[comp_id];
        std::vector<uint32_t>::iterator itEnd = (comp_id+1 != m_comp_cnt)? m_vVertexOrder.begin()+vBookmark[comp_id+1] : m_vVertexOrder.end();
        // solve component 
        // pass iterators to save memory 
        this->solve_component(itBgn, itEnd, comp_id);
	}
#if RECORD > 0
	if(m_db->use_stitch()){
		std::ofstream myfile,color_time,total_time;
		myfile.open ("record.txt", std::ofstream::app);
		color_time.open("color_w_stitch.txt",std::ofstream::app);
		total_time.open("total_w_stitch.txt",std::ofstream::app);
		myfile << t.format(2, "color time %ts(%p%), %ws real")<<"\n";
		myfile << total_timer.format(2, "total time %ts(%p%), %ws real")<<"\n";
		if(m_db->algo() == AlgorithmTypeEnum::ILP_GURBOI){
			color_time<<"\\\\\n"<<m_db->input_gds();
			total_time<<"\\\\\n"<<m_db->input_gds();
		}
		color_time<<t.format(2, "& %t  &  %w");
		total_time<<total_timer.format(2, "& %t  &  %w");
		myfile.close();
		color_time.close();
		total_time.close();
	}
	else{
		std::ofstream myfile,color_time,total_time;
		myfile.open ("record.txt", std::ofstream::app);
		color_time.open("color_wo_stitch.txt",std::ofstream::app);
		total_time.open("total_wo_stitch.txt",std::ofstream::app);
		myfile << t.format(2, "color time %ts(%p%), %ws real")<<"\n";
		myfile << total_timer.format(2, "total time %ts(%p%), %ws real")<<"\n";
		if(m_db->algo() == AlgorithmTypeEnum::ILP_GURBOI){
			color_time<<"\\\\\n"<<m_db->input_gds();
			total_time<<"\n"<<m_db->input_gds();
		}
		color_time<<t.format(2, "& %t  &  %w");
		total_time<<total_timer.format(2, "& %t  &  %w");
		myfile.close();
		color_time.close();
		total_time.close();
	}

#endif
	this->outStat();
	//this->writeJson();
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

void SimpleMPL::lgSimplification()
{
	std::vector<bool>().swap(in_DG);
	std::vector<bool>().swap(articulation_vec);
	in_DG.assign(m_db->vPatternBbox.size(), false);
	articulation_vec.assign(m_db->vPatternBbox.size(), false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_db->thread_num())
#endif
	for (uint32_t comp_id = 0; comp_id < m_comp_cnt; ++comp_id)
	{
		// construct a component 
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		// projection on each component
		this->lgSimplification(itBgn, itEnd, comp_id);
	}
}

void SimpleMPL::lgSimplification(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id)
{
#ifdef DEBUG_LIWEI_DONE
	std::cout << "Now in component " << comp_id << " projection.\n";
#endif
	(void)comp_id;
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	uint32_t simplify_strategy = graph_simplification_type::HIDE_SMALL_DEGREE;
#if LIWEI_BEFORE
	simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;
#else
	//simplify_strategy |= graph_simplification_type::BICONNECTED_COMPONENT;
#endif
	

	graph_simplification_type gs(dg, m_db->color_num());
	
	gs.set_isVDDGND(vdd_set);
	gs.simplify(simplify_strategy);
	std::vector<bool> projLocal(mGlobal2Local.size(), false);

	std::vector<vertex_descriptor> all_articulations;

	gs.get_articulations(all_articulations);
#ifdef DEBUG_LIWEI_DONE
	std::cout<<"LIWEI: articulations size is "<<all_articulations.size()<<", component number is "<<gs.num_component()<<std::endl;
#endif
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
			in_DG[it->first] = true;
		if (temp_set.find(it->second) != temp_set.end())
			articulation_vec[it->first] = true;
	}
}

void SimpleMPL::dgSimplColoring()
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
		std::vector<uint32_t>::iterator itBgn = m_vVertexOrder.begin() + vBookmark[comp_id];
		std::vector<uint32_t>::iterator itEnd = (comp_id + 1 != m_comp_cnt) ? m_vVertexOrder.begin() + vBookmark[comp_id + 1] : m_vVertexOrder.end();
		std::vector<std::vector<vertex_descriptor> > comp_vertex;
		std::map<vertex_descriptor, std::set<uint32_t> >& arti_point = m_mArtiPoints[comp_id];
		// Simplification on each component

		this->liweidgSimplColoring(itBgn, itEnd, comp_id, arti_point, comp_vertex);
	}
#ifdef _DGOUT
	return;
#endif 
}


/*
for (std::vector<vertex_descriptor>::iterator it = sub_vSimpl2Orig.begin(); it != sub_vSimpl2Orig.end(); it++)
{
	// std::cout << "*it : " << *it << std::endl;
	// std::cout << "vSimpl2Orig.at(*it) : " << vSimpl2Orig.at(*it) << std::endl;
	vertex_descriptor curpid = mLocal2Global.at(vSimpl2Orig.at(*it));
	// std::cout << curpid << std::endl;
				
	dgCompId[curpid] = globalCompId;
	if (isVDDGND[curpid])
	{
		vdd_multi_comp[curpid].push_back(globalCompId);
	}
}
globalCompId++;
*/

void SimpleMPL::dgSimplColoring(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::map<vertex_descriptor, std::set<uint32_t> >& m_ArtiPoint, std::vector<std::vector<vertex_descriptor> >& m_CompVertex)
{
	(void)m_ArtiPoint;
	(void)m_CompVertex;
	std::cout << "\n\n\n**************************************************************\n";
	std::cout << "Debug| dgSimplification comp " << comp_id << std::endl;
	uint32_t pattern_cnt = itEnd - itBgn;
	graph_type dg(pattern_cnt);
	std::vector<int8_t> vColor(pattern_cnt, -1);
	std::map<uint32_t, uint32_t> mGlobal2Local;
	std::set<vertex_descriptor> vdd_set;

	typedef lac::GraphSimplification<graph_type>   graph_simplification_type;
	this->construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, false);
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI: DG construct complete"<<std::endl;
	#endif
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
#ifdef DEBUG_LIWEI_DONE
		std::cout << "\n==================================================\n";
		std::cout << "Debug| dgSimplification sub_comp " << small_comp_id << std::endl;
#endif
		graph_type small_g;
		
		std::vector<int8_t>& small_comp_color = big_comp_color[small_comp_id];	///< used to store coloring results of current small component
		(void) small_comp_color;
		std::vector<vertex_descriptor>& vSimplSmall2OriBig = mSmall2Big[small_comp_id];		///< used to store the simplificaiton mapping information from current small component to original big component

		gs.simplified_graph_component(small_comp_id, small_g, vSimplSmall2OriBig);

		uint32_t temp = 0;
		std::set<vertex_descriptor> small_vdd_set;
		std::vector<int8_t> small_vColor(num_vertices(small_g), -1);
		///< construct small vColor
		///< also get vdd_set for current sub component
		for (std::vector<vertex_descriptor>::iterator it = vSimplSmall2OriBig.begin(); it != vSimplSmall2OriBig.end(); it++, temp++)
		{
			if (gs.whether_VDDGND(*it))
				small_vdd_set.insert(temp);
			small_vColor[temp] = vColor[*it];
		}

		///< ***************************************************///
		///< Now simplify current small component
		///< **************************************************///

		graph_simplification_type simplified_small_g(small_g, m_db->color_num());	///< construct small graph simplification object
		///< merge VDD nodes in sg, and conduct biconnected devision
		simplified_small_g.precolor(small_vColor.begin(), small_vColor.end());
		simplified_small_g.mergeVDD(small_vdd_set);
		simplified_small_g.simplify(graph_simplification_type::BICONNECTED_COMPONENT);
		std::vector<std::vector<int8_t> > small_comp_color_2 (simplified_small_g.num_component());	///< used to store the coloring results of sub small components, and recover small_comp_color;
		std::vector<std::vector<vertex_descriptor> > sub_small2small(simplified_small_g.num_component());	///< store the mapping information from sub small components to current small component
		
		double acc_obj_value = 0;

		for (uint32_t ss_comp_id = 0; ss_comp_id < simplified_small_g.num_component(); ++ss_comp_id)
		{
#ifdef DEBUG_LIWEI_DONE
			std::cout << "\n===============================\n";
			std::cout << "Debug| dgSimplification ss_comp " << ss_comp_id << std::endl;
			std::cout << "Debug| dgSimplification with nodes : ";
#endif
			graph_type sub_small_simpl_sg;
			
			std::vector<vertex_descriptor>& sub_vSimpl2OrigSmall = sub_small2small[ss_comp_id];
			simplified_small_g.simplified_graph_component(ss_comp_id, sub_small_simpl_sg, sub_vSimpl2OrigSmall);
			
			std::vector<int8_t>& sub_small_comp_color = small_comp_color_2[ss_comp_id];	///< used to store 
			sub_small_comp_color.assign(num_vertices(sub_small_simpl_sg), -1);

			typedef lac::Coloring<graph_type> coloring_solver_type;
			coloring_solver_type* pcs = create_coloring_solver(sub_small_simpl_sg);

			boost::graph_traits<graph_type>::vertex_iterator vi, vie;
			for (boost::tie(vi, vie) = vertices(sub_small_simpl_sg); vi != vie; vi++)
			{
				vertex_descriptor v = *vi;
				int8_t color = small_vColor[sub_vSimpl2OrigSmall[v]];
				if (color >= 0 && color < m_db->color_num())
				{
					pcs->precolor(v, color);
					sub_small_comp_color[v] = color; ///< necessary for merge_K4 coloring
				}
			}
			double obj_value1 = (*pcs)();
			double obj_value2 = std::numeric_limits<double>::max();
#ifdef DEBUG_NONINTEGERS
			if (obj_value1 >= 1 && boost::num_vertices(sub_small_simpl_sg) > 4 && (m_db->algo() == AlgorithmTypeEnum::SDP_CSDP) && (simplify_strategy & graph_simplification_type::MERGE_SUBK4) == 0) // merge K4 is not performed
				obj_value2 = merge_K4_coloring(sub_small_simpl_sg, sub_small_comp_color);
#endif
			if (obj_value1 < obj_value2)
			{
				acc_obj_value += obj_value1;
				///< collect  coloring results from simplified graph
				for (boost::tie(vi, vie) = vertices(sub_small_simpl_sg); vi != vie; vi++)
				{
					vertex_descriptor v = *vi;
					int8_t color = pcs->color(v);
					mplAssert(color >= 0 && color < m_db->color_num());
					sub_small_comp_color[v] = color;
				}
			}
			else // no need to update sub_small_comp_color, since it has already been updated by merge_K4_coloring
				acc_obj_value += obj_value2;
		}
	}
}


void SimpleMPL::liweidgSimplColoring(std::vector<uint32_t>::const_iterator itBgn, std::vector<uint32_t>::const_iterator itEnd, uint32_t comp_id, std::map<vertex_descriptor, std::set<uint32_t> >& m_ArtiPoint, std::vector<std::vector<vertex_descriptor> >& m_CompVertex)
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
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI: DG construct complete"<<std::endl;
	#endif
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
	#ifdef DEBUG_LIWEI_DONE
		if(pattern_cnt > vHiddenVertices.size())
			std::cout << "Debug| dgSimplification comp " << comp_id << ", with left vertex number "<< pattern_cnt-vHiddenVertices.size()<<", subcomponent number "<<gs.num_component()<<std::endl;
	#endif
	for (uint32_t small_comp_id = 0; small_comp_id < gs.num_component(); ++small_comp_id)
	{

		graph_type small_g;
		
		std::vector<int8_t>& small_comp_color = big_comp_color[small_comp_id];	///< used to store coloring results of current small component
		(void)small_comp_color;
		std::vector<vertex_descriptor>& vSimplSmall2OriBig = mSmall2Big[small_comp_id];		///< used to store the simplificaiton mapping information from current small component to original big component

		gs.simplified_graph_component(small_comp_id, small_g, vSimplSmall2OriBig);

		std::set<vertex_descriptor> small_vdd_set;
		std::vector<int8_t> small_vColor(num_vertices(small_g), -1);
#ifdef DEBUG_LIWEI_DONE
		std::cout << "Debug| dgSimplification sub_comp " << small_comp_id << " with vertices "<<num_vertices(small_g)<< std::endl;
#endif
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
		coloring_solver_type* pcs = create_coloring_solver(sg);

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
	std::vector<bool>().swap(isVDDGND);
	int threshold = 6000;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	isVDDGND.assign(vertex_num, false);
	for (uint32_t i = 0; i < vertex_num; i++)
	{
#if LIWEI_BEFORE
		if (boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) >= threshold)
		{
			isVDDGND[i] = true;
			// m_db->vPatternBbox[i]->color(m_db->color_num() - 1);
		}
#else
		if (isLongEnough(m_db->vPatternBbox[i]))
		{
			isVDDGND[i] = true;
			// m_db->vPatternBbox[i]->color(m_db->color_num() - 1);
		}
		if(boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL) == threshold ){
			isVDDGND[i] = false;
			std::cout<<"dirty coding for oracle cases"<<std::endl;
		}
#endif
#if 0
		else
		{
			coordinate_difference width = boost::polygon::delta(*m_db->vPatternBbox[i], gtl::HORIZONTAL);
			coordinate_difference height = boost::polygon::delta(*m_db->vPatternBbox[i], gtl::VERTICAL);
		}
#endif
	}
	return;
}

void SimpleMPL::calBoundaries(){
	std::vector<coordinate_type>().swap(boundaries);
	mplAssert(boundaries.empty());
	uint32_t vertex_num = m_db->vPatternBbox.size();
	coordinate_type tmp_bound[4] = {INT_MAX,0,INT_MAX ,0};
	for (uint32_t i = 0; i < vertex_num; i++){
		coordinate_type test0 = gtl::xl(*(m_db->vPatternBbox[i]));
		coordinate_type test1 = gtl::xh(*(m_db->vPatternBbox[i]));
		coordinate_type test2 = gtl::yl(*(m_db->vPatternBbox[i]));
		coordinate_type test3 = gtl::yh(*(m_db->vPatternBbox[i]));
		tmp_bound[0] = std::min(test0,tmp_bound[0]);
		tmp_bound[1] = std::max(test1,tmp_bound[1]);
		tmp_bound[2] = std::min(test2,tmp_bound[2]);
		tmp_bound[3] = std::max(test3,tmp_bound[3]);
	}
	boundaries.push_back(tmp_bound[0]);
	boundaries.push_back(tmp_bound[1]);
	boundaries.push_back(tmp_bound[2]);
	boundaries.push_back(tmp_bound[3]);
#if LIWEI_DEBUG
	std::cout<<"boundries of layout is calculated and the value is "<<boundaries[0]<<','\
																	<<boundaries[1]<<','\
																	<<boundaries[2]<<','\
																	<<boundaries[3]<<','\
	<<std::endl;
#endif
	return;
}

void SimpleMPL::projection()
{
	uint32_t vertex_num = m_db->vPatternBbox.size();
	std::vector<rectangle_pointer_type> rect_vec = m_db->polyrect_patterns();	///< original rectangle list
	std::vector<uint32_t> Poly_Rect_begin;
    if(is_Rec)
        {
			std::cout<<"LIWEI: RECTANGLE INPUT"<<std::endl;
			for (uint32_t i = 0; i < vertex_num; i++)
				Poly_Rect_begin.push_back(i);
		}
    else
	{
		Poly_Rect_begin = m_db->PolyRectBgnLoc();				///< original polygons mapping to rectangles
		std::cout<<"LIWEI: POLYGON INPUT"<<std::endl;
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
	std::vector<std::vector<uint32_t> >().swap(ori2new_polygon);	///< map from original polygons to new polygons
	std::vector<uint32_t>().swap(new2ori_polygon);					///< map from new polygons to old polygons
	ori2new_polygon.resize(vertex_num);
	
	uint32_t new_polygon_id = -1;						///< new polygon id
	uint32_t new_rectangle_id = -1;						///< new rectangle id
	uint32_t stitch_edge_number = 0;    
	//LIWEI: projection data preparetion finished"<<std::endl;
       
	for (uint32_t i = 0; i < vertex_num; i++)
	{
		uint32_t v = m_vVertexOrder[i];
		uint32_t comp_id = m_vCompId[v];
		rectangle_pointer_type const & pPattern = m_db->vPatternBbox[v];
		uint32_t pid = pPattern->pattern_id();
		uint32_t start_idx = Poly_Rect_begin[pid];
		uint32_t end_idx = Poly_Rect_end[pid];
		
		std::vector<uint32_t>& nei_vec = m_mAdjVertex[pid];
		//if (in_DG[pid] && isVDDGND[pid] == false) 
		if (in_DG[pid] && articulation_vec[pid] == false && isVDDGND[pid] == false) 
		{ 
			std::vector<rectangle_pointer_type> poss_nei_vec;
			for (std::vector<uint32_t>::iterator it = nei_vec.begin(); it != nei_vec.end(); it++)
			{

				if (in_DG[*it] == false) continue;
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
				splitRectangle(rect, split, poss_nei_vec);
				for (uint32_t i = 0; i < split.size(); i++)
					rect_split.push_back(std::make_pair(split[i], rect.pattern_id()));
				stitch_edge_number += split.size() - 1;
#if DEBUG_LIWEI_DONE
				if(split.size()>1){
					std::cout<<"rect "<<rect.pattern_id()<<" add stitch :";
					for (uint32_t i = 0; i < split.size(); i++)
					{	
						rectangle_pointer_type pRect = split[i];
						std::cout<<gtl::xl(*pRect)<<" ";
					}
						
					std::cout<<std::endl;
				}
				
#endif
			} ///< for

			uint32_t pivot = new_polygon_id;
			std::vector<uint32_t> new_polygon_id_list;
			// reconstruct polygons, also generate stitch relationships
			liweireconstruct_polygon(new_polygon_id, new_polygon_id_list, rect_split);
			mplAssert(new_polygon_id_list.size() == rect_split.size());

			for (uint32_t i = 0; i < rect_split.size(); i++)
			{
				rect_split[i].first->pattern_id(++new_rectangle_id);
				new_rect_vec.push_back(rect_split[i].first);
				new_Rect2ParentPoly.push_back(new_polygon_id_list[i]);
				
				// insert new polygon
				if (pivot != new_polygon_id_list[i])
				{
					pivot = new_polygon_id_list[i];
					new2ori_polygon.push_back(pid);
					ori2new_polygon[pid].push_back(pivot);
					new_vertex_order.push_back(pivot);
					new_vCompId_vec.push_back(comp_id);
				}
			} ///< for rect_list
		} ///< if, stitch generation core
		else
		{
			//LIWEI: NON-AVALIBLE vertex are being handled "<<v<<std::endl;

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
			new2ori_polygon.push_back(pid);
			new_vertex_order.push_back(new_polygon_id);
			new_vCompId_vec.push_back(comp_id);
			ori2new_polygon[pid].push_back(new_polygon_id);
			StitchRelation.push_back(std::vector<uint32_t>());
		} ///< else (in_DG)
	} ///< for all vertices

#if LIWEI_DEBUG
	int stitchCount = 0;
	mplAssert(StitchRelation.size() == new2ori_polygon.size());
	for(int i= 0;i<StitchRelation.size();i++){
		std::vector<uint32_t> tmp = StitchRelation[i];
		stitchCount += tmp.size();
	}
	std::cout<<"LIWEI: stitch_edge_number is "<<stitch_edge_number<<", StitchRelation size is "<< stitchCount<<std::endl;
#endif
	m_db->refresh(new_rect_vec, new_Rect2ParentPoly);

#ifdef DEBUG_LIWEI_DONE
	std::cout<<"LIWEI: refresh done"<<std::endl;
#endif
	updateConflictRelation(); 


	std::vector<uint32_t>().swap(m_vCompId);
	m_vCompId.swap(new_vCompId_vec);

	std::vector<uint32_t>().swap(m_vVertexOrder);
	m_vVertexOrder.swap(new_vertex_order);

	std::vector<uint32_t>().swap(vBookmark);
	vBookmark.resize(m_comp_cnt);
	for (uint32_t i = 0; i != m_vVertexOrder.size(); ++i)
	{
		if (i == 0 || m_vCompId[m_vVertexOrder[i - 1]] != m_vCompId[m_vVertexOrder[i]])
			vBookmark[m_vCompId[m_vVertexOrder[i]]] = i;
	}

	return;
}

void SimpleMPL::printVector(std::vector<uint32_t>& vector){
	for(uint32_t i = 0; i<vector.size();i++){
		std::cout<<vector[i]<<" ";
	}
}
void SimpleMPL::updateConflictRelation()
{
	// now update adjacency list, we still need original adjacency list.
	// here we use set to remove duplicated neighboring polygons
	std::vector<std::set<uint32_t> >new_mAdjVertex;
	uint32_t vertex_num = m_db->vPatternBbox.size();
	new_mAdjVertex.resize(m_db->vPatternBbox.size());
	uint32_t edge_num = 0;
	uint32_t dg_conflict = 0;
	uint32_t dg_illegal_conflict = 0;
	/*
	// generate new index information
	std::vector<uint32_t> new_poly_rect_begin = m_db->PolyRectBgnLoc();
	std::vector<uint32_t> new_poly_rect_end;
	
	mplAssert(vertex_num == m_db->vPatternBbox.size());

	new_poly_rect_end.resize(vertex_num);
	for (uint32_t i = 0; i < vertex_num - 1; i++)
		new_poly_rect_end[i] = new_poly_rect_begin[i + 1] - 1;
	new_poly_rect_end[vertex_num - 1] = new_rect_vec.size() - 1;
	*/

	//update in_DG
	std::vector<bool> tmp_in_DG;
	tmp_in_DG.assign(m_db->vPatternBbox.size(), false);
	for(uint32_t i = 0; i< tmp_in_DG.size();i++){
		if(in_DG[new2ori_polygon[i]])
			tmp_in_DG[i] = true;
	}
	in_DG.swap(tmp_in_DG);

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
			std::vector<uint32_t>::iterator itr = find(StitchRelation[pPattern->pattern_id()].begin(), StitchRelation[pPattern->pattern_id()].end(), pAdjPattern->pattern_id());
			//if (itr != StitchRelation[pPattern->pattern_id()].end())
			if(new2ori_polygon[pPattern->pattern_id()] == new2ori_polygon[pAdjPattern->pattern_id()])
			{
				#ifdef DEBUG_LIWEI_DONE
				if(in_DG[pPattern->pattern_id()]){
					dg_illegal_conflict ++;
					std::cout<<"LIWEI: in the same origin pattern"<<std::endl;
					std::cout<<gtl::xl(*(pPattern))<<" "<<gtl::yl(*(pPattern))<<std::endl;
					std::cout<<gtl::xl(*(pAdjPattern))<<" "<<gtl::yl(*(pAdjPattern))<<std::endl;
				}
				#endif
				continue;
			}
				
			mplAssert(itr == StitchRelation[pPattern->pattern_id()].end());
			mplAssert(pAdjPattern != pPattern);
			if (pAdjPattern != pPattern) // skip pPattern itself 
			{
				mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
				// we consider euclidean distance
				// use layoutdb_type::euclidean_distance to enable compatibility of both rectangles and polygons
				coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
				if (distance < m_db->coloring_distance)
				{
					if(in_DG[pPattern->pattern_id()] && in_DG[pAdjPattern->pattern_id()])	dg_conflict++;
					vAdjVertex.insert(pAdjPattern->pattern_id());
				}
			}
		}
		vAdjVertex.swap(vAdjVertex); // shrink to fit, save memory 
		edge_num += vAdjVertex.size();
	}
	std::cout << "After stitch insertion, it has " << dg_conflict << " dg_conflict. and dg_illegal_conflict\n"<<dg_illegal_conflict;

	///< update m_mAdjVertex
	std::vector<std::vector<uint32_t> >().swap(m_mAdjVertex);
	m_mAdjVertex.resize(new_mAdjVertex.size());
	for (uint32_t i = 0; i < new_mAdjVertex.size(); i++)
		m_mAdjVertex[i].insert(m_mAdjVertex[i].end(), new_mAdjVertex[i].begin(), new_mAdjVertex[i].end());

	/*
	// traverse all the original polygons to get the original neighbor list
	for (uint32_t i = 0, ie = m_mAdjVertex.size(); i < ie; i++)
	{
		// list all the possible polygon neighbors
		std::vector<uint32_t> poss_nei;
		bool split_flag = false;
		for (uint32_t j = 0, je = m_mAdjVertex[i].size(); j < je; j++)
		{
			uint32_t nei_id = m_mAdjVertex[i][j];
			for (uint32_t it = 0; it < ori2new_polygon[nei_id].size(); it++)
				poss_nei.push_back(ori2new_polygon[nei_id][it]);
		}
		// traverse all the newly-generated nodes in polygon i
		for (std::vector<uint32_t>::iterator it = ori2new_polygon[i].begin(); it != ori2new_polygon[i].end(); it++)
		{
			uint32_t start_idx = new_poly_rect_begin[*it];
			uint32_t end_idx = new_poly_rect_end[*it];
			for (std::vector<uint32_t>::iterator nei_poly = poss_nei.begin(); nei_poly != poss_nei.end(); nei_poly++)
			{
				if (*nei_poly > *it)
					continue;
				for (uint32_t now_rect = start_idx; now_rect <= end_idx; now_rect++)
				{
					uint32_t nei_start_idx = new_poly_rect_begin[*nei_poly];
					uint32_t nei_end_idx = new_poly_rect_end[*nei_poly];
					// traverse all rectangles in current newly-generated neighbor polygon
					for (uint32_t nei_rect = nei_start_idx; nei_rect <= nei_end_idx; nei_rect++)
					{
						uint32_t nei_rect_pid = new_rect_vec[nei_rect]->pattern_id();
						coordinate_difference distance = boost::geometry::distance(*new_rect_vec[now_rect], *new_rect_vec[nei_rect_pid]);
						if (distance < m_db->coloring_distance)
						{
							new_mAdjVertex[*it].insert(*nei_poly);
							new_mAdjVertex[*nei_poly].insert(*it);
							break;
						}// if
					} // for
				} // for
			} // for
		} // for
		
	} // for
*/
}
void SimpleMPL::findConflictRects(rectangle_pointer_type& prec, std::vector<rectangle_pointer_type>& conflict_rects){
	rectangle_type rect(*prec);
	gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
	gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
	for (rtree_type::const_query_iterator itq = m_db->tPatternBbox.qbegin(bgi::intersects(rect));
			itq != m_db->tPatternBbox.qend(); ++itq)
	{
		rectangle_pointer_type const& pAdjPattern = *itq;
		if(gtl::xl(*(pAdjPattern))<=gtl::xl(*(prec)) && gtl::xh(*(pAdjPattern))>=gtl::xh(*(prec)) \
			&& gtl::yl(*(pAdjPattern))<=gtl::yl(*(prec)) && gtl::yh(*(pAdjPattern))>=gtl::yh(*(prec)))
			continue;
		std::vector<rectangle_pointer_type> vPolyRectPattern = m_db->polyrect_patterns();
		std::vector<uint32_t> vPolyRectBeginId = m_db->PolyRectBgnLoc();
		uint32_t num_polyrects = vPolyRectPattern.size();
		uint32_t adj_id =  pAdjPattern->pattern_id();
		coordinate_difference distance = std::numeric_limits<coordinate_difference>::max();
   	 	uint32_t polyRectId1e = (adj_id+1 == vPolyRectBeginId.size())? num_polyrects : vPolyRectBeginId[adj_id+1];
		for (uint32_t polyRectId1 = vPolyRectBeginId[adj_id]; polyRectId1 != polyRectId1e; ++polyRectId1)
			distance = std::min(distance, (coordinate_difference)gtl::euclidean_distance(*(vPolyRectPattern[polyRectId1]), *prec) );
		if(distance < m_db->coloring_distance)
			conflict_rects.push_back(pAdjPattern);
	}
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI: find confilict patterning boxes with size "<<conflict_rects.size()<<std::endl;
	#endif

}

bool SimpleMPL::isLongEnough(rectangle_pointer_type& rec){
		if (boost::polygon::delta(*rec, gtl::HORIZONTAL) >= 0.64*(boundaries[1]-boundaries[0])||\
		boost::polygon::delta(*rec, gtl::VERTICAL) >= 0.64*(boundaries[3]-boundaries[2]))
			return true;
		else
			return false;
}

//TODO: FINISH this function
bool SimpleMPL::bIsBoundSlice(rectangle_pointer_type& prec1,rectangle_pointer_type& prec2){
	(void) prec1;
	(void) prec2;
	return false;
}
bool SimpleMPL::canIntroStitch(rectangle_pointer_type& prec1,rectangle_pointer_type& prec2){
	if(m_db->color_num() <= 1) return false;
	if(m_db->color_num() == 2)
	{
		std::cout<<"THE PROGRAM ONLY SUPPORTS triple coloring now"<<std::endl;
		return false;
	}
	mplAssert(m_db->color_num() == 3);
	std::vector<rectangle_pointer_type> vID1;
	std::vector<rectangle_pointer_type> vID2;
	findConflictRects(prec1,vID1);
	findConflictRects(prec2,vID2);
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI : find conflicts size "<<vID1.size()<<" "<<vID2.size()<<std::endl;
	#endif
	if(vID1.size() == 0 || vID2.size() == 0) return false;
	//TODO: use && and change recovery stratage
	//TODO: boundslice and share diff dirt neighbor 
	//remove the stitches on vdd
  	if (isLongEnough(prec1) || isLongEnough(prec2)) return false;
	//remove boundslice
	coordinate_difference threshold = m_db->coloring_distance / (m_db->color_num() + 2 );
	coordinate_difference hard_code = -3850;
	coordinate_difference hard_code_y = -1730;
	if(gtl::xl(*(prec1)) == hard_code && gtl::yl(*(prec1)) == hard_code_y){
		for (uint32_t i=0; i<vID1.size(); i++)
		{
			rectangle_pointer_type id = vID1[i];
			std::cout<<gtl::xl(*(id))<<std::endl;
			std::cout<<gtl::yl(*(id))<<std::endl;
			std::cout<<gtl::xh(*(id))<<std::endl;
			std::cout<<gtl::yh(*(id))<<std::endl;
		}
	}
	if(gtl::xl(*(prec2)) == hard_code && gtl::yl(*(prec2)) == hard_code_y){
		for (uint32_t i=0; i<vID2.size(); i++)
		{
			rectangle_pointer_type id = vID2[i];
			std::cout<<gtl::xl(*(id))<<std::endl;
			std::cout<<gtl::yl(*(id))<<std::endl;
			std::cout<<gtl::xh(*(id))<<std::endl;
			std::cout<<gtl::yh(*(id))<<std::endl;
		}
	}
	if((gtl::xh(*prec1) - gtl::xl(*prec1)) < threshold || (gtl::yh(*prec1) - gtl::yl(*prec1))< threshold)
	{
		if (true == bIsBoundSlice(prec1, prec2)) return false;
	}
	if(gtl::xh(*prec2) - gtl::xl(*prec2) < threshold || gtl::yh(*prec2) - gtl::yl(*prec2) < threshold)
	{
		if (true == bIsBoundSlice(prec2, prec1)) return false;
	}
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
	if (bFind1 || bFind2) return true;
	//LIWEI: bFind1 || bFind2 changes to bFind1 && bFind2, since I think each box should have their own special adj boxes
	//results show that this change is too strict that conflict number increases while 
	return false;	

}
void SimpleMPL::liweireconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t>& new_polygon_id_list, std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list)
{
	Graph G(rect_list.size());
	int stitchCount = 0;
	uint32_t start = polygon_id + 1;	
	std::vector<uint32_t>().swap(new_polygon_id_list);
	new_polygon_id_list.assign(rect_list.size(), std::numeric_limits<uint32_t>::max());
	std::vector<std::vector<uint32_t> > stitch_list;
	std::vector<std::vector<uint32_t> > rec_stitch_list;
	std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_temp;
	std::vector<uint32_t> new_polygon_id_temp;
	rec_stitch_list.resize(rect_list.size());
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI: pid "<<polygon_id<<", rect size "<<rect_list.size()<<std::endl;
	#endif
	for(uint32_t i=0;i<rect_list.size();i++){
		rectangle_pointer_type prec1 = rect_list[i].first;
		for(uint32_t j=i+1;j<rect_list.size();j++){
			rectangle_pointer_type prec2 = rect_list[j].first;
			coordinate_difference distance = boost::geometry::distance(*prec1, *prec2);
			#ifdef DEBUG_LIWEI_DONE
				std::cout<<"LIWEI: distance between "<<i<<" and "<<j<<" is "<<distance<<std::endl;
			#endif
			if(distance > 0) continue;
			bool canStitch = canIntroStitch(prec1,prec2);
			if(canStitch){
				#ifdef DEBUG_LIWEI_DONE
					std::cout<<"One stitsh found successfully"<<std::endl;
				#endif
				stitchCount ++;
				rec_stitch_list[i].push_back(j);
				rec_stitch_list[j].push_back(i);//if two recs with stitch relations locate in same poly, and warning should be arised
				continue;
			}
			else{
				if(rect_list[i].second == rect_list[j].second){
					std::cout<<"LIWEI: introStitch implementation is wrong"<<std::endl;
				}
			}
			add_edge(i, j, G);
		}
	}
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<new_polygon_id_list.size()<<" "<<num_vertices(G)<<", Graph build done"<<std::endl;
	#endif

	mplAssert(new_polygon_id_list.size() == num_vertices(G));
  	int component_num = connected_components(G, &new_polygon_id_list[0]);
	for(uint32_t i=0; i< new_polygon_id_list.size();i++){
		new_polygon_id_list[i] += start;
	}

	stitch_list.resize(component_num);
	polygon_id += component_num;
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<component_num<<", Graph calculate done"<<std::endl;
	#endif

	//transfrom the stitch list of rectangles into polys
	for (uint32_t i = 0; i < rect_list.size(); i++)
	{
		for (uint32_t j = i + 1; j < rect_list.size(); j++)
		{
			if(std::find(rec_stitch_list[i].begin(), rec_stitch_list[i].end(), j) != rec_stitch_list[i].end()){
				if(new_polygon_id_list[i] == new_polygon_id_list[j]){
					std::cout<<"LIWEI: warning: two recs with stitch relations locate in same poly: "<<start<<std::endl;
				}
				else{
					if(std::find(stitch_list[new_polygon_id_list[i] - start].begin(), stitch_list[new_polygon_id_list[i] - start].end(), new_polygon_id_list[j]) == stitch_list[new_polygon_id_list[i] - start].end())
						{
							stitch_list[new_polygon_id_list[i] - start].push_back(new_polygon_id_list[j]);    
							stitch_list[new_polygon_id_list[j] - start].push_back(new_polygon_id_list[i]);    
						}  
				}
			}
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
	StitchRelation.insert(StitchRelation.end(), stitch_list.begin(), stitch_list.end());
}	
void SimpleMPL::reconstruct_polygon(uint32_t& polygon_id, std::vector<uint32_t>& new_polygon_id_list, std::vector<std::pair<rectangle_pointer_type, uint32_t> >& rect_list)
{
	std::vector<std::vector<uint32_t> > stitch_list;
	uint32_t start = polygon_id + 1;
	std::vector<std::pair<rectangle_pointer_type, uint32_t> > rect_temp;
	std::vector<uint32_t> new_polygon_id_temp;

	uint32_t rect_list_size = rect_list.size();
	std::vector<bool> visited(rect_list_size, false);

	std::vector<uint32_t>().swap(new_polygon_id_list);
	new_polygon_id_list.assign(rect_list_size, std::numeric_limits<uint32_t>::max());

	// traverse all the rectangles to generate new IDs
	for (uint32_t i = 0; i< rect_list_size; i++)
	{
		if (new_polygon_id_list[i] == std::numeric_limits<uint32_t>::max())
		{
			if(visited[i])
				std::cout<<"VISITED POLYGON BUT ID IS MAX"<<i<<std::endl;
			polygon_id += 1; //TODO: here maybe a bug, += 1 should be behind
			new_polygon_id_list[i] = polygon_id;
		}
		visited[i] = true;
		for (uint32_t j = i + 1; j < rect_list_size; j++)
		{
			if (visited[j] == false && rect_list[i].second != rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				visited[j] = true;
				new_polygon_id_list[j] = new_polygon_id_list[i];
			}
		}
	}

#ifdef DEBUG_LIWEI_DONE
	for (uint32_t i = 0; i< rect_list_size; i++){
		for (uint32_t j = i + 1; j < rect_list_size; j++)
		{
			if (rect_list[i].second != rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				if(new_polygon_id_list[j] != new_polygon_id_list[i]){
					std::cout<<"ERROR: they should belong to same polygon: "<<rect_list[i].second<<" "<<rect_list[j].second<<std::endl;
					rectangle_pointer_type const & pPattern1 = m_db->polyrect_patterns()[rect_list[i].second];
					rectangle_pointer_type const & pPattern2 = m_db->polyrect_patterns()[rect_list[j].second];
					coordinate_type xl1 = gtl::xl(*pPattern1);
					coordinate_type yl1 = gtl::yl(*pPattern1);
					coordinate_type xl2 = gtl::xl(*pPattern2);
					coordinate_type yl2 = gtl::yl(*pPattern2);
					std::cout<<"ERROR: the x,y value are "<<xl1<<" "<<yl1<<" "<<xl2<<" "<<yl2<<std::endl;
			
				}
			}
		}
	}
#endif

	// for rectangles in the same original rectangle, we need to introduce a stitch between them, and record the stitch relationships.
	stitch_list.resize(polygon_id - start + 1);
	for (uint32_t i = 0; i < rect_list.size(); i++)
	{
		for (uint32_t j = i + 1; j < rect_list.size(); j++)
		{
			if (rect_list[i].second == rect_list[j].second && boost::geometry::distance(*(rect_list[i].first), *(rect_list[j].first)) == 0)
			{
				stitch_list[new_polygon_id_list[i] - start].push_back(new_polygon_id_list[j]);
				stitch_list[new_polygon_id_list[j] - start].push_back(new_polygon_id_list[i]);
			}
		}
	}
	// rectangles in the same polygon should be abutting
	for (uint32_t i = start, ie = polygon_id; i <= ie; ++i)
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
	StitchRelation.insert(StitchRelation.end(), stitch_list.begin(), stitch_list.end());
}

void SimpleMPL::splitRectangle(rectangle_type & pRect, std::vector<rectangle_pointer_type>& split, std::vector<rectangle_pointer_type> nei_Vec)
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
	else {
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
		
#ifdef DEBUG_LIWEI_DONE
		vInterSect.push_back(temp);
#else
		vInterSect.push_back(temp);
#endif 

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

	GenerateStitchPosition(pRect, vInterSect, vPossibleStitches, vstitches);

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

#ifdef STITCH_LEGAL
	if (hor)
	{
		for(int )
	}
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

void SimpleMPL::GenerateStitchPosition(const rectangle_type pRect, std::vector<rectangle_type> vInterSect,
	std::vector<coordinate_type>& vPossibleStitches, std::vector<coordinate_type>& vstitches)
{
	bool ishor = gtl::delta(pRect, gtl::HORIZONTAL) >= gtl::delta(pRect, gtl::VERTICAL);
	coordinate_type lower, upper;
	if (ishor) {
		lower = gtl::xl(pRect);
		upper = gtl::xh(pRect);
	}
	else {
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
			#ifdef DEBUG_LIWEI_DONE
				//std::cout<<"i is"<<i<<",vzeroids are "<< vzeroids << ", vStage is"<< vStages<<std::endl;
			#endif
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
		for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			int mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
			int diff = std::abs(vStages[i + 1].second - vStages[i - 1].second);
			double value = (double)mind + (double)diff*0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			coordinate_type pos = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
			vstitches.push_back(pos);
		}
	}
	std::sort(vstitches.begin(), vstitches.end());

#ifdef _DEBUG_PROJECTION
	if (lower == -2275 && upper == -2090)
	{
		std::cout << "\nDEBUG_PROJECTION| output the vzeroids : ";
		for (int i = 0; i < vzeroids.size(); i++)
			std::cout << vzeroids[i] << " ";
		std::cout << std::endl;
	}
	std::cout << "DEBUG_PROJECTION| output the vStages : ";
	for (uint32_t i = 0; i < vStages.size(); i++)
		std::cout << vStages[i].second << "[ " << vStages[i].first.first << "," << vStages[i].first.second << "] ;    ";
	std::cout << std::endl;

	std::cout << "DEBUG_PROJECTION| output the vstitches : ";
	for (uint32_t i = 0; i < vstitches.size(); i++)
		std::cout << vstitches[i] << " ";
	std::cout << std::endl;
#endif
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
	for (uint32_t i = 0; i < StitchRelation.size(); i++)
	{
		for(uint32_t j = 0; j < StitchRelation[i].size();j++){
			if(m_db->vPatternBbox[i]->color() != m_db->vPatternBbox[StitchRelation[i][j]]->color())
				stitch_count++;
		}
	}
#if RECORD > 0
	std::ofstream myfile;
	myfile.open ("record.txt", std::ofstream::app);
	myfile << "Stitch number: "<<stitch_count/2<<"\n";
	myfile << "Conflict number: "<<c_num<<"\n";
	myfile.close();
	std::ofstream result;
	if(m_db->use_stitch()){
		result.open("result_w_stitch.txt",std::ofstream::app);
		if(m_db->algo() == AlgorithmTypeEnum::ILP_GURBOI){
			result<<"\\\\\n"<<m_db->input_gds();
		}
		result<<"&"<<stitch_count/2<<"&"<<c_num<<"&"<<0.1*stitch_count/2 + c_num;
	}
	else{
		result.open("result_wo_stitch.txt",std::ofstream::app);
		if(m_db->algo() == AlgorithmTypeEnum::ILP_GURBOI){
			result<<"\\\\\n"<<m_db->input_gds();
		}
		result<<"&"<<c_num;
	}
	result.close();

#endif
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
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"number pf Edges with length zero are "<<zero_len_edge_num/2<<", m_db->coloring_distance is"<<m_db->coloring_distance<<", distance_value_count is "<< std::endl;
		for (std::vector<uint32_t>::iterator i = distance_value_count.begin(); i != distance_value_count.end(); ++i)
        	std::cout << *i << ' ';
	#endif
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
	std::cout<<"After component analysis, there are "<<m_comp_cnt<<" components in total."<<std::endl;
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
	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"For component "<<comp_id<<", there are "<<component_count<<" vertices in total."<<std::endl;
	#endif
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
    pcs->stitch_weight(m_db->parms.weight_stitch);
    pcs->color_num(m_db->color_num());
    pcs->threads(1); // we use parallel at higher level 

    return pcs;
}

void SimpleMPL::vdd_biconnected_component()
{

}

bool SimpleMPL::fast_color_trial(std::vector<int8_t>& vSubColor,SimpleMPL::graph_type const& sg){
	SimpleMPL::graph_type tmp_graph = sg;
	bool in_loop = true;
	std::vector<bool> in_graph;
	in_graph.assign(num_vertices(tmp_graph),true);
	std::stack<vertex_descriptor> S;
	bool debug = false;
	#ifdef DEBUG_LIWEI_DONE
		if(num_vertices(tmp_graph) == 5){
			debug = true;
			std::cout<<"LIWEI: fast color the 2 component!"<<std::endl;
		}
	#endif
	//std::cout<<"FST: graph enter done."<<std::endl;
	while(in_loop){
		in_loop = false; 
		boost::graph_traits<graph_type>::vertex_iterator vi, vie,next;
		boost::tie(vi, vie) = vertices(tmp_graph);
		for (next = vi; vi != vie; vi = next)
		{
			++next;
			vertex_descriptor v = *vi;
			mplAssert(uint32_t(v)<num_vertices(tmp_graph));
			if(in_graph[int(v)] == false)
				continue;
			int32_t conflictPreVDD_degree = 0;
			int32_t stitchPreVDD_degree = 0;
			//TODO: waiting for the response of LImbo
			if(debug)
				std::cout<<"LIWEI: now processing "<<v<<std::endl;
			boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
			boost::tie(vi2, vie2) = boost::adjacent_vertices(v, tmp_graph);
			for (next2 = vi2; vi2 != vie2; vi2 = next2)
			{
				++next2; 
				vertex_descriptor v2 = *vi2;
				mplAssert(uint32_t(v2)<num_vertices(tmp_graph));
				if(in_graph[int(v2)] == false)
					continue;
				if(debug)
					std::cout<<"liwei: PROCESS ADJ VERTEX "<<v2<<std::endl;
				// skip stitch edges 

				std::pair<edge_descriptor, bool> e12 = boost::edge(v, v2, tmp_graph);
				mplAssert(e12.second);
				if (boost::get(boost::edge_weight, tmp_graph, e12.first) < 0) 
					{	stitchPreVDD_degree += 1;
						continue;}
				conflictPreVDD_degree += 1;
				if(debug)
					std::cout<<"liwei: PROCESS ADJ VERTEX DONE"<<v2<<std::endl;
			}
			if(conflictPreVDD_degree < m_db->color_num() && stitchPreVDD_degree < 2)
			{
				S.push(v);
				in_graph[int(v)] = false;
				//boost::remove_vertex(v,tmp_graph);
				in_loop = true;
			if(debug)
				std::cout<<"LIWEI: now processing DONE, push "<<v<<std::endl;

			}	
	}}
	if(debug) std::cout<<"FST: graph remove done"<<std::count(in_graph.begin(),in_graph.end(),true)<<std::endl;
	if(std::count(in_graph.begin(),in_graph.end(),true) == 0){
		while(S.empty() == false){
			vertex_descriptor v = S.top();
			S.pop();
			if(debug)	std::cout<<"LIWEI: now calculate color of vertex "<<v<< ", this is "<<std::count(in_graph.begin(),in_graph.end(),true)<<std::endl;
			if(std::count(in_graph.begin(),in_graph.end(),true) == 0){
				in_graph[int(v)] = true;
				vSubColor[v] = 0;
				continue;
			}
			std::vector<bool> legal_color;
			legal_color.assign(m_db->color_num(),true);
			boost::graph_traits<SimpleMPL::graph_type>::adjacency_iterator vi2, vie2;
			
			for (boost::tie(vi2, vie2) = boost::adjacent_vertices(v, tmp_graph); vi2 != vie2; ++vi2)
			{
				if(in_graph[int(*vi2)] == false)	continue;
				// skip stitch edges 
				mplAssert(vSubColor[*vi2]>-1);
				std::pair<boost::graph_traits<SimpleMPL::graph_type>::edge_descriptor, bool> e12 = boost::edge(v, *vi2, tmp_graph);
				mplAssert(e12.second);
				if(debug)	std::cout<<"LIWEI: now calculate adj of vertex "<<v<<", adj is"<<*vi2<< ", the edge is "<<\
				(int)boost::get(boost::edge_weight, tmp_graph, e12.first)<<", adj color is"<<vSubColor[*vi2]<<std::endl;
				//For stitch edge, we assign the node color as same 
				if (boost::get(boost::edge_weight, tmp_graph, e12.first) < 0) 
					{	
						vSubColor[v] = vSubColor[*vi2];
					}
				else{
					legal_color[vSubColor[*vi2]] = false;
				} 
			} 
			//if the vSubColor still not set ( no stitch edge )
			if(vSubColor[v]>-1 && legal_color[vSubColor[v]] == true){
				mplAssert(vSubColor[v] > -1); 
			}
			else{
				for(uint32_t i = 0; i < legal_color.size();i++){
						if(legal_color[i] == true)
							{vSubColor[v] = i;break;}
					}
			}
			mplAssert(vSubColor[v] > -1);
			in_graph[int(v)] = true;
		}
		return true;
	}
	return false;
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
	//bool find_debug = false;
	// coordinate_difference hard_code = -3850;
	// coordinate_difference hard_code_y = -1730;

	// uint32_t hard_code_vertex = -1;
	// coordinate_difference hard_code2 = -3890;
	// coordinate_difference hard_code_y2 = -2230;
	 
	// uint32_t hard_code_vertex2 = -1;
	



	#ifdef DEBUG_LIWEI_DONE
		std::cout<<"LIWEI: for component "<<comp_id<<", size is"<<pattern_cnt<<", hiden size is"<<not_in_DG.size()<<std::endl;
	#endif
    // set max merge level, actually it only works when MERGE_SUBK4 is on 
    if (m_db->color_num() == 3)
        gs.max_merge_level(3);
    else if (m_db->color_num() == 4) // for 4-coloring, low level MERGE_SUBK4 works better 
        gs.max_merge_level(2);
    gs.simplify(simplify_strategy);
	
	// collect simplified information 
	std::stack<vertex_descriptor> vHiddenVertices = gs.hidden_vertices();
#ifdef DEBUG_LIWEI_DONE
    std::stack<vertex_descriptor> vHiddenVertices_debug = gs.hidden_vertices();
	while(!vHiddenVertices_debug.empty()){
		if(vHiddenVertices_debug.top() == hard_code_vertex){
			std::cout<<"LIWEI: is hiden!"<<std::endl;
			break;
		}
		if(vHiddenVertices_debug.top() == hard_code_vertex2){
			std::cout<<"LIWEI:2 is hiden!"<<std::endl;
			break;
		}
		vHiddenVertices_debug.pop();
	}
#endif

    // for debug, it does not affect normal run 
	// if (comp_id == m_db->dbg_comp_id() && simplify_strategy != graph_simplification_type::MERGE_SUBK4)
	// 	gs.write_simplified_graph_dot("graph_simpl");

	// in order to recover color from articulation points 
	// we have to record all components and mappings 
	// but graph is not necessary 
	std::vector<std::vector<int8_t> > mSubColor (gs.num_component());
	std::vector<std::vector<vertex_descriptor> > mSimpl2Orig (gs.num_component());
	DG_num += gs.num_component();
	double acc_obj_value = 0;

	std::vector<vertex_descriptor> all_articulations;

	gs.get_articulations(all_articulations);
#ifdef DEBUG_LIWEI_DONE
	if(all_articulations.size() > 0)
		std::cout<<"LIWEI: articulations size is "<<all_articulations.size()<<", component number is "<<gs.num_component()<<std::endl;
#endif
	
	for (uint32_t sub_comp_id = 0; sub_comp_id < gs.num_component(); ++sub_comp_id)
	{
		std::cout << std::endl;
		graph_type sg;
		std::vector<int8_t>& vSubColor =  mSubColor[sub_comp_id];
		std::vector<vertex_descriptor>& vSimpl2Orig = mSimpl2Orig[sub_comp_id];

		gs.simplified_graph_component(sub_comp_id, sg, vSimpl2Orig);
		
		vSubColor.assign(num_vertices(sg), -1);

#ifdef _OPENMP
#pragma omp critical(dgGlobal2Local)
#endif
		{
		for (std::map<uint32_t, uint32_t>::iterator it = dgGlobal2Local.begin(); it != dgGlobal2Local.end(); it++)
		{
			if (std::find(vSimpl2Orig.begin(), vSimpl2Orig.end(), it->second) != vSimpl2Orig.end())
			{
				dgCompId[it->first] = sub_comp_id + 1;
			}
		}
		}

		/***
		if(fast_color_trial(vSubColor,sg))
		{
			#ifdef DEBUG_LIWEI_DONE

				for(int i = 0; i< vSimpl2Orig.size();i++){
					if(vSimpl2Orig[i] == hard_code_vertex)
						std::cout<<"LIWEI: vertex 1: color by fast color trial "<<vSubColor[i]<<std::endl;
					if(vSimpl2Orig[i] == hard_code_vertex2)
					{
						for(int j = 0; j<vSimpl2Orig.size();j++){
							std::cout<<vSimpl2Orig[j]<<" ";
						}
						std::cout<<"LIWEI: vertex 2: color by fast color trial "<<vSubCol+or[i]<<std::endl;
					}
				}
				if(find_debug){
					std::cout<<"LIWEI: solved by fast color trial "<<sub_comp_id<<std::endl;
				}
			#endif
			continue;
		}
		
		#ifdef DEBUG_LIWEI_DONE
			std::cout<<"LIWEI: fct but failed"<<std::endl;
		#endif
		***/
    // for debug, it does not affect normal run 
		

		//if algorithm is Dancing Link, call it directly
		if(m_db->algo() == AlgorithmTypeEnum::DANCING_LINK){
			double cost = 0;
			if(m_db->use_stitch()){
				//cost = solve_by_dancing_link_with_stitch(sg,vSubColor, comp_id);
				cost = solve_by_dancing_link_with_one_stitch(sg,vSubColor, comp_id);
				//cost = solve_by_dancing_link_GPU(sg,vSubColor);
				}
			else{ cost = solve_by_dancing_link_with_stitch(sg,vSubColor, comp_id);}
	#if RECORD > 1
		if (comp_id == m_db->dbg_comp_id())
			{write_graph(sg, std::to_string(comp_id) +"_"+ std::to_string(sub_comp_id));
				std::ofstream myfile;
				myfile.open ("debug_component_color_results.txt", std::ofstream::app);
				myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<", obj_value: "<< cost<<"\n";
				for(uint32_t v = 0; v<vSubColor.size();v++){
					myfile<<"vertex "<<v<<" color is "<<(uint32_t)vSubColor[v]<<". \n";
				}
				myfile.close();
			}

			// std::ofstream myfile;
			// myfile.open ("small_results.txt", std::ofstream::app);
			// myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<", obj_value: "<< cost<<"\n";
			// myfile.close();
		if(cost != 0){
			std::ofstream myfile;
			myfile.open ("small_results.txt", std::ofstream::app);
			myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<", obj_value: "<< cost<<"\n";
			myfile.close();
		}

	#endif
				if (comp_id == m_db->dbg_comp_id()){
					for(uint32_t v = 0; v<vSubColor.size();v++){
						std::cout<<"vertex "<<v<<" color is "<<(uint32_t)vSubColor[v]<<". "<<std::endl;
					}
				}
			continue;
		}
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
		
		// if(obj_value1 != 0){
		// 	this->writeGraph(sg,std::to_string(comp_id),obj_value1);
		// 	std::cout<<"Write Component "<<comp_id<<", Sub component "<<sub_comp_id<<std::endl;
		// }

// #ifdef DEBUG
//         mplPrint(kDEBUG, "comp_id = %u, sub_comp_id = %u, %lu vertices, obj_value1 = %g\n", comp_id,sub_comp_id, num_vertices(sg), obj_value1); 
// #endif
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
                mplAssert(color >= 0 && color < m_db->color_num());
                vSubColor[v] = color;
				if (comp_id == m_db->dbg_comp_id()){
					std::cout<<"vertex "<<v<<" color is "<<(uint32_t)color<<". "<<std::endl;
				}
            }
        }
        else // no need to update vSubColor, as it is already updated by sub call 
            acc_obj_value += obj_value2;
			if (comp_id == m_db->dbg_comp_id()){
				for(uint32_t v = 0; v<vSubColor.size();v++){
					std::cout<<"vertex "<<v<<" color is "<<(uint32_t)vSubColor[v]<<". "<<std::endl;
				}
			}

#if RECORD > 1
	if(obj_value1 != 0){
		std::ofstream myfile;
		myfile.open ("small_results.txt", std::ofstream::app);
		myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<", obj_value: "<< obj_value1<<"\n";
		myfile.close();
	}
	if (comp_id == m_db->dbg_comp_id())
		{write_graph(sg, std::to_string(comp_id) +"_"+ std::to_string(sub_comp_id));
			std::ofstream myfile;
			myfile.open ("debug_component_color_results.txt", std::ofstream::app);
			myfile <<", comp_id: "<< comp_id<<", sub_comp_id: "<<sub_comp_id <<",num_vertices(sg) "<<num_vertices(sg)<<"\n";
			for(uint32_t v = 0; v<vSubColor.size();v++){
				myfile<<"vertex "<<v<<" color is "<<(uint32_t)vSubColor[v]<<". \n";
			}
			myfile.close();
		}

#endif	
		delete pcs;
	}

#ifdef _DGOUT
	return 1;
#endif
	// #ifdef DEBUG_LIWEI_DONE
	// if(find_debug)
	// 	std::cout<<"Color before recover is "<<(int)vColor[hard_code_vertex]<<" "<<(int)vColor[hard_code_vertex2]<<std::endl;
	// #endif
	// recover color assignment according to the simplification level set previously 
	// HIDE_SMALL_DEGREE needs to be recovered manually for density balancing 

	gs.recover(vColor, mSubColor, mSimpl2Orig); 

	// #ifdef DEBUG_LIWEI_DONE
	// if(find_debug)
	// 	std::cout<<"Color is "<<(int)vColor[hard_code_vertex]<<" "<<(int)vColor[hard_code_vertex2]<<std::endl;
	// #endif
	// recover colors for simplified vertices with balanced assignment 
	// recover hidden vertices with local balanced density control 
    RecoverHiddenVertexDistance(dg, itBgn, pattern_cnt, vColor, vHiddenVertices, m_vColorDensity, *m_db)();

	//mplPrint(kINFO, "Component %u solved", comp_id);

	//std::string name = std::to_string(comp_id) + "-" + std::to_string(sub_comp_id);
	//std::string name = std::to_string(comp_id);
	//this->writeJson(dg,(char*)name.c_str(),vColor);
    return acc_obj_value;
}


// void SimpleMPL::solve_by_dancing_link(SimpleMPL::graph_type& g,std::vector<int8_t>& color_vector){
// 	//Due to graph does not contain a parent&child node system
// 	// we simply redesign a node struct to achieve this kind of system
// 	mplAssert(num_vertices(g) == color_vector.size());
// 	DancingLink dl; 
// 	std::vector<Vertex*> node_list;
// 	node_list.resize(num_vertices(g));
// 	uint32_t vertex_numbers = 0;
// 	uint32_t edge_numbers = 0;
// 	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
// 	for (boost::tie(vi1, vie1) = boost::vertices(g); vi1 != vie1; ++vi1)
// 	{	
// 		vertex_descriptor v1 = *vi1;
// 		boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
// 		boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, g);
// 		for (next2 = vi2; vi2 != vie2; vi2 = next2)
// 		{
// 			++next2; 
// 			vertex_descriptor v2 = *vi2;
// 			if (v1 >= v2) continue;
// 			std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, g);
// 			mplAssert(e12.second);
			
// 			//if node source and node target is not created. New one/
// 			if(!node_list[v1]){
// 				Vertex* source_vertex = new Vertex;
// 				mplAssert(source_vertex->Conflicts.empty());
// 				node_list[v1] = source_vertex;
// 				source_vertex->Stitch_No = v1;
// 				vertex_numbers ++;
// 			}
// 			if(!node_list[v2]){
// 				Vertex* target_vertex = new Vertex;
// 				target_vertex->Stitch_No = v2;
// 				node_list[v2] = target_vertex;
// 				vertex_numbers ++;
// 			}

// 			//if two nodes are stitch relationships and both of them are parent
// 			if(boost::get(boost::edge_weight, g, e12.first) < 0){
// 				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
// 				if(node_list[v1]->parent != NULL && node_list[v1]->parent== node_list[v2]->parent) continue;
// 				Vertex* parent_vertex = new Vertex;
// 				mplAssert(parent_vertex->Childs.empty());
// 				parent_vertex->parentOf(node_list[v1]);
// 				parent_vertex->parentOf(node_list[v2]);
// 				vertex_numbers -= 1;
// 				continue;
// 			}
// 			node_list[v1]->Conflicts.insert(node_list[v2]);
// 		}
// 	}
// 	//need to store node list, which does not contain child node
// 	std::set<Vertex*> node_wo_stitch_list;
// 	int index = 1;
// 	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
// 		//if the node in parent node, means it has no stitch relations
// 		(*it)->updateConflicts();
// 		if((*it)->Is_Parent){
// 			node_wo_stitch_list.insert((*it));
// 			(*it)->No = index;
// 			index ++;
// 		}
// 		//else, add its parent node if it has not been added into node_wo_stitch_list
// 		else{
// 			if((*it)->parent->Is_Parent == false){
// 				std::cout<<(*it)->Stitch_No<<std::endl;
// 				std::cout<<(*it)->parent->Stitch_No<<std::endl;
// 			}
// 			mplAssert((*it)->parent->Is_Parent);
// 			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
// 				node_wo_stitch_list.insert((*it)->parent);
// 				(*it)->parent->updateConflicts();
// 				(*it)->parent->No = index;
// 				index ++;
// 			}

// 		}
// 	}
// 	mplAssert(node_wo_stitch_list.size() == vertex_numbers);


// 	std::vector<std::list<Edge_Simple> >   edge_list;
// 	edge_list.resize(vertex_numbers + 1);
// 	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it) {
// 		mplAssert((*it)->No != 0);
// 		for(std::set<Vertex*>::iterator itconflict = (*it)->Conflicts_in_LG.begin(); itconflict != (*it)->Conflicts_in_LG.end(); ++itconflict) {
// 			if((*itconflict)->Is_Parent){
// 				mplAssert((*itconflict)->No != 0);
// 				edge_numbers++;
// 				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->No, (int)edge_numbers });
// 				edge_list[(*itconflict)->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });	
// 			}
// 			else{
// 				mplAssert((*itconflict)->parent->Is_Parent);
// 				//avoid repeat (conflicts of stitch realation nodes may represent same conflict)
// 				edge_numbers++;	
// 				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->parent->No, (int)edge_numbers });
// 				edge_list[(*itconflict)->parent->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });
// 			}
// 		}
// 		}
// 	std::cout<<"EDGE list generated with size: "<<edge_numbers<<std::endl;

	
	
// 	int row_numbers = vertex_numbers * m_db->color_num() + 1;
// 	int col_numbers = edge_numbers * m_db->color_num() + vertex_numbers;
// 	DL_Init(dl, row_numbers, col_numbers);

// 	// Insert DL cells
// 	for (uint32_t it = 1; it < edge_list.size(); ++it)
// 	{
// 		for (int32_t i = 1; i <= m_db->color_num(); ++i)
// 		{
// 			Cell_Insert(dl,(it - 1)*m_db->color_num() + i,it);
// 			for (auto j = edge_list[it].begin(); j != edge_list[it].end(); ++j)
// 			{
// 				Cell_Insert(dl,(it - 1)*m_db->color_num() + i,vertex_numbers + (j->No - 1)*m_db->color_num() + i);
// 			}
// 		}
// 	}
// 	for (uint32_t i = 0; i < edge_numbers; ++i)
// 	{
// 		for (int32_t j = 1; j <= m_db->color_num(); ++j)
// 			Cell_Insert(dl,vertex_numbers * m_db->color_num() + 1,vertex_numbers + i * m_db->color_num() + j);
// 	}

// 	// construct verctors for calling API in DL_MPL.cpp
// 	std::vector<int> result_vec;
// 	std::pair<int, int>  conflict_pair;
// 	std::vector<int> Delete_the_Row_in_which_Col;
// 	std::vector<std::list<int> >  Order_of_Row_Deleted_in_Col;
// 	std::vector<int> MPLD_search_vector;
// 	std::vector<bool> col_cover_vector(vertex_numbers);
// 	col_cover_vector.assign(vertex_numbers,false);
// 	std::vector<int> row_select_vector;
// 	MPLD_search_vector = Simple_Order(edge_list.size());
// 	int depth = 1;
// 	std::vector<int8_t> color_results_wo_stitch (vertex_numbers);
// 	color_results_wo_stitch.assign(vertex_numbers,-1);
// 	Order_of_Row_Deleted_in_Col.resize(col_numbers + 1);
// 	Delete_the_Row_in_which_Col.resize(row_numbers + 1);
// 	bool result = MPLD_X_Solver(dl,color_results_wo_stitch, result_vec, conflict_pair, vertex_numbers, m_db->color_num(), Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, MPLD_search_vector, "decode_file",col_cover_vector,row_select_vector);
// 	mplAssert(result);
// 	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
// 		if((*it)->No == 0){
// 			mplAssert((*it)->Is_Parent == false);
// 			mplAssert((*it)->parent->Stitch_No == -1);
// 			mplAssert(color_results_wo_stitch[(*it)->parent->No -1]!= -1);
// 			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->parent->No -1];
// 		}
// 		else{
// 			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->No -1];
// 		}
// 	}
// }

void SimpleMPL::printGraph(SimpleMPL::graph_type& g){
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
			std::cout<<v1<<" "<<v2<<" "<<boost::get(boost::edge_weight, g, e12.first)<<std::endl;
		}
	}
}
double SimpleMPL::solve_by_dancing_link_with_stitch(SimpleMPL::graph_type& g,std::vector<int8_t>& color_vector, uint32_t comp_id ){
	// Due to graph does not contain a parent&child node system
	// we simply redesign a node struct to achieve this kind of system
	mplAssert(num_vertices(g) == color_vector.size());
	DancingLink dl; 
	std::vector<Vertex*> node_list;
	node_list.resize(num_vertices(g));
	uint32_t vertex_numbers = 0;
	uint32_t edge_numbers = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
			//std::cout<<v1<<" "<<v2<<" "<<boost::get(boost::edge_weight, g, e12.first)<<std::endl;
			//if node source and node target is not created. New one/
			if(!node_list[v1]){
				Vertex* source_vertex = new Vertex;
				mplAssert(source_vertex->Conflicts.empty());
				node_list[v1] = source_vertex;
				source_vertex->Stitch_No = v1;
				vertex_numbers ++;
			}
			if(!node_list[v2]){
				Vertex* target_vertex = new Vertex;
				target_vertex->Stitch_No = v2;
				node_list[v2] = target_vertex;
				vertex_numbers ++;
			}

			//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0){
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				if(node_list[v1]->parent != NULL && node_list[v1]->parent== node_list[v2]->parent) continue;
				Vertex* parent_vertex = new Vertex;
				mplAssert(parent_vertex->Childs.empty());
				parent_vertex->parentOf(node_list[v1]);
				parent_vertex->parentOf(node_list[v2]);
				vertex_numbers -= 1;
				continue;
			}
			node_list[v1]->Conflicts.insert(node_list[v2]);
		}
	}
	//need to store node list, which does not contain child node
	std::set<Vertex*> node_wo_stitch_list;
	int index = 1;
	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateConflicts();
		if((*it)->Is_Parent){
			node_wo_stitch_list.insert((*it));
			(*it)->No = index;
			index ++;
		}
		//else, add its parent node if it has not been added into node_wo_stitch_list
		else{
			if((*it)->parent->Is_Parent == false){
				std::cout<<(*it)->Stitch_No<<std::endl;
				std::cout<<(*it)->parent->Stitch_No<<std::endl;
			}
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateConflicts();
				(*it)->parent->No = index;
				index ++;
			}

		}
	}
	mplAssert(node_wo_stitch_list.size() == vertex_numbers);

	std::vector<std::list<Edge_Simple> >   edge_list;
	edge_list.resize(vertex_numbers + 1);
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it) {
		mplAssert((*it)->No != 0);
		for(std::set<Vertex*>::iterator itconflict = (*it)->Conflicts_in_LG.begin(); itconflict != (*it)->Conflicts_in_LG.end(); ++itconflict) {
			if((*itconflict)->Is_Parent){
				mplAssert((*itconflict)->No != 0);
				edge_numbers++;
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->No, (int)edge_numbers });
				edge_list[(*itconflict)->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });	
			}
			else{
				mplAssert((*itconflict)->parent->Is_Parent);
				//avoid repeat (conflicts of stitch realation nodes may represent same conflict)
				edge_numbers++;	
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->parent->No, (int)edge_numbers });
				edge_list[(*itconflict)->parent->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });
			}
		}
		}
	std::cout<<"EDGE list generated with size: "<<edge_numbers<<std::endl;
	std::vector<uint32_t> starting_indexes;
	uint32_t starting_index = 0;

	//Stitch Generation: recording strating index of each parent node block
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it){
		uint32_t childs_num = ((*it)->Childs).size();
		mplAssert(childs_num > 1 || childs_num == 0);
		starting_indexes.push_back(starting_index);
		if(childs_num > 1){
			starting_index += pow(m_db->color_num(),childs_num);
		}
	}

	//For debug
	if(comp_id == m_db->dbg_comp_id()){
		for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
			if((*it)->No == 0){
				std::cout<<"originial node "<<(*it)->Stitch_No<<" no-stitch graph node "<<(*it)->parent->No<<std::endl;
			}
			else{
				std::cout<<"originial node "<<(*it)->Stitch_No<<" no-stitch graph node "<<(*it)->No<<std::endl;
			}
		}
	}
	//Node structure done,  then build DL system
	int row_numbers = vertex_numbers * m_db->color_num() + 1 + starting_index;
	int col_numbers = edge_numbers * m_db->color_num() + vertex_numbers;
	std::cout<<"row/col is"<<row_numbers<< " "<<col_numbers<<std::endl;
	std::cout<<"vertex/edge number is"<<vertex_numbers<< " "<<edge_numbers<<std::endl;
	DL_Init(dl, row_numbers, col_numbers);

	// Insert DL cells
	
	for (uint32_t it = 1; it < edge_list.size(); ++it)
	{
		for (int32_t i = 1; i <= m_db->color_num(); ++i)
		{
			Cell_Insert(dl,(it - 1)*m_db->color_num() + i,it);
			for (auto j = edge_list[it].begin(); j != edge_list[it].end(); ++j)
			{
				Cell_Insert(dl,(it - 1)*m_db->color_num() + i,vertex_numbers + (j->No - 1)*m_db->color_num() + i);
			}
		}
	}
	for (uint32_t i = 0; i < edge_numbers; ++i)
	{
		for (int32_t j = 1; j <= m_db->color_num(); ++j)
			Cell_Insert(dl,vertex_numbers * m_db->color_num() + 1,vertex_numbers + i * m_db->color_num() + j);
	}


	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateDuplicateLGConflicts();
		if(!(*it)->Is_Parent)
		//else, add its parent node if it has not been added into node_wo_stitch_list
		{
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateDuplicateLGConflicts();
			}

		}
	}

	//Stitch Generation Step 2: Cell insertion
	uint32_t count = 0;
	for(std::set<Vertex*>::iterator v = node_wo_stitch_list.begin(); v != node_wo_stitch_list.end(); ++v){
		uint32_t childs_num = ((*v)->Childs).size();
		//std::cout<<count<<"th node with child num "<<childs_num<<std::endl;
		uint32_t child_count = 1;
		for(std::vector<Vertex*>::iterator child = ((*v)->Childs).begin(); child != ((*v)->Childs).end(); ++child){
			//c_conflit = child conflict
			for(std::set<Vertex*>::iterator c_conflit = ((*child)->Conflicts_in_LG).begin(); c_conflit != ((*child)->Conflicts_in_LG).end(); ++c_conflit){
				for(int i = 1;i<=pow(m_db->color_num(),childs_num);i++){
					int col_edge = -1;
					for(std::list<Edge_Simple>::iterator conflict_edge = edge_list[(*child)->parent->No].begin();conflict_edge != edge_list[(*child)->parent->No].end();++conflict_edge){
						if((*conflict_edge).target == (*c_conflit)->No){
							col_edge = (*conflict_edge).No;
							break;
						}
					}
					mplAssert(col_edge!=-1);
					//std::cout<<starting_indexes[count] + i <<" "<< 1 + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))<<std::endl;
					Cell_Insert(dl,vertex_numbers * m_db->color_num() + 1 + starting_indexes[count] + i, vertex_numbers + 1 + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num())));
				} 
			}
			child_count ++;
		}
		count++;
	}

	//Stitch Generation Step 3: Decode matrix generation. pair: 1st = stitch No, 2nd = color  
	count = -1;
	std::vector<std::vector<std::pair<uint32_t,uint32_t>>> decode_mat;
	for(std::set<Vertex*>::iterator v = node_wo_stitch_list.begin(); v != node_wo_stitch_list.end(); ++v){
		uint32_t childs_num = ((*v)->Childs).size();
		count++;
		if(childs_num == 0){continue;}
		for(int i = 1;i<=pow(m_db->color_num(),childs_num);i++){
			//std::cout<<i<<" "<<(*v)->No<<" "<<childs_num<<" "<<pow(m_db->color_num(),childs_num)<<std::endl;
			std::vector<std::pair<uint32_t,uint32_t>> row_decoder;
			uint32_t child_count =1;
			for(std::vector<Vertex*>::iterator child = ((*v)->Childs).begin(); child != ((*v)->Childs).end(); ++child){
				//std::cout<<(*child)->Stitch_No<<" "<<(i/(int)pow(m_db->color_num(),childs_num-count)%(int)(m_db->color_num())) <<std::endl;
				row_decoder.push_back(std::make_pair((*child)->Stitch_No,((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))));
				child_count ++;
			}
			//std::cout<<vertex_numbers * m_db->color_num() + 1 + starting_indexes[count] + i<<" "<<(*v)->No <<std::endl;
			Cell_Insert(dl,vertex_numbers * m_db->color_num() + 1 + starting_indexes[count] + i,(*v)->No);
			decode_mat.push_back(row_decoder);
		}
	}

	// construct verctors for calling API in DL_MPL.cpp
	//result_vec: row results of this iteration, should be empty if no optimal result (without conflict)
	//row_select_vector: partial row results, should be empty if there is optimal result
	//partial_last_rows: last removed rows of each column(node) 
	std::vector<int> result_vec;
	std::pair<int, int>  conflict_pair;
	std::set<std::pair<int, int> >  final_conflict;
	std::vector<int> final_result;
	std::vector<int> partial_last_rows;
	std::vector<int> last_rows;
	std::vector<int> partial_conflict_col_table;
	std::vector<int> conflict_col_table;
	std::vector<int> Delete_the_Row_in_which_Col;
	std::vector<std::list<int> >  Order_of_Row_Deleted_in_Col;
	std::vector<int> MPLD_search_vector;
	std::vector<bool> col_cover_vector(vertex_numbers);
	col_cover_vector.assign(vertex_numbers,false);
	std::vector<int> row_select_vector;
	MPLD_search_vector = BFS_Order(edge_list);
	int depth = 1;
	std::vector<int8_t> color_results_wo_stitch (vertex_numbers);
	color_results_wo_stitch.assign(vertex_numbers,-1);
	partial_conflict_col_table.assign(col_numbers + 1, 0);
	conflict_col_table.assign(col_numbers + 1, 0);
	partial_last_rows.assign(col_numbers + 1, 0);
	last_rows.assign(col_numbers + 1, 0);
	Order_of_Row_Deleted_in_Col.resize(col_numbers + 1);
	Delete_the_Row_in_which_Col.resize(row_numbers + 1);
	bool need_debug = false;
	if(comp_id == m_db->dbg_comp_id()){
		std::cout<<"debug comp_id"<< std::endl;
		need_debug = true;
	}
	// bool result = MPLD_X_Solver(dl,color_results_wo_stitch, result_vec, conflict_pair, (int)vertex_numbers, (int)m_db->color_num(), Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, 
	// depth, MPLD_search_vector, "decode_file",col_cover_vector,row_select_vector,partial_conflict_col_table, conflict_col_table,partial_last_rows,last_rows,need_debug);
	//result == false means that there is a exact conflict, then we will remove corresponding edge columns 
	// while(result==false){

	// 	final_conflict.insert(conflict_pair);
	// 	int col_edge = -1;
	// 	for(std::list<Edge_Simple>::iterator conflict_edge = edge_list[conflict_pair.first].begin();conflict_edge != edge_list[conflict_pair.first].end();++conflict_edge){
	// 		if((*conflict_edge).target == conflict_pair.second){
	// 			col_edge  = (*conflict_edge).No;
	// 			break;
	// 		}
	// 	}
	// 	if(col_edge == -1){
	// 		std::cout<<"bug found!"<<std::endl;
	// 	}
	// 	mplAssert(col_edge!=-1);
	// 	for(int i = 1;i<=m_db->color_num();i++){
	// 		int this_col = vertex_numbers + (col_edge -1)*(m_db->color_num()) + i;
	// 		//Cell *col = &dl.Col_Header_Table[this_col];
	// 		Remove_Single_Col(dl,this_col);
	// 	}
	// 	final_result.insert( final_result.end(), row_select_vector.begin(), row_select_vector.end() );
		
	// 	//NOTE that the order between the two parts (remove edges and recover partial results are important! Cause some rows won't be removed this time due to the exact conflict edge is removed)
	// 	for(auto i = row_select_vector.begin(); i != row_select_vector.end(); i++){
	// 		//if the row is selected ( in the partial result), recover the partial results
	// 			std::set<int> row_set;
	// 			std::set<int> col_set;
	// 			Select_All_Rows_Cols(dl, (*i), row_set, col_set);
	// 			store_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col,conflict_col_table,last_rows);
	// 			Remove_Rows_Cols(dl, row_set, col_set);

	// 	}
	// 	row_select_vector.clear();
	// 	Delete_the_Row_in_which_Col.clear();
	// 	Order_of_Row_Deleted_in_Col.clear();
	// 	partial_conflict_col_table.swap(conflict_col_table);
	// 	partial_last_rows.swap(last_rows);
	// 	depth = 1;
	// 	result = MPLD_X_Solver(dl,color_results_wo_stitch, result_vec, conflict_pair, vertex_numbers, m_db->color_num(), 
	// 	Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, 
	// 	MPLD_search_vector, "decode_file",col_cover_vector,row_select_vector,
	// 	partial_conflict_col_table, conflict_col_table,partial_last_rows,last_rows,need_debug);
	// 	//std::cout<<starting_indexes[count] + i <<" "<< 1 + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))<<std::endl;
	// 	//remove corresponding edge of conflict
	// }
	conflict_col_table.clear();
	bool result = Efficient_MPLD_X_Solver(dl,color_results_wo_stitch, result_vec, conflict_pair, (int)vertex_numbers, (int)m_db->color_num(), Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, 
	depth, MPLD_search_vector, "decode_file",row_select_vector,partial_conflict_col_table, conflict_col_table,need_debug);
	int itration = 1;
	if(comp_id == m_db->dbg_comp_id()){
		std::cout<<"debug comp_id"<< std::endl;
	}
	while(result==false && row_select_vector.size() < vertex_numbers){
		itration++;
		std::cout<<"itration "<<itration<<", comp_id "<<comp_id<<std::endl;
		final_conflict.insert(conflict_pair);
		int col_edge = -1;
		for(std::list<Edge_Simple>::iterator conflict_edge = edge_list[conflict_pair.first].begin();conflict_edge != edge_list[conflict_pair.first].end();++conflict_edge){
			if((*conflict_edge).target == conflict_pair.second){
				col_edge  = (*conflict_edge).No;
				break;
			}
		}
		if(col_edge==-1){
			std::cout<<"bug found!"<<std::endl;
		}
		mplAssert(col_edge!=-1);
		mplAssert(result_vec.size() == 0);
		for(int i = 1;i<=m_db->color_num();i++){
			int this_col = vertex_numbers + (col_edge -1)*(m_db->color_num()) + i;
			//Cell *col = &dl.Col_Header_Table[this_col];
			Remove_Single_Col(dl,this_col);
		}
		final_result.insert( final_result.end(), row_select_vector.begin(), row_select_vector.end() );
		//NOTE that the order between the two parts (remove edges and recover partial results are important! Cause some rows won't be removed this time due to the exact conflict edge is removed)
		for(int i = 0; i < row_select_vector.size(); i++){
			//if the row is selected ( in the partial result), recover the partial results
				std::set<int> row_set;
				std::set<int> col_set;
				Select_All_Rows_Cols(dl, row_select_vector[i], row_set, col_set);
				efficient_store_intermediate_process(dl, partial_conflict_col_table[i], row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
				Remove_Rows_Cols(dl, row_set, col_set);

		}
		row_select_vector.clear();
		partial_conflict_col_table.swap(conflict_col_table);
		conflict_col_table.clear();
		depth = 1;
		result = Efficient_MPLD_X_Solver(dl,color_results_wo_stitch, result_vec, conflict_pair, vertex_numbers, m_db->color_num(), 
		Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, 
		MPLD_search_vector, "decode_file",row_select_vector,
		partial_conflict_col_table, conflict_col_table,need_debug);
		//std::cout<<starting_indexes[count] + i <<" "<< 1 + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))<<std::endl;
		//remove corresponding edge of conflict
	}
	if(result_vec.size() != 0){
		final_result.insert( final_result.end(), result_vec.begin(), result_vec.end() );
		mplAssert(final_result.size() == vertex_numbers);
	}
	else final_result.insert( final_result.end(), row_select_vector.begin(), row_select_vector.end() );
	for (auto i = final_result.begin(); i != final_result.end(); i++)
	{
		// std::cout << *i << std::endl;
		if((*i) < vertex_numbers * m_db->color_num() + 1){
		int No = ((*i) - 1) / m_db->color_num() + 1;
		int mask = (*i + 2) % m_db->color_num();
		color_results_wo_stitch[No-1] = mask; }
		if((*i) > (vertex_numbers * (uint32_t)(m_db->color_num()) + (uint32_t)1)){
			std::vector<std::pair<uint32_t,uint32_t>>& row_decoder = decode_mat[(*i)- vertex_numbers * m_db->color_num() -2];
			for(std::vector<std::pair<uint32_t,uint32_t>>::iterator it = row_decoder.begin(); it != row_decoder.end(); ++it) {
				color_vector[(*it).first] = (*it).second;
			}
		}
	}


	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		if(color_vector[(*it)->Stitch_No]!= -1){continue;}
		if((*it)->No == 0){
			mplAssert((*it)->Is_Parent == false);
			mplAssert((*it)->parent->Stitch_No == -1);
			mplAssert(color_results_wo_stitch[(*it)->parent->No -1]!= -1);
			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->parent->No -1];
		}
		else{
			mplAssert(color_results_wo_stitch[(*it)->No -1]!=-1);
			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->No -1];
		}
	}
	double cost = new_calc_cost(g,color_vector);
	std::cout<<"cost is "<<cost<<std::endl;
	return cost;
}




double SimpleMPL::solve_by_dancing_link_with_one_stitch(SimpleMPL::graph_type& g,std::vector<int8_t>& color_vector, uint32_t comp_id ){
	// Due to graph does not contain a parent&child node system
	// we simply redesign a node struct to achieve this kind of system
	mplAssert(num_vertices(g) == color_vector.size());
	DancingLink dl; 
	std::vector<Vertex*> node_list;
	node_list.resize(num_vertices(g));
	uint32_t vertex_numbers = 0;
	uint32_t edge_numbers = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
			//std::cout<<v1<<" "<<v2<<" "<<boost::get(boost::edge_weight, g, e12.first)<<std::endl;
			//if node source and node target is not created. New one/
			if(!node_list[v1]){
				Vertex* source_vertex = new Vertex;
				mplAssert(source_vertex->Conflicts.empty());
				node_list[v1] = source_vertex;
				source_vertex->Stitch_No = v1;
				vertex_numbers ++;
			}
			if(!node_list[v2]){
				Vertex* target_vertex = new Vertex;
				target_vertex->Stitch_No = v2;
				node_list[v2] = target_vertex;
				vertex_numbers ++;
			}

			//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0){
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				if(node_list[v1]->parent != NULL && node_list[v1]->parent== node_list[v2]->parent) continue;
				Vertex* parent_vertex = new Vertex;
				mplAssert(parent_vertex->Childs.empty());
				parent_vertex->parentOf(node_list[v1]);
				parent_vertex->parentOf(node_list[v2]);
				vertex_numbers -= 1;
				continue;
			}
			node_list[v1]->Conflicts.insert(node_list[v2]);
		}
	}
	//need to store node list, which does not contain child node
	std::set<Vertex*> node_wo_stitch_list;
	int index = 1;
	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateConflicts();
		if((*it)->Is_Parent){
			node_wo_stitch_list.insert((*it));
			(*it)->No = index;
			index ++;
		}
		//else, add its parent node if it has not been added into node_wo_stitch_list
		else{
			if((*it)->parent->Is_Parent == false){
				std::cout<<(*it)->Stitch_No<<std::endl;
				std::cout<<(*it)->parent->Stitch_No<<std::endl;
			}
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateConflicts();
				(*it)->parent->No = index;
				index ++;
			}

		}
	}
	mplAssert(node_wo_stitch_list.size() == vertex_numbers);

	std::vector<std::list<Edge_Simple> >   edge_list;
	edge_list.resize(vertex_numbers + 1);
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it) {
		mplAssert((*it)->No != 0);
		for(std::set<Vertex*>::iterator itconflict = (*it)->Conflicts_in_LG.begin(); itconflict != (*it)->Conflicts_in_LG.end(); ++itconflict) {
			if((*itconflict)->Is_Parent){
				mplAssert((*itconflict)->No != 0);
				edge_numbers++;
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->No, (int)edge_numbers });
				edge_list[(*itconflict)->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });	
			}
			else{
				mplAssert((*itconflict)->parent->Is_Parent);
				//avoid repeat (conflicts of stitch realation nodes may represent same conflict)
				edge_numbers++;	
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->parent->No, (int)edge_numbers });
				edge_list[(*itconflict)->parent->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });
			}
		}
		}
	std::cout<<"EDGE list generated with size: "<<edge_numbers<<std::endl;
	std::vector<uint32_t> starting_indexes;
	uint32_t starting_index = 0;

	//Stitch Generation: recording strating index of each parent node block
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it){
		uint32_t childs_num = ((*it)->Childs).size();
		mplAssert(childs_num > 1 || childs_num == 0);
		starting_indexes.push_back(starting_index);
		if(childs_num > 1){
			starting_index += 9*(childs_num-1);
		}
	}
	bool need_debug = false;
	//For debug
	if(comp_id == m_db->dbg_comp_id()){
		need_debug = true;
		for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
			if((*it)->No == 0){
				std::cout<<"originial node "<<(*it)->Stitch_No<<" no-stitch graph node "<<(*it)->parent->No<<std::endl;
			}
			else{
				std::cout<<"originial node "<<(*it)->Stitch_No<<" no-stitch graph node "<<(*it)->No<<std::endl;
			}
		}
	}
	//BUILD DL MATRIX WITHOUT STITCH
	//Node structure done,  then build DL system
	int row_numbers = vertex_numbers * m_db->color_num() + 1;
	int col_numbers = edge_numbers * m_db->color_num() + vertex_numbers;
	std::cout<<"DL1: row/col is"<<row_numbers<< " "<<col_numbers<<std::endl;
	std::cout<<"DL1: vertex/edge number is"<<vertex_numbers<< " "<<edge_numbers<<std::endl;
	DL_Init(dl, row_numbers, col_numbers);

	// Insert DL cells
	
	for (uint32_t it = 1; it < edge_list.size(); ++it)
	{
		for (int32_t i = 1; i <= m_db->color_num(); ++i)
		{
			//Insert elements representing nodes color ( first edge_num cols)
			Cell_Insert(dl,(it - 1)*m_db->color_num() + i,it);
			for (auto j = edge_list[it].begin(); j != edge_list[it].end(); ++j)
			{
				//Insert elements representing edge conflict
				Cell_Insert(dl,(it - 1)*m_db->color_num() + i,vertex_numbers + (j->No - 1)*m_db->color_num() + i);
			}
		}
	}
	for (uint32_t i = 0; i < edge_numbers; ++i)
	{
		for (int32_t j = 1; j <= m_db->color_num(); ++j)
		//Insert singelon row
			Cell_Insert(dl,vertex_numbers * m_db->color_num() + 1,vertex_numbers + i * m_db->color_num() + j);
	}


	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateDuplicateLGConflicts();
		if(!(*it)->Is_Parent)
		//else, add its parent node if it has not been added into node_wo_stitch_list
		{
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateDuplicateLGConflicts();
			}

		}
	}

	//Firstly solve dl without stitches
	std::vector<int> selected_rows_by_dl1;
	selected_rows_by_dl1 = core_solve_dl(dl,edge_list,row_numbers, col_numbers,(int)vertex_numbers,(int)m_db->color_num(),need_debug);
	//decode_mat for transforming the selected row to represented color. First: node id, Second: color;
	std::vector<std::vector<std::pair<uint32_t,uint32_t>>> decode_mat;
	//cost 1 is the cost of non-stitch solver
	decode_row_results(selected_rows_by_dl1, color_vector,(int)vertex_numbers, (int)m_db->color_num(), decode_mat,node_list);
	double cost1 = new_calc_cost(g,color_vector);
	std::cout<<"cost1 is "<<cost1<<std::endl;
	if(cost1 > 0){

		//if cost1 is larger than 1, which means that there is one more conflicts, than we insert stitch in the node
		//1 st dim: num_parent node
		//2 nd dim: num_stitch of this parent node
		//3 rd pair: deivided child node set by this stitch
		std::vector<std::vector<std::pair<std::set<uint32_t>,std::set<uint32_t>>>> encode_mat;
		for(auto i = 0; i<vertex_numbers; i++){
			std::vector<std::pair<std::set<uint32_t>,std::set<uint32_t>>> small_encode_mat;
			encode_mat.push_back(small_encode_mat);
		}

		//I traverse all of the stitch edges for getting the encode matrix
		boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
				//if two nodes are stitch relationships
				if(boost::get(boost::edge_weight, g, e12.first) < 0){
					std::set<uint32_t> set1;
					set1.insert(v1);
					std::set<uint32_t> set2;
					set2.insert(v2);
					std::pair<std::set<uint32_t>,std::set<uint32_t>> encode_mat_per_stitch = std::make_pair(set1,set2);
					push_adj_into_set(v1,set1,set2);
					push_adj_into_set(v2,set2,set1);

				}
			}
		}
		//Decode matrix generation
		count = -1;
		for(auto parent_no = 0; i<vertex_numbers; parent_no++){
			uint32_t childs_num = ((node_wo_stitch_list[parent_no])->Childs).size();
			count++;
			if(childs_num == 0){continue;}
			mplAssert(childs_num > 1);
			for(auto stitch_no = 0; stitch_no < childs_num - 1; stitch_no++){
				std::vector<std::pair<uint32_t,uint32_t>> row_decoder;
			}
		}
		//BUILD DL MATRIX WITH STITCH
		//Node structure done,  then build DL system
		 row_numbers = vertex_numbers * m_db->color_num() + 1 + ;
		 col_numbers = edge_numbers * m_db->color_num() + vertex_numbers;
		std::cout<<"DL1: row/col is"<<row_numbers<< " "<<col_numbers<<std::endl;
		std::cout<<"DL1: vertex/edge number is"<<vertex_numbers<< " "<<edge_numbers<<std::endl;
		DL_Init(dl, row_numbers, col_numbers);
	}
	return cost1;
}




double SimpleMPL::solve_by_dancing_link_GPU(SimpleMPL::graph_type& g,std::vector<int8_t>& color_vector){
	//Due to graph does not contain a parent&child node system
	// we simply redesign a node struct to achieve this kind of system
	mplAssert(num_vertices(g) == color_vector.size());
	DancingLink dl; 
	std::vector<Vertex*> node_list;
	node_list.resize(num_vertices(g));
	uint32_t vertex_numbers = 0;
	uint32_t edge_numbers = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
			std::cout<<v1<<" "<<v2<<" "<<boost::get(boost::edge_weight, g, e12.first)<<std::endl;
			//if node source and node target is not created. New one/
			if(!node_list[v1]){
				Vertex* source_vertex = new Vertex;
				mplAssert(source_vertex->Conflicts.empty());
				node_list[v1] = source_vertex;
				source_vertex->Stitch_No = v1;
				vertex_numbers ++;
			}
			if(!node_list[v2]){
				Vertex* target_vertex = new Vertex;
				target_vertex->Stitch_No = v2;
				node_list[v2] = target_vertex;
				vertex_numbers ++;
			}

			//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0){
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				if(node_list[v1]->parent != NULL && node_list[v1]->parent== node_list[v2]->parent) continue;
				Vertex* parent_vertex = new Vertex;
				mplAssert(parent_vertex->Childs.empty());
				parent_vertex->parentOf(node_list[v1]);
				parent_vertex->parentOf(node_list[v2]);
				vertex_numbers -= 1;
				continue;
			}
			node_list[v1]->Conflicts.insert(node_list[v2]);
		}
	}
	//need to store node list, which does not contain child node
	std::set<Vertex*> node_wo_stitch_list;
	int index = 1;
	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateConflicts();
		if((*it)->Is_Parent){
			node_wo_stitch_list.insert((*it));
			(*it)->No = index;
			index ++;
		}
		//else, add its parent node if it has not been added into node_wo_stitch_list
		else{
			if((*it)->parent->Is_Parent == false){
				std::cout<<(*it)->Stitch_No<<std::endl;
				std::cout<<(*it)->parent->Stitch_No<<std::endl;
			}
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateConflicts();
				(*it)->parent->No = index;
				index ++;
			}

		}
	}
	mplAssert(node_wo_stitch_list.size() == vertex_numbers);


	std::vector<std::list<Edge_Simple> >   edge_list;
	edge_list.resize(vertex_numbers + 1);
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it) {
		mplAssert((*it)->No != 0);
		for(std::set<Vertex*>::iterator itconflict = (*it)->Conflicts_in_LG.begin(); itconflict != (*it)->Conflicts_in_LG.end(); ++itconflict) {
			if((*itconflict)->Is_Parent){
				mplAssert((*itconflict)->No != 0);
				edge_numbers++;
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->No, (int)edge_numbers });
				edge_list[(*itconflict)->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });	
			}
			else{
				mplAssert((*itconflict)->parent->Is_Parent);
				//avoid repeat (conflicts of stitch realation nodes may represent same conflict)
				edge_numbers++;	
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->parent->No, (int)edge_numbers });
				edge_list[(*itconflict)->parent->No].push_back(Edge_Simple{ (*it)->No, (int)edge_numbers });
			}
		}
		}
	std::cout<<"EDGE list generated with size: "<<edge_numbers<<std::endl;
	std::vector<uint32_t> starting_indexes;
	uint32_t starting_index = 0;

	//Stitch Generation Step 1: recording strating index of each parent node block
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it){
		uint32_t childs_num = ((*it)->Childs).size();
		mplAssert(childs_num > 1 || childs_num == 0);
		starting_indexes.push_back(starting_index);
		if(childs_num > 1){
			starting_index += pow(m_db->color_num(),childs_num);
		}
	}

	//int row_numbers = vertex_numbers * m_db->color_num() + 1 + starting_index;
	int row_numbers = vertex_numbers * m_db->color_num() + starting_index;
	int col_numbers = edge_numbers * m_db->color_num() + vertex_numbers;
	std::cout<<"row/col is"<<row_numbers<< " "<<col_numbers<<std::endl;
	std::cout<<"vertex/edge number is"<<vertex_numbers<< " "<<edge_numbers<<std::endl;


	//Step 1 of Haoyu GPU version .initilize several matrix/arrays required by Haoyu
    int **dl_matrix;

    dl_matrix = new int *[row_numbers];
    for (int i = 0; i < row_numbers; i++)
    {
        dl_matrix[i] = new int[col_numbers];
		std::fill((dl_matrix[i]), (dl_matrix[i]) + col_numbers, 0);
    }

	int *deleted_cols = new int[col_numbers];
	int *col_group = new int[col_numbers];
	int *results = new int[row_numbers];

	//Step 2 of Haoyu GPU version . Assign values in col_group several
	for(uint32_t i = 0; i<(uint32_t)col_numbers; i++){
		if(i<(uint32_t)vertex_numbers){
			col_group[i] = -1 - i;
		}
		else{
			col_group[i] = (i - vertex_numbers)/m_db->color_num() + 1;
		}
		std::cout<<col_group[i]<<" ";
	}
	std::cout<<"END"<<std::endl;


	//Step 3 of Haoyu GPU version . Insert ones in dl matrix
	
	for (uint32_t it = 1; it < edge_list.size(); ++it)
	{
		for (int i = 1; i <= m_db->color_num(); ++i)
		{
			dl_matrix[(it - 1)*m_db->color_num() + i - 1][it -1 ] = 1;
			for (auto j = edge_list[it].begin(); j != edge_list[it].end(); ++j)
			{
				dl_matrix[(it - 1)*m_db->color_num() + i -1][vertex_numbers + (j->No - 1)*m_db->color_num() + i -1] = 1;
			}
		}
	}
	//ones representing singelon row
	// for (uint32_t i = 0; i < edge_numbers; ++i)
	// {
	// 	for (int j = 1; j <= m_db->color_num(); ++j)
	// 		dl_matrix[vertex_numbers * m_db->color_num() + 1 - 1][vertex_numbers + i * m_db->color_num() + j -1 ] = 1;
	// }


	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		(*it)->updateDuplicateLGConflicts();
		if(!(*it)->Is_Parent)
		//else, add its parent node if it has not been added into node_wo_stitch_list
		{
			mplAssert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateDuplicateLGConflicts();
			}

		}
	}

	uint32_t count = 0;
	for(std::set<Vertex*>::iterator v = node_wo_stitch_list.begin(); v != node_wo_stitch_list.end(); ++v){
		uint32_t childs_num = ((*v)->Childs).size();
		//std::cout<<count<<"th node with child num "<<childs_num<<std::endl;
		uint32_t child_count = 1;
		for(std::vector<Vertex*>::iterator child = ((*v)->Childs).begin(); child != ((*v)->Childs).end(); ++child){
			//c_conflit = child conflict
			for(std::set<Vertex*>::iterator c_conflit = ((*child)->Conflicts_in_LG).begin(); c_conflit != ((*child)->Conflicts_in_LG).end(); ++c_conflit){
				for(int i = 1;i<=pow(m_db->color_num(),childs_num);i++){
					int col_edge = -1;
					for(std::list<Edge_Simple>::iterator conflict_edge = edge_list[(*child)->parent->No].begin();conflict_edge != edge_list[(*child)->parent->No].end();++conflict_edge){
						if((*conflict_edge).target == (*c_conflit)->No){
							col_edge = (*conflict_edge).No;
							break;
						}
					}
					mplAssert(col_edge!=-1);
					//dl_matrix[vertex_numbers * m_db->color_num()  + starting_indexes[count] + i][vertex_numbers + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))] = 1;
					dl_matrix[vertex_numbers * m_db->color_num()  + starting_indexes[count] + i - 1][vertex_numbers + (col_edge -1)*(m_db->color_num()) + ((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))] = 1;
				} 
			}
			child_count ++;
		}
		count++;
	}

	//Stitch Generation Step 3: Decode matrix generation. pair: 1st = stitch No, 2nd = color  
	count = -1;
	std::vector<std::vector<std::pair<uint32_t,uint32_t>>> decode_mat;
	for(std::set<Vertex*>::iterator v = node_wo_stitch_list.begin(); v != node_wo_stitch_list.end(); ++v){
		uint32_t childs_num = ((*v)->Childs).size();
		count++;
		if(childs_num == 0){continue;}
		for(int i = 1;i<=pow(m_db->color_num(),childs_num);i++){
			//std::cout<<i<<" "<<(*v)->No<<" "<<childs_num<<" "<<pow(m_db->color_num(),childs_num)<<std::endl;
			std::vector<std::pair<uint32_t,uint32_t>> row_decoder;
			uint32_t child_count =1;
			for(std::vector<Vertex*>::iterator child = ((*v)->Childs).begin(); child != ((*v)->Childs).end(); ++child){
				//std::cout<<(*child)->Stitch_No<<" "<<(i/(int)pow(m_db->color_num(),childs_num-count)%(int)(m_db->color_num())) <<std::endl;
				row_decoder.push_back(std::make_pair((*child)->Stitch_No,((i-1)/(int)pow(m_db->color_num(),childs_num-child_count)%(int)(m_db->color_num()))));
				child_count ++;
			}
			//std::cout<<vertex_numbers * m_db->color_num() + 1 + starting_indexes[count] + i<<" "<<(*v)->No <<std::endl;
			//dl_matrix[vertex_numbers * m_db->color_num() + starting_indexes[count] + i][(*v)->No-1] = 1;
			dl_matrix[vertex_numbers * m_db->color_num() + starting_indexes[count] + i - 1][(*v)->No-1] = 1;
			decode_mat.push_back(row_decoder);
		}
	}

	//print dl_matrix for debug
#if RECORD > 1
	std::ofstream myfile;
	myfile.open ("matrix_record.txt", std::ofstream::app);
	myfile<<row_numbers<<" " <<col_numbers<<endl;
    for(int x=0;x<row_numbers;x++)  // loop 3 times for three lines
    {
        for(int y=0;y<col_numbers;y++)  // loop for the three elements on the line
        {
            myfile<<dl_matrix[x][y]<<" ";  // display the current element out of the array
        }
    myfile<<endl;  // when the inner loop is done, go to a new line
    }
	myfile<<"========================================================="<<endl;
	myfile.close();
#endif

	mc_solver(dl_matrix, results, deleted_cols, col_group, vertex_numbers, row_numbers, col_numbers);
	std::vector<int8_t> color_results_wo_stitch (vertex_numbers);
	color_results_wo_stitch.assign(vertex_numbers,-1);
	for (int i = 0; i <row_numbers ; i++)
	{
		if(results[i] != 0){
			// if(i+1 > (int)vertex_numbers * m_db->color_num() +1){
			// 	std::vector<std::pair<uint32_t,uint32_t>>& row_decoder = decode_mat[i+1- vertex_numbers * m_db->color_num() -2];
			if(i+1 > (int)vertex_numbers * m_db->color_num() ){
				std::vector<std::pair<uint32_t,uint32_t>>& row_decoder = decode_mat[i+1- vertex_numbers * m_db->color_num() -1];
				for(std::vector<std::pair<uint32_t,uint32_t>>::iterator it = row_decoder.begin(); it != row_decoder.end(); ++it) {
					color_vector[(*it).first] = (*it).second;
				}
			}
			else{
				int No = (i) / m_db->color_num() + 1;
				int mask = (i) % m_db->color_num();
				color_results_wo_stitch[No-1] = mask;
			}
		}
		// std::cout << *i << std::endl;

	}


	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		if(color_vector[(*it)->Stitch_No]!= -1){continue;}
		if((*it)->No == 0){
			mplAssert((*it)->Is_Parent == false);
			mplAssert((*it)->parent->Stitch_No == -1);
			mplAssert(color_results_wo_stitch[(*it)->parent->No -1]!= -1);
			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->parent->No -1];
		}
		else{
			mplAssert(color_results_wo_stitch[(*it)->No -1]!=-1);
			color_vector[(*it)->Stitch_No] = color_results_wo_stitch[(*it)->No -1];
		}
	}
	double cost = calc_cost(g,color_vector);
	std::cout<<"cost is "<<cost<<std::endl;
	return cost;
}


void SimpleMPL::iterative_mark(SimpleMPL::graph_type& g,std::vector<uint32_t>& parent_node_ids, SimpleMPL::vertex_descriptor& v1){
	boost::graph_traits<graph_type>::adjacency_iterator vi2, vie2,next2;
	boost::tie(vi2, vie2) = boost::adjacent_vertices(v1, g);
	for (next2 = vi2; vi2 != vie2; vi2 = next2)
	{	
		++next2; 
		vertex_descriptor v2 = *vi2;
		std::pair<edge_descriptor, bool> e12 = boost::edge(v1, v2, g);
		mplAssert(e12.second);
	//if two nodes are stitch relationships 
		if(boost::get(boost::edge_weight, g, e12.first) < 0){
			if(parent_node_ids[(uint32_t)v2] == -1){
				parent_node_ids[(uint32_t)v2] = parent_node_ids[(uint32_t)v1];
				iterative_mark(g,parent_node_ids,v2);
			}
			else{
			mplAssert(parent_node_ids[(uint32_t)v2] == parent_node_ids[(uint32_t)v1]);
		}
		}
	}
}
double SimpleMPL::new_calc_cost(SimpleMPL::graph_type& g,std::vector<int8_t>& vColor){
	double cost = 0;
	std::vector<uint32_t> parent_node_ids;
	parent_node_ids.assign(boost::num_vertices(g),-1);
	uint32_t parent_node_id = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
	for (boost::tie(vi1, vie1) = boost::vertices(g); vi1 != vie1; ++vi1)
	{	
		vertex_descriptor v1 = *vi1;
		if(parent_node_ids[(uint32_t)v1] == -1){
			parent_node_ids[(uint32_t)v1] = parent_node_id;
			iterative_mark(g,parent_node_ids,v1);
			parent_node_id++;
		}
	}


	//calculate conflict by parent node
	std::vector<std::vector<bool>> conflict_mat;
	for(uint32_t i = 0; i < parent_node_id; i++){
		std::vector<bool> row_mat;
		row_mat.assign(parent_node_id,false);
		conflict_mat.push_back(row_mat);
	}

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
		//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0){
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				cost += (vColor[v1] != vColor[v2])*m_db->parms.weight_stitch;
			}
			else{
				if(vColor[v1] == vColor[v2]){
					if(conflict_mat[parent_node_ids[(uint32_t)v1]][parent_node_ids[(uint32_t)v2]] == false){
						cost += 1;
						conflict_mat[parent_node_ids[(uint32_t)v1]][parent_node_ids[(uint32_t)v2]] = true;
						conflict_mat[parent_node_ids[(uint32_t)v2]][parent_node_ids[(uint32_t)v1]] = true;
					}
				}
			}
		}
	}






	if (num_vertices(g) == 21 && cost == 0.2){
		std::cout<<"BUG HERE!"<<std::endl;
	}
	return cost;
}
double SimpleMPL::calc_cost(SimpleMPL::graph_type& g,std::vector<int8_t> const& vColor)  
{
	mplAssert(vColor.size() == boost::num_vertices(g));
	double cost = 0;
	boost::graph_traits<graph_type>::vertex_iterator vi1, vie1;
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
		//if two nodes are stitch relationships and both of them are parent
			if(boost::get(boost::edge_weight, g, e12.first) < 0){
				//std::cout<<"Stitch:"<<v1<<" "<<v2<<std::endl;
				cost += (vColor[v1] != vColor[v2])*m_db->parms.weight_stitch;
			}
			else{
				cost += (vColor[v1] == vColor[v2]);
			}
		}
	}
	return cost;
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
		if (isVDDGND[v])
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
			for (std::vector<uint32_t>::const_iterator it = StitchRelation[v].begin(); it != StitchRelation[v].end(); ++it)
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
#ifdef _DGOUT
	return 1;
#endif
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
    // if (acc_obj_value != std::numeric_limits<uint32_t>::max())
    //    mplAssertMsg(acc_obj_value == component_conflict_num, "%u != %u", acc_obj_value, component_conflict_num);

	if (m_db->verbose())
		mplPrint(kDEBUG, "Component %u has %u patterns...%u conflicts\n", comp_id, (uint32_t)(itEnd-itBgn), component_conflict_num);

	return component_conflict_num;
}

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
	construct_component_graph(itBgn, pattern_cnt, dg, mGlobal2Local, vColor, vdd_set, flag);
#ifdef _OPENMP
#pragma omp critical(dgGlobal2Local)
#endif
	{
		std::map<uint32_t, uint32_t>().swap(dgGlobal2Local);
		dgGlobal2Local.insert(mGlobal2Local.begin(), mGlobal2Local.end());
	}


    // for debug, it does not affect normal run 
	// if (comp_id == m_db->dbg_comp_id())
    //     write_graph(dg, "graph_init");

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
    acc_obj_value = solve_graph_coloring(comp_id, dg, itBgn, pattern_cnt, simplify_strategy, vColor, vdd_set);

#ifdef _DGOUT
	return 1;
#endif

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
			std::set<uint32_t> parent_indexes;
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

void SimpleMPL::find_all_stitches(uint32_t vertex, std::vector<uint32_t>& stitch_vec){
	for(uint32_t stit = 0; stit < StitchRelation[vertex].size() ;stit++){
		if(std::find(stitch_vec.begin(), stitch_vec.end(), StitchRelation[vertex][stit]) == stitch_vec.end()){
			stitch_vec.push_back(StitchRelation[vertex][stit]);
			find_all_stitches(StitchRelation[vertex][stit],stitch_vec);
		}
	}
}
uint32_t SimpleMPL::conflict_num() 
{
	m_vConflict.clear();
	std::vector<rectangle_pointer_type> const& vPatternBbox = m_db->vPatternBbox;
	for (uint32_t v = 0; v != vPatternBbox.size(); ++v)
	{
		std::vector<uint32_t> stitch_vec;
		find_all_stitches(v,stitch_vec);
		int8_t color1 = vPatternBbox[v]->color();
		if (color1 >= 0 && color1 <= m_db->color_num())
		{
			std::set<uint32_t> parent_indexes;
			for (std::vector<uint32_t>::const_iterator itAdj = m_mAdjVertex[v].begin(); itAdj != m_mAdjVertex[v].end(); ++itAdj)
			{
				uint32_t u = *itAdj;
				// std::cout<<m_db->vParentPolygonId[u]<<std::endl;
				// if(parent_indexes.find(m_db->vParentPolygonId[u])!=parent_indexes.end()){continue;}
				// else{parent_indexes.insert(m_db->vParentPolygonId[u]);}
                if (v < u) // avoid duplicate 
                {
                    int8_t color2 = vPatternBbox[u]->color();
                    if (color2 >= 0 && color2 <= m_db->color_num())
                    {
                        if (color1 == color2) 
							{
								bool need_count = true;
								for(uint32_t j = 0; j < stitch_vec.size();j++){
									if(vPatternBbox[stitch_vec[j]]->color() == color1){
										if(stitch_vec[j] >=  v)	continue; //avoid duplicate, for each 
										need_count = false;
										break;
									}
								}
								if(need_count)	m_vConflict.push_back(std::make_pair(v, u));
							}
                            
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

void SimpleMPL::print_welcome() const
{
  mplPrint(kNONE, "\n\n");
  mplPrint(kNONE, "=======================================================================\n");
  mplPrint(kNONE, "                      OpenMPL - Version 1.1                          \n");
  mplPrint(kNONE, "                                by                                   \n");  
  mplPrint(kNONE, "                Yibo Lin, Bei Yu, Qi Sun and  David Z. Pan           \n");
  mplPrint(kNONE, "               ECE Department, University of Texas at Austin         \n");
  mplPrint(kNONE, "               CSE Department, Chinese University of Hong Kong       \n");
  mplPrint(kNONE, "                         Copyright (c) 2018                          \n");
  mplPrint(kNONE, "            Contact Authors:  {yibolin,dpan}@cerc.utexas.edu         \n");
  mplPrint(kNONE, "                              {byu, qsun}@cse.cuhk.edu.hk            \n");
  mplPrint(kNONE, "=======================================================================\n");
}

SIMPLEMPL_END_NAMESPACE

