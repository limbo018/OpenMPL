/**
 * @file   DL_MPL.cpp
 * @author Wei Li 
 * @date   Oct 2019
 */
#include <algorithm>
#include <vector>
#include <boost/timer/timer.hpp>
#include"DL_MPL.h"

SIMPLEMPL_BEGIN_NAMESPACE

std::vector<std::list<Edge_Simple> >  Read_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers)
{
	std::ifstream filein(filename.c_str());
	// std::cout << "filename : " << filename << std::endl;
	filein >> vertex_numbers >> edge_numbers; // >> mask_numbers;
	/*
	std::cout << "vertex_numbers : " << vertex_numbers << std::endl;
	std::cout << "edge_numbers : " << edge_numbers << std::endl;
	std::cout << "mask_numbers : " << mask_numbers << std::endl;
	*/

	int source, target;
	std::vector<std::list<Edge_Simple> > edge_list;
	edge_list.resize(vertex_numbers + 1);
	for (int i = 1; i <= edge_numbers; i++)
	{
		filein >> source >> target;
		edge_list[source].push_back(Edge_Simple{ target, i });
		edge_list[target].push_back(Edge_Simple{ source, i });
	}
	filein.close();
	return edge_list;
}

std::vector<std::list<Edge_Simple> >  Read_Stitch_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers)
{
	std::ifstream filein(filename.c_str());
	int total_vertex_numbers;
	// std::cout << "filename : " << filename << std::endl;
	filein >> total_vertex_numbers; // >> mask_numbers;
	/*
	std::cout << "vertex_numbers : " << vertex_numbers << std::endl;
	std::cout << "edge_numbers : " << edge_numbers << std::endl;
	std::cout << "mask_numbers : " << mask_numbers << std::endl;
	*/


	int source, target;
	double weight;
	std::vector<Vertex*> node_list;
	node_list.resize(total_vertex_numbers);
	vertex_numbers = 0;
	edge_numbers = 0;
	while (filein >> source >> target >> weight){
		//if node source and node target is not created. New one/
		if(!node_list[source]){
			Vertex* source_vertex = new Vertex;
			assert(source_vertex->Conflicts.empty());
			node_list[source] = source_vertex;
			vertex_numbers ++;
		}
		if(!node_list[target]){
			Vertex* target_vertex = new Vertex;
			node_list[target] = target_vertex;
			vertex_numbers ++;
		}

		//if two nodes are stitch relationships and both of them are parent
		if(weight < 0){
			Vertex* parent_vertex = new Vertex;
			assert(parent_vertex->Childs.empty());
			parent_vertex->parentOf(node_list[source]);
			parent_vertex->parentOf(node_list[target]);
			vertex_numbers -= 1;
			continue;
		}
		node_list[source]->Conflicts.insert(node_list[target]);
	}
	filein.close();

	//need to store node list, which does not contain child node
	std::set<Vertex*> node_wo_stitch_list;
	int index = 1;
	for(std::vector<Vertex*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
		//if the node in parent node, means it has no stitch relations
		if((*it)->Is_Parent){
			node_wo_stitch_list.insert((*it));
			(*it)->updateConflicts();
			(*it)->No = index;
			index ++;
		}
		//else, add its parent node if it has not been added into node_wo_stitch_list
		else{
			assert((*it)->parent->Is_Parent);
			if(node_wo_stitch_list.find((*it)->parent) == node_wo_stitch_list.end()){
				node_wo_stitch_list.insert((*it)->parent);
				(*it)->parent->updateConflicts();
				(*it)->parent->No = index;
				index ++;
			}

		}
	}
    mplAssert(node_wo_stitch_list.size() == (unsigned int)vertex_numbers);

	std::vector<std::list<Edge_Simple> >  edge_list;
	edge_list.resize(vertex_numbers + 1);
	for(std::set<Vertex*>::iterator it = node_wo_stitch_list.begin(); it != node_wo_stitch_list.end(); ++it) {
		assert((*it)->No != 0);
		for(std::set<Vertex*>::iterator itconflict = (*it)->Conflicts_in_LG.begin(); itconflict != (*it)->Conflicts_in_LG.end(); ++itconflict) {
			if((*itconflict)->Is_Parent){
				assert((*itconflict)->No != 0);
				edge_numbers++;
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->No, edge_numbers });
				edge_list[(*itconflict)->No].push_back(Edge_Simple{ (*it)->No, edge_numbers });	
			}
			else{
				assert((*itconflict)->parent->Is_Parent);
				//avoid repeat (conflicts of stitch realation nodes may represent same conflict)
				edge_numbers++;	
				edge_list[(*it)->No].push_back(Edge_Simple{ (*itconflict)->parent->No, edge_numbers });
				edge_list[(*itconflict)->parent->No].push_back(Edge_Simple{ (*it)->No, edge_numbers });
			}
		}
		}
	std::cout<<"EDGE list generated with size: "<<edge_numbers<<std::endl;
	return edge_list;
}

std::vector<int> BFS_Order_no_stitch_first(std::vector<std::list<Edge_Simple> >  & edge_list,DancingLink & dl ){
	std::vector<int> result_vector;
	std::queue<int> intermediate_queue;
	std::set<int> nonexistent;
	std::vector<int> node_degree;
	node_degree.assign(edge_list.size(),0); //edge_list.size() == vertex_number + 1
	//int first_vertex = find_max_degree_node(edge_list);
	calcualte_degree_of_each_node(edge_list,node_degree);
	int max_degree = find_max_degree(edge_list);
	std::vector<int> rows_num;
	rows_num.assign(edge_list.size(),0); 
	int max_rows = 0;
	int min_rows = INT_MAX;
	for(auto i = 1; i<(int)edge_list.size();i++){
		rows_num[i] = dl.Col_Header_Table[i].Children_Number;
		if(max_rows < dl.Col_Header_Table[i].Children_Number){max_rows = dl.Col_Header_Table[i].Children_Number;}
		if(min_rows > dl.Col_Header_Table[i].Children_Number){min_rows = dl.Col_Header_Table[i].Children_Number;}
	}
	// std::cout<<"rows_num ";
	// for (auto i = rows_num.begin(); i != rows_num.end(); ++i)
	// 	std::cout << *i << ' ';
	// std::cout<<std::endl;
	//we calculate the case of each node as the root and return the smallest calculation costs
	int smallest_cost = INT_MAX;
	int smallest_cost_node = -1;
	std::vector<int> tmp_result_vector;
	for(auto root_node = 1; root_node < (int)edge_list.size(); root_node++){
		tmp_result_vector.clear();
		nonexistent.clear();
		intermediate_queue.push(root_node);
		tmp_result_vector.push_back(root_node);
		while (!intermediate_queue.empty())
		{
			int next = intermediate_queue.front();
			intermediate_queue.pop();
			if (nonexistent.find(next) == nonexistent.end())
			{
				nonexistent.insert(next);
				tmp_result_vector.push_back(next);
				// for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++)
				// 	intermediate_queue.push(i->target);
				//push the nodes to the queue by small-degree firstly order
				for(int row_num = min_rows; row_num <= max_rows; row_num++){
					for(int d = 1; d<= max_degree;d++){
						for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++){
							if(node_degree[i->target] == d && rows_num[i->target] == row_num){
								intermediate_queue.push(i->target);
							}
						}
					}
			}

				continue;
			}
			else
				continue;
		}
		//calculate the cost of each node : calculation function is \sum rows_num[node] * (node_num - index) (largest rows_num, we want to put it into the final location)
		int total_cost = 0;
		for(auto index = 1; index <edge_list.size();index++){
			total_cost += rows_num[tmp_result_vector[index]] * (edge_list.size() - index);
		}
		if(total_cost < smallest_cost){
			result_vector= tmp_result_vector;
			smallest_cost = total_cost;
			smallest_cost_node = root_node;
		}
		// std::cout<<"total_cost " <<total_cost<<std::endl;
		// std::cout<<"TMP RESULT IS"<<std::endl;
		// for (auto i = tmp_result_vector.begin(); i != tmp_result_vector.end(); ++i)
    	// 	std::cout << *i << ' ';
		// std::cout<<std::endl;
	}
	return result_vector;
}
std::vector<int> BFS_Order(std::vector<std::list<Edge_Simple> >  & edge_list)
{
	std::vector<int> result_vector;
	std::queue<int> intermediate_queue;
	std::set<int> nonexistent;
	std::vector<int> node_degree;
	node_degree.assign(edge_list.size(),0); //edge_list.size() == vertex_number + 1
	//int first_vertex = find_max_degree_node(edge_list);
	int first_vertex = calcualte_degree_of_each_node(edge_list,node_degree);
	intermediate_queue.push(first_vertex);
	result_vector.push_back(first_vertex);

	int max_degree = find_max_degree(edge_list);
	while (!intermediate_queue.empty())
	{
		int next = intermediate_queue.front();
		intermediate_queue.pop();
		if (nonexistent.find(next) == nonexistent.end())
		{
			nonexistent.insert(next);
			result_vector.push_back(next);
			// for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++)
			// 	intermediate_queue.push(i->target);
			//push the nodes to the queue by small-degree firstly order
			for(int d = 1; d<= max_degree;d++){
				for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++){
					if(node_degree[i->target] == d){
						intermediate_queue.push(i->target);
					}
				}
			}

			continue;
		}
		else
			continue;
	}
	return result_vector;
}

std::vector<int> BFS_Order_max_first(std::vector<std::list<Edge_Simple> >  & edge_list)
{
	std::vector<int> result_vector;
	std::queue<int> intermediate_queue;
	std::set<int> nonexistent;
	std::vector<int> node_degree;
	node_degree.assign(edge_list.size(),0); //edge_list.size() == vertex_number + 1
	int first_vertex = find_max_degree_node(edge_list);
	//int first_vertex = calcualte_degree_of_each_node(edge_list,node_degree);
	intermediate_queue.push(first_vertex);
	result_vector.push_back(first_vertex);
	//int first_vertex = find_max_degree_node(edge_list);
	calcualte_degree_of_each_node(edge_list,node_degree);
	int max_degree = find_max_degree(edge_list);
	while (!intermediate_queue.empty())
	{
		int next = intermediate_queue.front();
		intermediate_queue.pop();
		if (nonexistent.find(next) == nonexistent.end())
		{
			nonexistent.insert(next);
			result_vector.push_back(next);
			for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++)
				intermediate_queue.push(i->target);
			//push the nodes to the queue by max-degree firstly order


			// for(int d = max_degree; d>= 1;d--){
			// 	for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++){
			// 		if(node_degree[i->target] == d){
			// 			intermediate_queue.push(i->target);
			// 		}
			// 	}
			// }

			// for(int d = 1; d<= max_degree;d++){
			// 	for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++){
			// 		if(node_degree[i->target] == d){
			// 			intermediate_queue.push(i->target);
			// 		}
			// 	}
			// }
			continue;
		}
		else
			continue;
	}
	return result_vector;
}
int  find_max_degree_node(std::vector<std::list<Edge_Simple> >  & edge_list){
	int max_degree = 0;
	int vertex_node = 1;
	for ( int it = 1; it < edge_list.size(); ++it){
		if(edge_list[it].size()> max_degree){
			max_degree = edge_list[it].size();
			vertex_node = it;
		}
	}
	return vertex_node;
}

int  find_max_degree(std::vector<std::list<Edge_Simple> >  & edge_list){
	int max_degree = 0;
	int vertex_node = 1;
	for ( int it = 1; it < edge_list.size(); ++it){
		if(edge_list[it].size()> max_degree){
			max_degree = edge_list[it].size();
			vertex_node = it;
		}
	}
	return max_degree;
}

int calcualte_degree_of_each_node(std::vector<std::list<Edge_Simple> >  & edge_list, std::vector<int> & node_degree){
	int min_degree = std::numeric_limits<int>::max();
	int vertex_node = 1;
	for ( int it = 1; it < edge_list.size(); ++it){
		if(edge_list[it].size()< min_degree){
			min_degree = edge_list[it].size();
			vertex_node = it;
		}
		node_degree[it] = edge_list[it].size();
	}
	return vertex_node;
}
std::vector<int> Simple_Order(int size)
{
	std::vector<int> result_vector;
	result_vector.reserve(size);
	result_vector.push_back(1);
	for (int i = 1; i < size; i++)
	{
		result_vector.push_back(i);
	}
	return result_vector;
}

int Sorting_Queue(DancingLink & dl)
{
	int col_no = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::max();
	for (auto c = dl.DL_Header.Right; c != &dl.DL_Header; c = c->Right)
	{
		if (c->Children_Number < max)
		{
			max = c->Children_Number;
			col_no = c->Col;
		}
	}
	return col_no;
}

void Convert_to_Exat_Cover(int & row_numbers, int & col_numbers, std::string infilename, std::string Exact_Cover_Filename, bool whether_BFS,
	int & vertex_numbers, int & edge_numbers, int & mask_numbers, std::vector<int> & MPLD_search_vector)
{
	std::ofstream fileout(Exact_Cover_Filename);

	std::vector<std::list<Edge_Simple> >  edge_list = Read_Stitch_Graph_File(infilename, vertex_numbers, edge_numbers); //, mask_numbers);

	if (whether_BFS)
		MPLD_search_vector = BFS_Order(edge_list);
	else
		MPLD_search_vector = Simple_Order(edge_list.size());
	row_numbers = vertex_numbers * mask_numbers + 1;
	col_numbers = edge_numbers * mask_numbers + vertex_numbers;
	fileout << row_numbers << std::endl;
	fileout << col_numbers << std::endl;
	int count = 0;
	int total_edge = 0;
	for (auto it = edge_list.begin(); it != edge_list.end(); it++)
	{
		if(it == edge_list.begin()){continue;}
		std::cout<<it->size()<<std::endl;
		total_edge += it->size();
 		int temp = (it->size() + 1) * mask_numbers;
		count += temp;
	}
	std::cout<<total_edge<<std::endl;
	assert(total_edge == 2 * edge_numbers );
	count += edge_numbers*mask_numbers;

	fileout << count << std::endl;
	for (unsigned int it = 1; it < edge_list.size(); ++it)
	{
		for (int i = 1; i <= mask_numbers; ++i)
		{
			fileout << (it - 1)*mask_numbers + i << " " << it << std::endl;
			for (auto j = edge_list[it].begin(); j != edge_list[it].end(); ++j)
			{
				fileout << (it - 1)*mask_numbers + i << " " << vertex_numbers + (j->No - 1)*mask_numbers + i << std::endl;
			}
		}
	}
	for (int i = 0; i < edge_numbers; ++i)
	{
		for (int j = 1; j <= mask_numbers; ++j)
			fileout << vertex_numbers * mask_numbers + 1 << " " << vertex_numbers + i * mask_numbers + j << std::endl;
	}
	fileout.close();
	return;
}

int Next_Column(std::vector<int> & MPLD_search_vector, int depth)
{
	return MPLD_search_vector[depth];
}

int Next_Column_stitch(DancingLink & dl,std::vector<int> & MPLD_search_vector)
{
	bool not_find_last_uncovered_col = true;
	int last_col = 0;
	for(unsigned int i = 0; i< MPLD_search_vector.size(); i++){
		Cell *col = &dl.Col_Header_Table[MPLD_search_vector[i]];
		//if the column(vertex) has been removed(selected)
		//if(col_cover_vector[MPLD_search_vector[i]])	continue;
		if(col->Left->Right != col) continue;
		else{
			mplAssert(col->InDLX);
			if(not_find_last_uncovered_col){
				last_col = MPLD_search_vector[i];
				not_find_last_uncovered_col = false;
			}
		}
		
		if(col-> Children_Number == 1){
			return MPLD_search_vector[i];
		}
	}
	assert(not_find_last_uncovered_col == false);
	return last_col;
}
bool Vertices_All_Covered(DancingLink & dl, int& vertex_numbers)
{
	if (dl.DL_Header.Right->Col > vertex_numbers)
		return true;
	else
		return false;
}

void store_intermediate_process(DancingLink & dl, int this_col, std::set<int> & row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col,std::vector<int>  &  conflict_col_table,std::vector<int>  &last_rows)
{
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it) {
		Delete_the_Row_in_which_Col[*row_it] = this_col;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			{	//std::cout<<"row deleted is "<<*row_it<<"in column "<<this_col<<std::endl;
							Order_of_Row_Deleted_in_Col[row_ele->Col].push_back(*row_it);
			conflict_col_table[row_ele->Col] = this_col;
			last_rows[row_ele->Col] = (*row_it);}
	}
}

void efficient_store_intermediate_process(DancingLink & dl, int this_col, std::set<int> & row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col)
{
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it) {
		Delete_the_Row_in_which_Col[*row_it] = this_col;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			{	//std::cout<<"row deleted is "<<*row_it<<"in column "<<this_col<<std::endl;
							Order_of_Row_Deleted_in_Col[row_ele->Col].push_back(*row_it);}
	}
}

void recover_intermediate_process(DancingLink & dl, int this_col, std::set<int>& row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col,std::vector<int>  &  conflict_col_table,std::vector<int>  &last_rows)
{
	(void)this_col;
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it)
	{
		 Delete_the_Row_in_which_Col[*row_it] = 0;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			{Order_of_Row_Deleted_in_Col[row_ele->Col].pop_back();
			conflict_col_table[row_ele->Col] = Order_of_Row_Deleted_in_Col[row_ele->Col].back();
			last_rows[row_ele->Col] =  0;}
	}
}

void efficient_recover_intermediate_process(DancingLink & dl, int this_col, std::set<int>& row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col)
{
	(void)this_col;
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it)
	{
		 Delete_the_Row_in_which_Col[*row_it] = 0;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			{Order_of_Row_Deleted_in_Col[row_ele->Col].pop_back();}
	}
}

bool MPLD_X_Solver(DancingLink & dl, std::vector<int8_t>& color_vector,std::vector<int> & result_vec, std::pair<int, int>  & conflict_pair, 
			int vertex_numbers, int mask_numbers,
			std::vector<int> & Delete_the_Row_in_which_Col,
			std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col, int depth, std::vector<int> & MPLD_search_vector, const char* result_file,
			std::vector<bool>& col_cover_vector,std::vector<int>& row_select_vector,
			std::vector<int> & partial_conflict_col_table,
			std::vector<int>  &  conflict_col_table,std::vector<int>  &partial_last_rows,std::vector<int>  &last_rows)
{
	// If there is no columns left or all the verteices are covered, then the algorithm terminates.
	if (dl.DL_Header.Right == &dl.DL_Header || Vertices_All_Covered(dl, vertex_numbers))
	{
		//Decode_OpenMPL(vertex_numbers, mask_numbers,color_vector,result_vec, conflict_pair, result_file);
		return true;
	}
	/* 
	if (dl.DL_Header.Right->Col > vertex_numbers )
	{
		this_col = Sorting_Queue(dl);
	}
	else
	{
	*/
	Cell *col;
	//TODO: this depth should be a bug cause sometimes it is not by order (need a bool vector to record the col cover information)
	int this_col = Next_Column_stitch(dl,MPLD_search_vector);
	
	col = &dl.Col_Header_Table[this_col];
	LR_remove(*col);
	col_cover_vector[this_col] = true;
	if (col->Children_Number == 0)
	{
		// int last_row = Order_of_Row_Deleted_in_Col[this_col].back();
		// int conflict_col = Delete_the_Row_in_which_Col[last_row];
		//std::cout<<"last row is "<<last_row<<std::endl;
		//mplAssert(conflict_col == conflict_col_table[this_col]);
		conflict_pair = std::make_pair(conflict_col_table[this_col], this_col);
		//Also,  record the partial results (rows) in the last detected conflict
		row_select_vector.assign(result_vec.begin(),result_vec.end());
		row_select_vector.push_back(last_rows[this_col]);
		//We should should push back the last row(coloring result of this conflict col)
		partial_last_rows.assign(last_rows.begin(), last_rows.end());
		partial_conflict_col_table.assign(conflict_col_table.begin(), conflict_col_table.end());
		// if (MPLD_X_Solver(dl, color_vector,result_vec, conflict_pair, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, MPLD_search_vector, result_file))
		// 	return true;
	}
	for (Cell *j = col->Down; j != col; j = j->Down)
	{
		std::set<int> row_set;
		std::set<int> col_set;
		result_vec.push_back(j->Row);
		Select_All_Rows_Cols(dl, j->Row, row_set, col_set);
		store_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col,conflict_col_table,last_rows);

		Remove_Rows_Cols(dl, row_set, col_set);
		if (MPLD_X_Solver(dl, color_vector,result_vec, conflict_pair, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, 
		MPLD_search_vector, result_file,col_cover_vector,row_select_vector,partial_conflict_col_table,conflict_col_table,partial_last_rows,last_rows))
			return true;

		result_vec.pop_back();
		recover_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col,conflict_col_table,last_rows);
		Recover_Rows_Cols(dl, row_set, col_set);
	}
	LR_recover(*col);
	col_cover_vector[this_col] = false;
	
	return false;
}

bool Efficient_MPLD_X_Solver(DancingLink & dl,std::vector<int8_t>& color_vector, std::vector<int> & result_vec, std::pair<int, int>  & conflict_pair, 
			int vertex_numbers, int mask_numbers,
			std::vector<int> & Delete_the_Row_in_which_Col,
			std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col, int depth, std::vector<int> & MPLD_search_vector, const char* result_file,
			std::vector<int> & partial_row_results, std::vector<int> & partial_col_results,std::vector<int> & col_results)
{
	if (dl.DL_Header.Right == &dl.DL_Header || Vertices_All_Covered(dl, vertex_numbers))
	{
		//Decode_OpenMPL(vertex_numbers, mask_numbers,color_vector,result_vec, conflict_pair, result_file);
		return true;
	}

	Cell *col;
	int this_col = Next_Column_stitch(dl,MPLD_search_vector);
#ifdef DEBUG_DANCINGLINKCOLORING
    mplPrint(kDEBUG, "BEGIN: this_col : %d depth : %d\n", this_col, depth);
#endif

	col = &dl.Col_Header_Table[this_col];
	LR_remove(*col);
	col_results.push_back(this_col);
	if (col->Children_Number == 0)
	{
		int last_row = Order_of_Row_Deleted_in_Col[this_col].back();
		int conflict_col = Delete_the_Row_in_which_Col[last_row];
		conflict_pair = std::make_pair(conflict_col, this_col);
		partial_row_results.assign(result_vec.begin(),result_vec.end());
		partial_row_results.push_back(last_row);
		partial_col_results.assign(col_results.begin(), col_results.end());
	}
	
	for (Cell *j = col->Down; j != col; j = j->Down)
	{
		std::set<int> row_set;
		std::set<int> col_set;
		result_vec.push_back(j->Row);
		Select_All_Rows_Cols(dl, j->Row, row_set, col_set);
		efficient_store_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
		Remove_Rows_Cols(dl, row_set, col_set);
		if (Efficient_MPLD_X_Solver(dl,color_vector, result_vec, conflict_pair, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, 
		MPLD_search_vector, result_file,partial_row_results,partial_col_results,col_results))
			return true;

		result_vec.pop_back();
		efficient_recover_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
		Recover_Rows_Cols(dl, row_set, col_set);
	}
	LR_recover(*col);
#ifdef DEBUG_DANCINGLINKCOLORING
    mplPrint(kDEBUG, "END: this_col : %d depth : %d\n", this_col, depth);
#endif
	col_results.pop_back();
	return false;
}

bool Efficient_MPLD_X_Solver_v2(DancingLink & dl, std::vector<int> & result_vec, std::pair<int, int>  & conflict_pair, 
			int vertex_numbers, std::vector<int> & Delete_the_Row_in_which_Col,
			std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col, int depth, std::vector<int> & MPLD_search_vector,
			std::vector<int> & partial_row_results, std::vector<int> & partial_col_results,std::vector<int> & col_results,
			std::vector<std::vector<int>> &early_quit_count,bool & early_quit)
{
	if (dl.DL_Header.Right == &dl.DL_Header || Vertices_All_Covered(dl, vertex_numbers))
	{
		return true;
	}

	Cell *col;
	int this_col = Next_Column_stitch(dl,MPLD_search_vector);

	col = &dl.Col_Header_Table[this_col];
	LR_remove(*col);
	col_results.push_back(this_col);
	if (col->Children_Number == 0)
	{
		int last_row = Order_of_Row_Deleted_in_Col[this_col].back();
		int conflict_col = Delete_the_Row_in_which_Col[last_row];
		conflict_pair = std::make_pair(conflict_col, this_col);
		partial_row_results.assign(result_vec.begin(),result_vec.end());
		partial_row_results.push_back(last_row);
		partial_col_results.assign(col_results.begin(), col_results.end());
		early_quit_count[conflict_col][this_col] ++;
		early_quit_count[this_col][conflict_col] ++;
		if(early_quit_count[this_col][conflict_col] > 500 || early_quit_count[conflict_col][this_col]> 500){
			early_quit = true;
			LR_recover(*col);
			col_results.pop_back();
			return false;
		}
	}
	
	for (Cell *j = col->Down; j != col; j = j->Down)
	{
		std::set<int> row_set;
		std::set<int> col_set;
		result_vec.push_back(j->Row);
		Select_All_Rows_Cols(dl, j->Row, row_set, col_set);
		efficient_store_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
		Remove_Rows_Cols(dl, row_set, col_set);
		if (Efficient_MPLD_X_Solver_v2(dl, result_vec, conflict_pair, vertex_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, 
		MPLD_search_vector,partial_row_results,partial_col_results,col_results,early_quit_count,early_quit))
			return true;
		else{
			if(early_quit){
				result_vec.pop_back();
				efficient_recover_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
				Recover_Rows_Cols(dl, row_set, col_set);
				LR_recover(*col);
				col_results.pop_back();
#ifdef DEBUG_DANCINGLINKCOLORING
                mplPrint(kDEBUG, "END by early stop: this_col : %d depth : %d\n", this_col, depth);
#endif
				return false;
			}	
		}

		result_vec.pop_back();
		efficient_recover_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
		Recover_Rows_Cols(dl, row_set, col_set);
	}
	LR_recover(*col);
	col_results.pop_back();
	return false;
}


/***
 * The core function to solve one dl.
 * Typically, this function may run MPLD_X_Solver  multiple times due to the existence of some exact conflicts
 * Return : selected row 
 * **/
std::vector<int> core_solve_dl(DancingLink & dl, std::vector<std::list<Edge_Simple> > & edge_list,  int  row_numbers,  int  col_numbers,
 int  vertex_numbers, int mask_number){
	//the total conflict pairs 
	// boost::timer::cpu_timer dancing_link_timer;
	// //dancing_link_timer.start();
	std::pair<int, int>  conflict_pair;
	//the conflict pair if MPLD_X_Solver has NO a conflict-free solution
	std::set<std::pair<int, int> >  final_conflict;
	//selected result rows if MPLD_X_Solver has a conflict-free solution
	std::vector<int> selected_rows;
	//the final selected cols;
	std::vector<int> selected_cols;
	//the partial selected rows(coloring results actually) if MPLD_X_Solver has NO a conflict-free solution
	std::vector<int> partial_selected_rows;
	//the partial selected cols(selected )  if MPLD_X_Solver has NO a conflict-free solution
	std::vector<int> partial_selected_cols;
	//the final selected rows
	std::vector<int> final_result;
	//middle information for recover partial results
	std::vector<int> Delete_the_Row_in_which_Col;
	std::vector<std::list<int> >  Order_of_Row_Deleted_in_Col;
	Order_of_Row_Deleted_in_Col.resize(col_numbers + 1);
	Delete_the_Row_in_which_Col.resize(row_numbers + 1);
	std::vector<int> MPLD_search_vector;
	//MPLD_search_vector = BFS_Order(edge_list);
	//MPLD_search_vector = Simple_Order(edge_list.size());
	MPLD_search_vector = BFS_Order_no_stitch_first(edge_list,dl);
	//MPLD_search_vector = BFS_Order_max_first(edge_list);
	// for(auto i = 0; i<MPLD_search_vector.size(); i++){
	// 	std::cout<<MPLD_search_vector[i]<<" ";
	// }
	// std::cout<<std::endl;
	int depth = 1;
	std::vector<std::vector<int>> early_quit_count;
	for(auto i = 0; i <= vertex_numbers; i++){
		std::vector<int> row_early_quite_count;
		row_early_quite_count.assign(vertex_numbers+1,0);
		early_quit_count.push_back(row_early_quite_count);
	}
	bool early_quit = false;
	bool result = Efficient_MPLD_X_Solver_v2(dl, selected_rows, conflict_pair, vertex_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, 
		MPLD_search_vector,partial_selected_rows,partial_selected_cols,selected_cols,early_quit_count, early_quit);
	int iteration = 1;
	//std::cout<< "IN MPL:::::::::::::::::::::first try solve " << dancing_link_timer.format(6)<<std::endl;
	//dancing_link_timer.start();
	while(result==false && partial_selected_rows.size()< vertex_numbers){
		mplAssert(partial_selected_rows.size() == partial_selected_cols.size());
		iteration++;
#ifdef DEBUG_DANCINGLINKCOLORING
        mplPrint(kDEBUG, "iteration %d\n", iteration);
#endif
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
		mplAssert(selected_rows.size() == 0);
		for(int i = 1;i<=mask_number;i++){
			int this_col = vertex_numbers + (col_edge -1)*(mask_number) + i;
			//Cell *col = &dl.Col_Header_Table[this_col];
			Remove_Single_Col(dl,this_col);
		}
		//std::cout<< "IN MPL:::::::::::::::::::::remove conflict cols " << dancing_link_timer.format(6)<<std::endl;
		//dancing_link_timer.start();
		final_result.insert( final_result.end(), partial_selected_rows.begin(), partial_selected_rows.end() );
		//NOTE that the order between the two parts (remove edges and recover partial results are important! Cause some rows won't be removed this time due to the exact conflict edge is removed)
		for(int i = 0; i < partial_selected_rows.size(); i++){
			//if the row is selected ( in the partial result), recover the partial results
				std::set<int> row_set;
				std::set<int> col_set;
				Select_All_Rows_Cols(dl, partial_selected_rows[i], row_set, col_set);
				efficient_store_intermediate_process(dl, partial_selected_cols[i], row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
				Remove_Rows_Cols(dl, row_set, col_set);
		}
		//std::cout<< "IN MPL:::::::::::::::::::::recover partial results " << dancing_link_timer.format(6)<<std::endl;
		//dancing_link_timer.start();
		partial_selected_rows.clear();
		partial_selected_cols.swap(selected_cols);
		selected_cols.clear();
		depth = 1;
		for (auto &v: early_quit_count) {
			std::fill(v.begin(), v.end(), 0);
		}
		early_quit = false;
		result = Efficient_MPLD_X_Solver_v2(dl, selected_rows, conflict_pair, vertex_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, 
		MPLD_search_vector,partial_selected_rows,partial_selected_cols,selected_cols,early_quit_count, early_quit);
		//std::cout<< "IN MPL:::::::::::::::::::::second try solve " << dancing_link_timer.format(6)<<std::endl;
		//dancing_link_timer.start();
	}
	if(selected_rows.size() != 0){
		final_result.insert( final_result.end(), selected_rows.begin(), selected_rows.end() );
		mplAssert(final_result.size() == vertex_numbers);
	}
	else final_result.insert( final_result.end(), partial_selected_rows.begin(), partial_selected_rows.end() );
	return final_result;
}
/***
 * The function to decode and calculate the cost of the results of DL.
 * Args: selected rows
 * Return : selected row 
 * **/

void decode_row_results(std::vector<int> & final_result, std::vector<int8_t>& color_vector, int vertex_number,
int  mask_number, std::vector<std::vector<std::pair<uint32_t,uint32_t>>>& decode_mat, std::vector<Vertex*> & node_list ){
	std::vector<int8_t> color_results_wo_stitch (vertex_number);
	for (auto i = final_result.begin(); i != final_result.end(); i++)
	{
		// if the selected row is in the range of parent node (no stitch)
		// std::cout << *i << std::endl;
		if((*i) < vertex_number * mask_number + 1){
		int No = ((*i) - 1) / mask_number + 1;
		int mask = (*i + 2) % mask_number;
		color_results_wo_stitch[No-1] = mask; }
		// std::cout<<"selected rows: "<<(*i)<<std::endl;
		//if the selected_rows represent some stitches rows 
		if((*i) > (vertex_number * (uint32_t)(mask_number) + (uint32_t)1)){
			std::vector<std::pair<uint32_t,uint32_t>>& row_decoder = decode_mat[(*i)- vertex_number * mask_number -2];
			for(std::vector<std::pair<uint32_t,uint32_t>>::iterator it = row_decoder.begin(); it != row_decoder.end(); ++it) {
				color_vector[(*it).first] = (*it).second;
				// std::cout<<"Corresponding color "<<(*it).first<<" "<<(*it).second<<std::endl;
			}
		}
	}
	// std::cout<<"Finished"<<std::endl;

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
	return;
}

void decode_row_results_wo_skeleton(std::vector<int> & final_result, std::vector<int8_t>& color_vector, int vertex_number,
int  mask_number, std::vector<std::vector<std::pair<uint32_t,uint32_t>>>& decode_mat, std::vector<Vertex*> & node_list ){
	std::vector<int8_t> color_results_wo_stitch (vertex_number);
	for (auto i = final_result.begin(); i != final_result.end(); i++)
	{
		// if the selected row is in the range of parent node (no stitch)
		// std::cout << *i << std::endl;
		if((*i) < vertex_number * mask_number + 1){
		int No = ((*i) - 1) / mask_number + 1;
		int mask = (*i + 2) % mask_number;
		color_results_wo_stitch[No-1] = mask; }

		//if the selected_rows represent some stitches rows 
		if((*i) > (vertex_number * (uint32_t)(mask_number))){
			std::vector<std::pair<uint32_t,uint32_t>>& row_decoder = decode_mat[(*i)- vertex_number * mask_number -1];
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
	return;
}
void Decode(int vertex_numbers, int mask_numbers, std::vector<int> result_vec, std::pair<int, int> conflict_pair, std::string filename)
{
	(void)conflict_pair;
	(void)filename;
	(void)vertex_numbers;
	std::map<int, int> Final_Color;
	std::ofstream fileout(filename);
	fileout << "Solution : " << std::endl;
	for (auto i = result_vec.begin(); i != result_vec.end(); i++)
	{
		// std::cout << *i << std::endl;
		int No = ((*i) - 1) / mask_numbers + 1;
		int mask = (*i + 2) % mask_numbers + 1;
		Final_Color.insert(std::make_pair(No, mask));
	}
	for (auto i = Final_Color.begin(); i != Final_Color.end(); i++)
	{
		fileout << "Vertex " << i->first << " \t : " << i->second << std::endl;
	}
	fileout << "=============================================" << std::endl;

	// if (conflict_pair!=nullptr)
	// {
	// 	fileout << "Conflicts: " << std::endl;
	// 	fileout << "No\t Source\t Target" << std::endl;
	// 	fileout << count++ << "\t " << conflict_pair->first << "\t " <<conflict_pair->second << std::endl;
	// 	fileout << "=============================================" << std::endl;
	// }
	
	fileout.close();
} 

void Decode_OpenMPL(int vertex_numbers, int mask_numbers, std::vector<int8_t>& color_vector,std::vector<int> result_vec, std::pair<int, int> conflict_pair,  char* filename)
{
	(void)filename;
	(void)conflict_pair;
	for (auto i = result_vec.begin(); i != result_vec.end(); i++)
	{
		//std::cout << *i << std::endl;
		if((*i) < vertex_numbers * mask_numbers + 1){
		int No = ((*i) - 1) / mask_numbers + 1;
		int mask = (*i + 2) % mask_numbers;
		color_vector[No-1] = mask; }
	}

	// if (!conflict_pair.empty())
	// {
	// 	std::cout << "Conflicts: " << std::endl;
	// 	std::cout << "No\t Source\t Target" << std::endl;
	// 	int count = 1;
	// 	for (auto i = conflict_pair.begin(); i != conflict_pair.end(); i++)
	// 		{color_vector[i->second-1] = color_vector[i->first-1];
	// 		std::cout << count++ << "\t " << i->first << "\t " << i->second << std::endl;}
	// 	std::cout << "=============================================" << std::endl;
	// }
} 

void MPLD_Solver(std::string Graph_Filename, std::string Exact_Cover_Filename, bool whether_BFS, int mask_numbers,  const char* result_file)
{
	(void)result_file;
	int vertex_numbers;
	int edge_numbers;
	int row_numbers;
	int col_numbers;
	std::vector<int8_t> color_vec;
	std::vector<int> result_vec;
	std::pair<int, int>  conflict_pair;
	std::vector<int> Delete_the_Row_in_which_Col;
	std::vector<std::list<int> >  Order_of_Row_Deleted_in_Col;
	std::vector<int> MPLD_search_vector;
	DancingLink dl;
	std::cout<<"Ready to generate cover matrix"<<std::endl;
	Convert_to_Exat_Cover(row_numbers, col_numbers, Graph_Filename, Exact_Cover_Filename, whether_BFS, vertex_numbers, edge_numbers, mask_numbers, MPLD_search_vector);
	Order_of_Row_Deleted_in_Col.resize(col_numbers + 1);
	Delete_the_Row_in_which_Col.resize(row_numbers + 1);
	DL_Load(dl, Exact_Cover_Filename);
	//int depth = 1;
	//MPLD_X_Solver(dl, color_vec,result_vec, conflict_pair, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, MPLD_search_vector, result_file,col_cover_vector);
}

SIMPLEMPL_END_NAMESPACE
