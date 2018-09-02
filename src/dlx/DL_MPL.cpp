#include"DL_MPL.h"


std::vector<std::list<Edge_Simple>> Read_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers)
{
	std::ifstream filein(filename);
	// std::cout << "filename : " << filename << std::endl;
	filein >> vertex_numbers >> edge_numbers; // >> mask_numbers;
	/*
	std::cout << "vertex_numbers : " << vertex_numbers << std::endl;
	std::cout << "edge_numbers : " << edge_numbers << std::endl;
	std::cout << "mask_numbers : " << mask_numbers << std::endl;
	*/

	int source, target;
	std::vector<std::list<Edge_Simple>> edge_list;
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

std::vector<int> BFS_Order(std::vector<std::list<Edge_Simple>> & edge_list)
{
	std::vector<int> result_vector;
	std::queue<int> intermediate_queue;
	std::set<int> nonexistent;
	intermediate_queue.push(1);
	result_vector.push_back(1);
	int next;
	while (!intermediate_queue.empty())
	{
		next = intermediate_queue.front();
		intermediate_queue.pop();
		if (nonexistent.find(next) == nonexistent.end())
		{
			nonexistent.insert(next);
			result_vector.push_back(next);
			for (auto i = edge_list[next].begin(); i != edge_list[next].end(); i++)
				intermediate_queue.push(i->target);
			continue;
		}
		else
			continue;
	}
	return result_vector;
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
	int col_no;
	int max = MAX;
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

	std::vector<std::list<Edge_Simple>> edge_list = Read_Graph_File(infilename, vertex_numbers, edge_numbers); //, mask_numbers);

	if (whether_BFS)
		MPLD_search_vector = BFS_Order(edge_list);
	else
		MPLD_search_vector = Simple_Order(edge_list.size());
	row_numbers = vertex_numbers * mask_numbers + 1;
	col_numbers = edge_numbers * mask_numbers + vertex_numbers;
	fileout << row_numbers << std::endl;
	fileout << col_numbers << std::endl;
	int count = 0;
	for (auto it = edge_list.begin(); it != edge_list.end(); it++)
	{
		int temp = (it->size() + 1) * mask_numbers;
		count += temp;
	}
	count += edge_numbers*mask_numbers;

	fileout << count << std::endl;
	for (int it = 1; it < edge_list.size(); ++it)
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

bool Vertices_All_Covered(DancingLink & dl, int vertex_numbers)
{
	if (dl.DL_Header.Right->Col > vertex_numbers)
		return true;
	else
		return false;
}

void store_intermediate_process(DancingLink & dl, int this_col, std::set<int> & row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col)
{
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it) {
		Delete_the_Row_in_which_Col[*row_it] = this_col;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			Order_of_Row_Deleted_in_Col[row_ele->Col].push_back(*row_it);
	}
}

void recover_intermediate_process(DancingLink & dl, int this_col, std::set<int> row_set,
								std::vector<int> & Delete_the_Row_in_which_Col,
								std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col)
{
	for (auto row_it = row_set.begin(); row_it != row_set.end(); ++row_it)
	{
		// Delete_the_Row_in_which_Col[*row_it] = 0;
		for (auto row_ele = dl.Row_Header_Table[*row_it].Right; row_ele != &dl.Row_Header_Table[*row_it]; row_ele = row_ele->Right)
			Order_of_Row_Deleted_in_Col[row_ele->Col].pop_back();
	}
}

bool MPLD_X_Solver(DancingLink & dl, std::vector<int> & result_vec, std::set<std::pair<int, int>> & conflict_set, 
			int vertex_numbers, int mask_numbers,
			std::vector<int> & Delete_the_Row_in_which_Col,
			std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col, int depth, std::vector<int> & MPLD_search_vector, std::string result_file)
{
	// If there is no columns left or all the verteices are covered, then the algorithm terminates.
	if (dl.DL_Header.Right == &dl.DL_Header || Vertices_All_Covered(dl, vertex_numbers))
	{
		Decode(vertex_numbers, mask_numbers, result_vec, conflict_set, result_file);
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
	// std::cout << "depth : " << depth << std::endl;
	int this_col = Next_Column(MPLD_search_vector, depth);
	// std::cout << "this_col : " << this_col << std::endl;
	Cell *col = &dl.Col_Header_Table[this_col];
	LR_remove(*col);
	if (col->Children_Number == 0)
	{
		int last_row = Order_of_Row_Deleted_in_Col[this_col].back();
		int conflict_col = Delete_the_Row_in_which_Col[last_row];
		conflict_set.insert(std::make_pair(conflict_col, this_col));
		if (MPLD_X_Solver(dl, result_vec, conflict_set, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, MPLD_search_vector, result_file))
			return true;
	}
	for (Cell *j = col->Down; j != col; j = j->Down)
	{
		std::set<int> row_set;
		std::set<int> col_set;
		result_vec.push_back(j->Row);
		Select_All_Rows_Cols(dl, j->Row, row_set, col_set);
		store_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);

		Remove_Rows_Cols(dl, row_set, col_set);
		if (MPLD_X_Solver(dl, result_vec, conflict_set, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth + 1, MPLD_search_vector, result_file))
			return true;

		result_vec.pop_back();
		recover_intermediate_process(dl, this_col, row_set, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col);
		Recover_Rows_Cols(dl, row_set, col_set);
	}
	LR_recover(*col);
	return false;
}

void Decode(int vertex_numbers, int mask_numbers, std::vector<int> result_vec, std::set<std::pair<int, int>> conflict_set, std::string filename)
{
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

	if (!conflict_set.empty())
	{
		fileout << "Conflicts: " << std::endl;
		fileout << "No\t Source\t Target" << std::endl;
		int count = 1;
		for (auto i = conflict_set.begin(); i != conflict_set.end(); i++)
			fileout << count++ << "\t " << i->first << "\t " << i->second << std::endl;
		fileout << "=============================================" << std::endl;
	}
	
	fileout.close();
} 

void MPLD_Solver(std::string Graph_Filename, std::string Exact_Cover_Filename, bool whether_BFS, int mask_numbers, std::string result_file)
{
	int vertex_numbers;
	int edge_numbers;
	int row_numbers;
	int col_numbers;
	std::vector<int> result_vec;
	std::set<std::pair<int, int>> conflict_set;
	std::vector<int> Delete_the_Row_in_which_Col;
	std::vector<std::list<int>> Order_of_Row_Deleted_in_Col;
	std::vector<int> MPLD_search_vector;
	DancingLink dl;
	
	Convert_to_Exat_Cover(row_numbers, col_numbers, Graph_Filename, Exact_Cover_Filename, whether_BFS, vertex_numbers, edge_numbers, mask_numbers, MPLD_search_vector);
	Order_of_Row_Deleted_in_Col.resize(col_numbers + 1);
	Delete_the_Row_in_which_Col.resize(row_numbers + 1);
	DL_Load(dl, Exact_Cover_Filename);
	int depth = 1;
	MPLD_X_Solver(dl, result_vec, conflict_set, vertex_numbers, mask_numbers, Delete_the_Row_in_which_Col, Order_of_Row_Deleted_in_Col, depth, MPLD_search_vector, result_file);
}

