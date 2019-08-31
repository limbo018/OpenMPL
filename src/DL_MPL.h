#include"DL_struct.h"
#include "Msg.h"
#include<vector>
#include<queue>

/*
 * @brief: This function read the graph file and return the relevant information to help solving the MPLD problem.
 *
 * @param	filename:		the name of the file storing the graph 
 * @param	vertex_numbers:	number of vertices in the graph
 * @param	edge_numbers:	number of edges in the graph
 * @param	mask_numbers:	number of masks in the MPLD problem
 *
 * @return	std::map<int, std::list<Edge_Simple>>
 *			A map storing the graph. The index element of the map is source vertex, and the correspongding value 
 *			is a struct object storing the target vertex and edge No. For a undirected graph, each edge has two duplicates in the map,
 *			while their indices are source and target respectively.
 */
std::vector<std::list<Edge_Simple> > Read_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers);


/*
 * @brief: This function read the graph file which contains stitch candidate information and return the relevant information to help solving the MPLD problem.
 *
 * @param	filename:		the name of the file storing the graph 
 * @param	vertex_numbers:	number of vertices in the graph
 * @param	edge_numbers:	number of edges in the graph
 * @param	mask_numbers:	number of masks in the MPLD problem
 *
 * @return	std::map<int, std::list<Edge_Simple>>
 *			A map storing the graph. The index element of the map is source vertex, and the correspongding value 
 *			is a struct object storing the target vertex and edge No. For a undirected graph, each edge has two duplicates in the map,
 *			while their indices are source and target respectively.
 */
std::vector<std::list<Edge_Simple> > Read_Stitch_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers);




/*
 * @brief:	This function calculates and returns the maximal degree. 
 * @param	edge_list:		a map storing the graph information obtained from Read_Graph_File()
 * 
 * @return	int max_degree
 */

int  find_max_degree(std::vector<std::list<Edge_Simple> >  & edge_list);
/*
 * @brief:	This function calculates and returns the node of maximal degree. 
 * @param	edge_list:		a map storing the graph information obtained from Read_Graph_File()
 * 
 * @return	int node: the node with maximal degree 
 */
int find_max_degree_node(std::vector<std::list<Edge_Simple> >  & edge_list);



/*
 * @brief:	This function calculates and returns the degree_of_each_node
 * @param	edge_list:		a map storing the graph information obtained from Read_Graph_File()
 * @param node_degree: a vector indicating each node's degree
 * @return	int node: the node with minimal degree 
 */
int calcualte_degree_of_each_node(std::vector<std::list<Edge_Simple> >  & edge_list, std::vector<int> & node_degree);
/*
 * @brief:	This function calculates and returns the column covering order by BFS traversal. 
 * @param	edge_list:		a map storing the graph information obtained from Read_Graph_File()
 * 
 * @return	std::vector<int>
 *			a vector in which the elements are stored in order of BFS, each element means the column No. 
 */


std::vector<int> BFS_Order(std::vector<std::list<Edge_Simple> > & edge_list);

/*
 * @brief:	This function returns a simple order for next-column selecting in the order of vertex index.
 */
std::vector<int> Simple_Order(int size);

/*
 * @brief:	This function returns the column with least children.
 * @param	dl:		the dancinglink instance
 * 
 * @return	int
 *			the column No.
 */
int Sorting_Queue(DancingLink & dl);

/*
 * @brief:	This function consturct a exact cover problem according the graph and search rule.
 *			This function will call Read_Graph_File() to read the graph file. The exact cover 
 *			problem will be outputed into a file.
 * @param	infilename:		the name of the file storing the graph
 * @param	Exact_Cover_Filename:	the output file name
 * @param	whether_BFS:	whether adopt the BFS strategy
 * @param	vertex_numbers:	the number of vertices in the graph
 * @param	edge_numbers:	the number of edges in the graph
 * @param	mask_numbers:	the number of mask in the MPLD problem
 * 
 */
void Convert_to_Exat_Cover(int & row_numbers, int & col_numbers, std::string infilename, std::string Exact_Cover_Filename, bool whether_BFS,
	int & vertex_numbers, int & edge_numbers, int & mask_numbers, std::vector<int> & MPLD_search_vector);

/*
 * @brief:	Fetch the corresponding column from MPLD_search_vector
 * @param	
 * @return	int
 *			the column No.
 */
int Next_Column(std::vector<int> & MPLD_search_vector, int depth);
/*
 * @brief:	Fetch the corresponding column from MPLD_search_vector in stitch_enabled version
 * @param	
 * @return	int
 *			the column No.
 */
int Next_Column_stitch(DancingLink & dl,std::vector<int> & MPLD_search_vector, std::vector<bool>& col_cover_vector);
/*
 * @brief:	Whether all the corresponding columns of vertices in the graph are covered is one of the termination conditions.
 *			If all are covered, return true, else false.
 *			
 * @param	dl:				the dancinglink instance
 * @param	vertex_numbers:	number of vertices in the graph
 *
 */
bool Vertices_All_Covered(DancingLink & dl, int& vertex_numbers);

/*
 * @brief:	In order to find the conflict edges in the graph, we employ two intermediate maps to store the operation process.
 * 
 * @param	this_col:		the No. of the column we should choose in this step
 * @param	row_set:		a set containing all the rows we should cover in this step
 *
 * @param	Delete_the_Row_in_which_Col:	a map whose index is the No. of a row and the corresponding value is the No. of column,
 *											which means this row is covered due to the selection of the column.
 * @param	Order_of_Row_Deleted_in_Col:	a map whose index is the No. of a column, the corresponding value is a list storing 
 *											the children's row No. in the order of removement.
 */
void store_intermediate_process(DancingLink & dl, int this_col, std::set<int> & row_set,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col,std::vector<int>  &  conflict_col_table,std::vector<int>  &last_rows);
void efficient_store_intermediate_process(DancingLink & dl, int this_col, std::set<int> & row_set,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col);

/*
 * @brief:	This function is the reverse process of store_intermediate_process()
 */
void recover_intermediate_process(DancingLink & dl, int this_col, std::set<int>& row_set,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col,std::vector<int>  &  conflict_col_table,std::vector<int>  &last_rows);
void efficient_recover_intermediate_process(DancingLink & dl, int this_col, std::set<int>& row_set,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col);
/*
 * @brief:	The recursive implementation of Algortihm X*. 
 *			If all vertices in the graph are covered or the dlx has no column headers, then the function terminated.
 *			Whether this problem has a valid solution relies on the conflict_pair. If there are conflict edges in conflict_pair, 
 *			then no valid solution has been found, while the items in result_vec are still valid.
 * @param	dl:		the dancinglink instance
 * @param	result_vec:		a vector storing the valid coloring till now.
 * @param	conflict_pair:	a set storing the conflict edges till now.
 * @param	vertex_numbers:	the number of vertices in the graph.
 * @param	Delete_the_Row_in_which_Col:	a map whose index is the No. of a row and the corresponding value is the No. of column,
 *											which means this row is covered due to the selection of the column.
 * @param	Order_of_Row_Deleted_in_Col:	a map whose index is the No. of a column, the corresponding value is a list storing
 *											the children's row No. in the order of removement.
 * @param	result_file:	the file storing the final result
 */
bool MPLD_X_Solver(DancingLink & dl,std::vector<int8_t>& color_vector, std::vector<int> & result_vec, std::pair<int, int> & conflict_pair,
	int  vertex_numbers, int  mask_numbers,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col, 
	int depth, std::vector<int> & MPLD_search_vector,  const char*  result_file,std::vector<bool>& col_cover_vector,std::vector<int>& row_select_vector,
				std::vector<int> & partial_conflict_col_table,
			std::vector<int>  &  conflict_col_table,std::vector<int>  &partial_last_rows,std::vector<int>  &last_rows, bool & need_debug);

bool Efficient_MPLD_X_Solver(DancingLink & dl,std::vector<int8_t>& color_vector, std::vector<int> & result_vec, std::pair<int, int> & conflict_pair,
	int  vertex_numbers, int  mask_numbers,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int> > & Order_of_Row_Deleted_in_Col, 
	int depth, std::vector<int> & MPLD_search_vector,  const char*  result_file,std::vector<int>& row_select_vector,
				std::vector<int> & partial_conflict_col_table,
			std::vector<int>  &  conflict_col_table,bool & need_debug);


bool Efficient_MPLD_X_Solver_v2(DancingLink & dl, std::vector<int> & result_vec, std::pair<int, int>  & conflict_pair, 
			int vertex_numbers, std::vector<int> & Delete_the_Row_in_which_Col,
			std::vector<std::list<int> >  & Order_of_Row_Deleted_in_Col, int depth, std::vector<int> & MPLD_search_vector,
			std::vector<int> & partial_row_results, std::vector<int> & partial_col_results,std::vector<int> & col_results, bool & need_debug);

//The core function to solve one dl.
std::vector<int> core_solve_dl(DancingLink & dl, std::vector<std::list<Edge_Simple> > & edge_list,  int  row_numbers,  int  col_numbers,
 int  vertex_numbers, int mask_number,  bool & need_debug);

//The function to decode and calculate the cost of the results of DL.
void decode_row_results(std::vector<int> & final_result, std::vector<int8_t> & color_vector, int  vertex_number,
int  mask_number, std::vector<std::vector<std::pair<uint32_t,uint32_t>>>& decode_mat, std::vector<Vertex*> & node_list );

/*
 * @brief:	This function decode the result and convert the exact cover problem to origin coloring problem.
 * @param	vertex_numbers:		the number of vertices in the graph
 * @param	mask_numbers:		the number of masks in the MPLD problem
 * @param	result_vec:			a vector storing the all the chosen rows
 * @param	filename:			the file storing the result
 */
void Decode(int vertex_numbers, int mask_numbers, std::vector<int> result_vec, std::pair<int, int>  conflict_pair, std::string filename);


void Decode_OpenMPL(int vertex_numbers, int mask_numbers, std::vector<int8_t>& color_vector,std::vector<int> result_vec, std::pair<int, int> conflict_pair,  const char*  filename);

/*
 * @brief:	MPLD_Solver
 * 
 * @param	Graph_Filename:			the file storing the origin graph
 *									The data in it should be in the following format:
 *									- 1st line: vertex_numbers
 *									- 2nd line: edge_numbers
 *									- 3rd line: mask_numbers
 *									- following lines: source - target
 * @param	Exact_Cover_Filename:	the intermediate file storing the exact cover problem
 *									The data in it should be in the following format:
 *									- 1st line: row_numbers
 *									- 2nd line:	col_numbers
 *									- 3rd line: number of 1s in the dlx
 *									- following lines: the location of 1, [row_No, col_No]
 * @param	whether_BFS:			whether we use the BFS order
 * @param	mask_numbers:			the number of masks
 * @param	result_file:			the file storing the result
 */
void MPLD_Solver(std::string Graph_Filename, std::string Exact_Cover_Filename, bool whether_BFS, int mask_numbers, std::string result_file);