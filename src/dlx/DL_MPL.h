#include"DL_struct.h"
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
std::vector<std::list<Edge_Simple>> Read_Graph_File(std::string filename, int & vertex_numbers, int & edge_numbers);

/*
 * @brief:	This function calculates and returns the column covering order by BFS traversal. 
 * @param	edge_list:		a map storing the graph information obtained from Read_Graph_File()
 * 
 * @return	std::vector<int>
 *			a vector in which the elements are stored in order of BFS, each element means the column No. 
 */

std::vector<int> BFS_Order(std::vector<std::list<Edge_Simple>> & edge_list);

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
 * @brief:	Whether all the corresponding columns of vertices in the graph are covered is one of the termination conditions.
 *			If all are covered, return true, else false.
 *			
 * @param	dl:				the dancinglink instance
 * @param	vertex_numbers:	number of vertices in the graph
 *
 */
bool Vertices_All_Covered(DancingLink & dl, int vertex_numbers);

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
	std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col);

/*
 * @brief:	This function is the reverse process of store_intermediate_process()
 */
void recover_intermediate_process(DancingLink & dl, int this_col, std::set<int> row_set,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col);

/*
 * @brief:	The recursive implementation of Algortihm X*. 
 *			If all vertices in the graph are covered or the dlx has no column headers, then the function terminated.
 *			Whether this problem has a valid solution relies on the conflict_set. If there are conflict edges in conflict_set, 
 *			then no valid solution has been found, while the items in result_vec are still valid.
 * @param	dl:		the dancinglink instance
 * @param	result_vec:		a vector storing the valid coloring till now.
 * @param	conflict_set:	a set storing the conflict edges till now.
 * @param	vertex_numbers:	the number of vertices in the graph.
 * @param	Delete_the_Row_in_which_Col:	a map whose index is the No. of a row and the corresponding value is the No. of column,
 *											which means this row is covered due to the selection of the column.
 * @param	Order_of_Row_Deleted_in_Col:	a map whose index is the No. of a column, the corresponding value is a list storing
 *											the children's row No. in the order of removement.
 * @param	result_file:	the file storing the final result
 */
bool MPLD_X_Solver(DancingLink & dl, std::vector<int> & result_vec, std::set<std::pair<int, int>> & conflict_set,
	int vertex_numbers, int mask_numbers,
	std::vector<int> & Delete_the_Row_in_which_Col,
	std::vector<std::list<int>> & Order_of_Row_Deleted_in_Col, 
	int depth, std::vector<int> & MPLD_search_vector, std::string result_file);


/*
 * @brief:	This function decode the result and convert the exact cover problem to origin coloring problem.
 * @param	vertex_numbers:		the number of vertices in the graph
 * @param	mask_numbers:		the number of masks in the MPLD problem
 * @param	result_vec:			a vector storing the all the chosen rows
 * @param	filename:			the file storing the result
 */
void Decode(int vertex_numbers, int mask_numbers, std::vector<int> result_vec, std::set<std::pair<int, int>>  conflict_set, std::string filename);

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