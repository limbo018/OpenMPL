#include <stdlib.h>
#include <algorithm>
#include <iostream>
//exact cover matrix
//cols are vertices and all possible conflicts
//row are all vertices with possible coloring options
/* a b c ab1 ab2 ac1 ac2
a1|1     1       1
a2|1         1       1
b1| 1    1   
b2| 1        1
c1|   1          1
c2|   1              1
 */
//int **dl_matrix;

//[0, 1, 1, 2, 1, 2, 3]
// equal to the number of rows in dl_matrix
// 0 means the corresponding rows are not removed
// i>0 means the corresponding rows are removed in the ith iteration
//int *selected_rows;

//[0, 1, 2, 3]
// equal to the number of vertices
// 0 means the column is not covered
// i>0 means the column is covered in the ith iteration
//int *covered_cols;

// similarly defined above
//int *deleted_rows;
//int *deleted_cols;

// rows selected in which iteration
//int selected_row_id
//int *results;

// global iterator in dancing link search

// int search_depth

//

//Operation of delete rows and columns
void delete_rows_and_columns(
    int **dl_matrix,
    int *deleted_rows,
    int *deleted_cols,
    int search_depth,
    int selected_row_id,
    int total_dl_matrix_row_num,
    int total_dl_matrix_col_num);

//initialize vectors to #value
void init_vectors(int *vec, int vec_length, int value);

//
int get_largest_value(int *vec, int vec_length);

//
bool check_existance_of_candidate_rows(int *deleted_rows, int *row_group, int search_depth, int total_dl_matrix_row_num);

// row_grouping rows by different vertices e.g. [1 1 2 2 3 3 4 4 4 4]
void get_vertex_row_group(int *row_group, int **dl_matrix, int vertex_num, int total_dl_matrix_row_num);

//select one row
int select_row(int *deleted_rows, int *row_group, int search_depth, int total_dl_matrix_row_num);

//
void recover_deleted_rows(int *deleted_rows, int search_depth, int total_dl_matrix_row_num);

//
void recover_deleted_cols(int **dl_matrix, int *deleted_cols, int search_depth, int total_dl_matrix_col_num);

//
void recover_results(int *results, int search_depth, int total_dl_matrix_row_num);

//
int get_conflict_node_id(int *deleted_rows, int *row_group, int search_depth, int total_dl_matrix_row_num);

//
int get_conflict_col(int **dl_matrix, int *deleted_rows, int *deleted_cols, int *row_group, int conflict_node_id, int search_depth, int vertex_num, int total_dl_matrix_row_num, int total_dl_matrix_col_num);

//
void mc_solver(int **dl_matrix, int *results, int *deleted_cols, int *col_group, int vertex_num, int total_dl_matrix_row_num, int total_dl_matrix_col_num);
