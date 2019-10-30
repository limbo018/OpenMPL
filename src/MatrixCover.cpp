#include "MatrixCover.h"

void delete_rows_and_columns(int **dl_matrix, int *deleted_rows,
                             int *deleted_cols, int search_depth,
                             int selected_row_id, int total_dl_matrix_row_num,
                             int total_dl_matrix_col_num) {
  for (int i = 0; i < total_dl_matrix_col_num; i++) {
    if (dl_matrix[selected_row_id][i] == 1 &&
        deleted_cols[i] ==
            0) { // we only delete rows that are not deleted or removed
      deleted_cols[i] = search_depth;
      for (int j = 0; j < total_dl_matrix_row_num; j++) {
        if (dl_matrix[j][i] == 1 && deleted_rows[j] == 0) {
          deleted_rows[j] = search_depth;
        }
      }
    }
  }
}

void init_vectors(int *vec, int vec_length, int value = 0) {
  for (int i = 0; i < vec_length; i++) {
    vec[i] = value;
  }
}

int get_largest_value(int *vec, int vec_length) {
  int max_value = *std::max_element(vec, vec + vec_length);
  for (int i = vec_length - 1; i >= 0; i--) {
    if (vec[i] == max_value) {
      return i;
    }
  }
}

void init_vectors_reserved(int *vec, int vec_length, int value = 0) {
  for (int i = 0; i < vec_length; i++) {
    if (vec[i] != -1) {
      vec[i] = value;
    }
  }
}

bool check_existance_of_candidate_rows(int *deleted_rows, int *row_group,
                                       int search_depth,
                                       int total_dl_matrix_row_num) {
  for (int i = 0; i < total_dl_matrix_row_num; i++) {
    // std::cout << deleted_rows[i] << ' ' << row_group[i] << std::endl;
    if (deleted_rows[i] == 0 && row_group[i] == search_depth) {
      // std::cout<<"Candidate Row Found...."<<std::endl;
      return true;
    }
  }
  return false;
}

// instead we have to row_group nodes
void get_vertex_row_group(int *row_group, int **dl_matrix, int vertex_num,
                          int total_dl_matrix_row_num) {
  for (int i = 0; i < vertex_num; i++) {
    for (int j = 0; j < total_dl_matrix_row_num; j++) {
      row_group[j] += dl_matrix[j][i] * (i + 1);
    }
  }
}

int select_row(int *deleted_rows, int *row_group, int search_depth,
               int total_dl_matrix_row_num) {
  for (int i = 0; i < total_dl_matrix_row_num; i++) {
    if (deleted_rows[i] == 0 && row_group[i] == search_depth) {
      return i;
    }
  }
}

void recover_deleted_rows(int *deleted_rows, int search_depth,
                          int total_dl_matrix_row_num) {
  for (int i = 0; i < total_dl_matrix_row_num; i++) {
    if (abs(deleted_rows[i]) > search_depth ||
        deleted_rows[i] == search_depth) {
      deleted_rows[i] = 0;
    }
  }
}

// To be optimized
void recover_deleted_cols(int *deleted_cols, int search_depth,
                          int total_dl_matrix_col_num) {
  for (int i = 0; i < total_dl_matrix_col_num; i++) {
    if (deleted_cols[i] >= search_depth) {
      deleted_cols[i] = 0;
    }
  }
}

void recover_results(int *results, int search_depth,
                     int total_dl_matrix_row_num) {
  for (int i = 0; i < total_dl_matrix_row_num; i++) {
    if (results[i] == search_depth) {
      results[i] = 0;
    }
  }
}

// need to optimized to map on GPU array
int get_conflict_node_id(int *deleted_rows, int *row_group, int search_depth,
                         int total_dl_matrix_row_num) {
  int conflict_node_id = 0;
  for (int i = 0; i < total_dl_matrix_row_num; i++) {
    if (row_group[i] == search_depth + 1) {
      if (conflict_node_id < deleted_rows[i] && deleted_rows[i] < search_depth+1) {
        conflict_node_id = deleted_rows[i];
      }
    }
  }
  return conflict_node_id;
}

int get_conflict_col(int **dl_matrix, int *deleted_rows, int *deleted_cols,
                     int *row_group, int conflict_node_id, int search_depth,
                     int vertex_num, int total_dl_matrix_row_num,
                     int total_dl_matrix_col_num) {
  int conflict_col_id = 0;
  int idxa = 0;
  int idxb = 0;
  for (int i = 0; i < total_dl_matrix_row_num;
       i++) { // find the conflict edge that connects current node and the most
              // closest node.
    if (deleted_rows[i] == -conflict_node_id) {
      idxa = i;
    }
    if (row_group[i] == search_depth + 1 &&
        deleted_rows[i] == conflict_node_id) {
      idxb = i;
    }
    // if (row_group[i] == search_depth + 1 &&
    //   deleted_rows[i] == conflict_node_id) {
    //  for (int j = total_dl_matrix_col_num - 1; j > vertex_num; j--) {
    //    if (dl_matrix[i][j] * deleted_cols[j] == conflict_node_id) {
    //      conflict_col_id = j;
    //    }
    //  }
    //}
  }
#ifdef BENCHMARK
  std::cout<< "conflict rows are"<<std::endl;
  print_vec(dl_matrix[idxa], total_dl_matrix_col_num);
  print_vec(dl_matrix[idxb], total_dl_matrix_col_num);
#endif
  for (int j = vertex_num; j < total_dl_matrix_col_num; j++) {
    if (dl_matrix[idxa][j] == dl_matrix[idxb][j] && dl_matrix[idxa][j]==1) {
      if(conflict_col_id<j){
        conflict_col_id = j;
      }

    }
  }
  return conflict_col_id;
}

void remove_cols(int *deleted_cols, int *col_group, int conflict_col_id,
                 int total_dl_matrix_col_num) {
  for (int i = 0; i < total_dl_matrix_col_num; i++) {
    if (col_group[i] == col_group[conflict_col_id]) {
      deleted_cols[i] = -1;
    }
  }
}

void print_vec(int *vec, int vec_length) {
  for (int i = 0; i < vec_length; i++) {
    std::cout << vec[i] << ' ';
  }
  std::cout << std::endl;
}

void print_matrix(int **mat, int row, int col) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      std::cout << mat[i][j] << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void mc_solver(int **dl_matrix, int *results, int *deleted_cols, int *col_group,
               int vertex_num, int total_dl_matrix_row_num,
               int total_dl_matrix_col_num) {
  // to be refreshed if one conflict reaches many counts
  int search_depth = 0;
  int selected_row_id = 0;
  int *conflict_count = new int[total_dl_matrix_col_num];

  int *deleted_rows = new int[total_dl_matrix_row_num];

  int *vertices_covered = new int[vertex_num];
  int *row_group = new int[total_dl_matrix_row_num];
  // int *col_group = new int[total_dl_matrix_col_num];
  int selected_row_id_in_previous_search;
  int conflict_node_id;
  int conflict_col_id;
  int hard_conflict_threshold = 500;
  bool token;

  // init lots of vectors
  init_vectors(conflict_count, total_dl_matrix_col_num);
  init_vectors(deleted_cols, total_dl_matrix_col_num);
  init_vectors(deleted_rows, total_dl_matrix_row_num);
  init_vectors(results, total_dl_matrix_row_num);
  init_vectors(row_group, total_dl_matrix_row_num);
  get_vertex_row_group(row_group, dl_matrix, vertex_num,
                       total_dl_matrix_row_num);

#ifdef BENCHMARK
  char tmp;
  print_vec(row_group, total_dl_matrix_row_num);
  print_vec(col_group, total_dl_matrix_col_num);
  print_matrix(dl_matrix, total_dl_matrix_row_num, total_dl_matrix_col_num);
#endif

  for (search_depth = 1; search_depth <= vertex_num;) {
    token=check_existance_of_candidate_rows(
            deleted_rows, row_group, search_depth,
            total_dl_matrix_row_num);
#ifdef BENCHMARK
    std::cout << "Deleted Cols" << std::endl;
    print_vec(deleted_cols, total_dl_matrix_col_num);
    std::cout << "Deleted Rows" << std::endl;
    print_vec(deleted_rows, total_dl_matrix_row_num);
    std::cout << "Conflict Count" << std::endl;
    print_vec(conflict_count, total_dl_matrix_col_num);
    std::cout << "Results" << std::endl;
    print_vec(results, total_dl_matrix_row_num);
    std::cin >> tmp;
    std::cout << "search depth is " << search_depth << std::endl;
    std::cout << "existence of candicate row: " << token << std::endl;
#endif

    if (token) { // check if there are candidate rows
                                        // existing, if no, do backtrace
      selected_row_id =
          select_row(deleted_rows, row_group, search_depth,
                     total_dl_matrix_row_num); // select row and add to results
#ifdef BENCHMARK
      std::cout << "selected id is " << selected_row_id << std::endl;
#endif
      results[selected_row_id] = search_depth;
      delete_rows_and_columns(
          dl_matrix, deleted_rows, deleted_cols, search_depth, selected_row_id,
          total_dl_matrix_row_num,
          total_dl_matrix_col_num); // delete covered rows and columns
      deleted_rows[selected_row_id] = -search_depth;
      search_depth++; // next step
      // print_vec(deleted_cols, total_dl_matrix_col_num);
      // print_vec(deleted_rows, total_dl_matrix_row_num);
      // print_vec(conflict_count, total_dl_matrix_col_num);
      // print_vec(results, total_dl_matrix_row_num);

    } else { // do backtrace
      search_depth--;
#ifdef BENCHMARK
      std::cout << "search depth is" << search_depth << std::endl;
#endif      
      if (search_depth > 0) {
        conflict_node_id = get_conflict_node_id(
            deleted_rows, row_group, search_depth, total_dl_matrix_row_num);
        if(conflict_node_id>0){
          conflict_col_id = get_conflict_col(
              dl_matrix, deleted_rows, deleted_cols, row_group, conflict_node_id,
              search_depth, vertex_num, total_dl_matrix_row_num,
              total_dl_matrix_col_num); // get conflict edge
#ifdef BENCHMARK
          std::cout << "conflict node id is " << conflict_node_id << std::endl;
          std::cout << "conflict col id is " << conflict_col_id << std::endl;
          //if(conflict_col_id==0){std::cin>>tmp;}
#endif
          conflict_count[conflict_col_id]++; // update conflict edge count
          recover_deleted_rows(deleted_rows, search_depth,
                              total_dl_matrix_row_num); // recover deleted rows
                                                        // previously selected
                                                        // rows
          recover_deleted_cols(deleted_cols, search_depth,
                              total_dl_matrix_col_num); // recover deleted cols
                                                        // except afftected by
                                                        // previously selected
                                                        // rows
          recover_results(results, search_depth,
                          total_dl_matrix_row_num); // recover results

          if (conflict_count[conflict_col_id] > hard_conflict_threshold) {
            search_depth = 1;
            init_vectors(conflict_count, total_dl_matrix_col_num);
            init_vectors_reserved(deleted_cols, total_dl_matrix_col_num);
            init_vectors(deleted_rows, total_dl_matrix_row_num);
            init_vectors(results, total_dl_matrix_row_num);
            remove_cols(deleted_cols, col_group, conflict_col_id,
                        total_dl_matrix_col_num);
            deleted_cols[conflict_col_id] = -1;
          }
        } else {
          recover_deleted_rows(deleted_rows, search_depth,
                              total_dl_matrix_row_num); // recover deleted rows
                                                        // previously selected
                                                        // rows
          recover_deleted_cols(deleted_cols, search_depth,
                              total_dl_matrix_col_num); // recover deleted cols
                                                        // except afftected by
                                                        // previously selected
                                                        // rows
          recover_results(results, search_depth,
                          total_dl_matrix_row_num); // recover results
        }
      } else { // if all vertices are gone through, directly remove the edge
               // with largest conflict count.
        search_depth = 1;
        conflict_col_id =
            get_largest_value(conflict_count, total_dl_matrix_col_num);
        init_vectors(conflict_count, total_dl_matrix_col_num);
        init_vectors_reserved(deleted_cols, total_dl_matrix_col_num);
        init_vectors(deleted_rows, total_dl_matrix_row_num);
        init_vectors(results, total_dl_matrix_row_num);
        remove_cols(deleted_cols, col_group, conflict_col_id,
                    total_dl_matrix_col_num);
      }
    }
  }
#ifdef BENCHMARK
  print_matrix(dl_matrix, total_dl_matrix_row_num, total_dl_matrix_col_num);
#endif
  delete[] deleted_rows;
  delete[] row_group;
  delete[] vertices_covered;
  delete[] conflict_count;
}
