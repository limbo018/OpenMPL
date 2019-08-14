#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdexcept>
#include<set>
#include<exception>
#include<list>
#include<map>
#include<limits.h>
#include <assert.h>

struct Cell {
	Cell *Left, *Right, *Up, *Down;
	int Row = 0;
	int Col = 0;
	bool InDLX = true;
	int Children_Number;
};

struct DancingLink{  
	Cell DL_Header;
	Cell* Row_Header_Table;
	Cell* Col_Header_Table;
	int Row_Number_All = 0;
	int Col_Number_All = 0;
	int Row_Number_Now = 0;
	int Col_Number_Now = 0;
};

struct Edge {
	int source;
	int target;
	int No;
};

struct Vertex {
	bool Is_Parent = true;
	Vertex* parent = NULL;
	int No = 0;
	int Stitch_No = -1;
	std::set<Vertex*> Conflicts;
	std::set<Vertex*> Conflicts_in_LG;
	std::vector<Vertex*> Childs;
	void parentOf(Vertex* child){
		/*
		three cases of child:
		1. parent and real node
		2. parent and fake node (represent one feature in layout)
		3. child and real node
		 */
		
		//std::cout<<this->Stitch_No<<" "<<child->Stitch_No<<std::endl;
		// case 3
		if(child->Is_Parent == false){
			//std::cout<<"case 3"<<std::endl;
			assert(child->Childs.empty());
			assert(!child->parent->Childs.empty());
			this->parentOf(child->parent);
			child->parent = this;
		}
		else{
			//insert conflicts into this conflicts if no repeatation
			//case 2
			if(!child->Childs.empty()){
				//std::cout<<"case 2"<<std::endl;
				this->Childs.insert(this->Childs.end(),child->Childs.begin(),child->Childs.end());
				for(std::vector<Vertex*>::iterator it = child->Childs.begin(); it != child->Childs.end(); ++it) {
					assert((*it)->Is_Parent == false);
					assert((*it)->parent == child);
					(*it)->parent = this;
				}
				delete child;
			}
			//case 1
			else{
				//std::cout<<"case 1"<<std::endl;
				child->Is_Parent = false;
				child->parent = this;
				this->Childs.push_back(child);
			}

		}

	}
	//this function is needed cause latter conflicts will be generated when new line reads in
	void updateConflicts(){
		for(std::vector<Vertex*>::iterator it = this->Childs.begin(); it != this->Childs.end(); ++it) {
			this->Conflicts.insert((*it)->Conflicts.begin(),(*it)->Conflicts.end());
		}
		for(std::set<Vertex*>::iterator it = this->Conflicts.begin(); it != this->Conflicts.end(); ++it) {
			if((*it)->Is_Parent){
				this->Conflicts_in_LG.insert((*it));
			}
			else{
				this->Conflicts_in_LG.insert((*it)->parent);
			}
			// if(this->Is_Parent){
			// 	(*it)->Conflicts_in_LG.insert(this);
			// }
			// else{
			// 	(*it)->Conflicts_in_LG.insert(this->parent);
			// }
		}
	}
	void updateDuplicateLGConflicts(){
		for(std::set<Vertex*>::iterator it = this->Conflicts.begin(); it != this->Conflicts.end(); ++it) {
			if(this->Is_Parent){
				(*it)->Conflicts_in_LG.insert(this);
			}
			else{
				(*it)->Conflicts_in_LG.insert(this->parent);
			}
		}
	}
};


struct Edge_Simple {
	int target;
	int No;
};


// Set a vertex 
// Set the Up and Down links to Cell 'c' itself.
void UD_self(Cell &c);

// Set the Left and Right links to Cell 'c' itself.
void LR_self(Cell & c);


// Remove the links from the Up and Down neighbors to Cell 'c'
void UD_remove(Cell & c);

// Remove the links from the Left and Right neighbors to Cell 'c'
void LR_remove(Cell & c);

// Recover the links from the Up and Down neighbors to Cell 'c'
void LR_recover(Cell & c);

// Recover the links from the Up and Down neighbors to Cell 'c'
void UD_recover(Cell & c);

// Print the solution of this dlx problem.
void Print_Result(std::vector<int> result_vec);

/* 
 * @brief : Initialize the dancing links problem.
 * @param dl	dlx object
 * @param row	number of rows in the dlx problem
 * @param col	number of columns in the dlx problem
 *
*/
void DL_Init(DancingLink & dl, int row, int col);

/*
* Insert a cell into the specified location (row and column).
* @param row	the row label
* @param col	the column label
* @param dl		the dancing link problem instance 
*/
void Cell_Insert(DancingLink & dl, int row, int col);

/*
 * @brief : Finish the construction of DLX through calling Cell_Insert().
 *			Read data from the given file and insert data into the DaningLink object.
 * 
 * @param filename	The given file. The data in it should be in the following format:
 *					- The 1st line : [Row , Col] . Row : the number of rows in the dlx. 
 *												   Col : the number of columns in the dlx.
 *					- The 2nd line : [N] .	N : the number of cells in the dlx.
 *					- The following [N] lines : [row_No, col_No] . The location of each cell. The index starts at 1.
 */
void DL_Load(DancingLink &dl, std::string filename);

/*
 * @brief : Select the column with least cells and cover it in the next iteration.
 * @param dl	the dancing link problem instance
 *
 * @return
 *		- the Column No.
*/
int Select_Next_Column_Sorting(DancingLink & dl);

/*
 * @brief : Select the column simply by get the header's right neighbor.
 */
int Select_Next_Column_Simply(DancingLink & dl);

/*
 * @brief : Due the selection of one row in an iteration, we should remove some relevant colums and rows meanwhile.
 *			For easy understanding and implementations, firstly, we determine the rows and columns we should remove.
 *
 * @oaram dl			the dancing link problem instance
 * @param target_row	the row we select in this iteration
 * @param row_set		a set containing all rows we should remove. We use <set> because each element in a set must be unique.
 * @param col_set		a set containing all cols we should remove. We use <set> because each element in a set must be unique.
 * 
 */
void Select_All_Rows_Cols(DancingLink & dl, int target_row, std::set<int> & row_set, std::set<int> & col_set);

/*
 * @brief : Remove one row from the dlx problem in the present iteration temporarily, 
 *			according to the row No. (it may be recovered in the subsequent iterations.)
 *
 * @param dl	the dancing link problem instance
 * @param row	the No. of the row we remove in this iteration.
*/
void Remove_Single_Row(DancingLink & dl, int row);

/*
 * @brief : Remove one row from the dlx problem in the present iteration temporarily, 
 *			according to the row No. (it may be recovered in the subsequent iterations.)
 *
 * @param dl	the dancing link problem instance
 * @param row	the No. of the row we remove in this iteration.
*/
void Remove_Single_Col(DancingLink & dl, int col);

/*
 * @brief : Remove all the rows we select in one iteration, also modify the headers of the corresponding columns.
 *			This function relys on the result of 'Select_All_Rows_Cols()'.
 *			This function calls 'Remove_Single_Row()' repeatedly to remove all the rows.
 * 
 * @param dl			the dancing link problem instance.
 * @param row_set		a set containing all rows we should remove. We use <set> because each element in a set must be unique.
 * @param col_set		a set containing all cols we should remove. We use <set> because each element in a set must be unique.
 */
void Remove_Rows_Cols(DancingLink & dl, std::set<int> row_set, std::set<int> col_set);

/*
 * @brief : Recover one row. It's the reverse process of 'Remove_Single_Row()'.
 *
 * @param dl	the dancing link problem instance
 * @param row	the No. of the row we remove
 */
void Recover_Single_Row(DancingLink & dl, int row);

/*
 * @brief : Recover all the rows we select in the previous iteration, also recover the corresponding column headers.
 *			It's the reverse process of 'Remove_Rows_Cols()'.
 */
void Recover_Rows_Cols(DancingLink & dl, std::set<int> row_set, std::set<int> col_set);

/*
 * @brief : The DLX solver, solves the dlx problem recursively.
 * 
 * @param dl							the dlx problem instance
 * @param result_vec					a vector containing the solution in present recursive step.
 * @param Select_Next_Column_Method		a function pointer, used to configure how to select the next column
 *
 * @return
 *			bool	whether the dlx instance has valid solution in present step.
 */
bool DancingLinksSolver(DancingLink & dl, std::vector<int> & result_vec, int(*Select_Next_Column_Method)(DancingLink&));