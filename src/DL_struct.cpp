#include"DL_struct.h"

void UD_self(Cell &c) {
	c.Up = &c;
	c.Down = &c;
}

void LR_self(Cell & c) {
	c.Left = c.Right = &c;
}

void UD_remove(Cell & c)
{
	if(c.InDLX){
	c.Up->Down = c.Down;
	c.Down->Up = c.Up;
	c.InDLX = false;
	}

}

void LR_remove(Cell & c)
{
	// std::cout<<c.InDLX<<std::endl;
	// std::cout<<c.Children_Number<<std::endl;
	// std::cout<<c.Row<<","<<c.Col<<std::endl;
	// std::cout<<c.Left->Row<<","<<c.Left->Col<<","<<std::endl;
	// std::cout<<c.Right->Row<<","<<c.Left->Col<<","<<std::endl;
	if(c.InDLX){
		c.Left->Right = c.Right;
		c.Right->Left = c.Left;
		c.InDLX = false;
	}
}

void LR_recover(Cell & c) {
	if(c.InDLX == false){
		c.Left->Right = &c;
		c.Right->Left = &c;
		c.InDLX = true;
	}
}

void UD_recover(Cell & c) {
	if(c.InDLX == false){
		c.Up->Down = &c;
		c.Down->Up = &c;
		c.InDLX = true;
	}

}

void Print_Result(std::vector<int> & result_vec) {
	std::vector<int>::iterator it;
	std::cout << "Result : ";
	for (it = result_vec.begin(); it != result_vec.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
}

void DL_Init(DancingLink & dl, int & row, int & col)
{
	dl.Row_Number_Now = dl.Row_Number_All = row;
	dl.Col_Number_All = dl.Col_Number_Now = col;

	dl.DL_Header.Col = 0;
	dl.DL_Header.Row = 0;
	
	dl.Col_Header_Table = new Cell[col + 1];
	for (int i = 1; i < col; i++)	dl.Col_Header_Table[i].Right = &dl.Col_Header_Table[i + 1];
	for (int i = 2; i <= col; i++)	dl.Col_Header_Table[i].Left = &dl.Col_Header_Table[i - 1];
	for (int i = 1; i <= col; i++)
	{
		UD_self(dl.Col_Header_Table[i]);
		dl.Col_Header_Table[i].Children_Number = 0;
		dl.Col_Header_Table[i].Col = i;
	}
	dl.DL_Header.Right = &dl.Col_Header_Table[1];
	dl.DL_Header.Left = &dl.Col_Header_Table[col];
	dl.Col_Header_Table[col].Right = &dl.DL_Header;
	dl.Col_Header_Table[1].Left = &dl.DL_Header;
	
	dl.Row_Header_Table = new Cell[row + 1];
	for (int i = 1; i < row; i++)	dl.Row_Header_Table[i].Down = &dl.Row_Header_Table[i + 1];
	for (int i = 2; i <= row; i++)	dl.Row_Header_Table[i].Up = &dl.Row_Header_Table[i - 1];
	for (int i = 1; i <= row; i++)
	{
		LR_self(dl.Row_Header_Table[i]);
		dl.Row_Header_Table[i].Children_Number = 0;
		dl.Row_Header_Table[i].Row = i;
	}
	
	dl.DL_Header.Down = &dl.Row_Header_Table[1];
	dl.DL_Header.Up = &dl.Row_Header_Table[row];
	dl.Row_Header_Table[row].Down = &dl.DL_Header;
	dl.Row_Header_Table[1].Up = &dl.DL_Header;
	
	/*
	std::cout << &dl << std::endl;
	std::cout << &dl.Col_Header_Table[3] << std::endl;
	std::cout << &dl.Col_Header_Table[3].Children_Number << std::endl;
	*/
}

/*
void print(DancingLink & dl)
{
	for (Cell *it = dl.DL_Header.Right; it != &dl.DL_Header; it = it->Right) {
		std::cout << "\ncol : " << (*it).Col << std::endl;
		for (Cell *j = (*it).Down; j != it; j = j->Down)
			std::cout << "Row : " << (*j).Row << " ";
	}
	std::cout << std::endl;
}
*/

void Cell_Insert(DancingLink & dl, uint32_t  row, uint32_t  col) {
	Cell *r = &dl.Row_Header_Table[row];
	for (Cell *i = r->Right; i != r; i = i->Right) {
		if(i->Col == col) return;
	}
	Cell *c = new Cell;
	c->Row = (int)row;
	c->Col = (int)col;
	assert(row>0);
	assert(col>0);
	// std::cout << "Row : " << row << " Col : " << col << std::endl;
	// c.Col_Header = &dl.Col_Header_Table[col];
	// std::cout << "&dl\t" << &dl << std::endl;
	// std::cout << "&dl.Col_Header_Table[" << col << "] \t" << &dl.Col_Header_Table[col] << std::endl;
	// std::cout << "&dl.Col_Header_Table[" << col << "].Children_Number\t" << &dl.Col_Header_Table[col].Children_Number << std::endl;
	// std::cout << "children_number : \t\t" << dl.Col_Header_Table[col].Children_Number << std::endl;
	
	dl.Col_Header_Table[col].Children_Number++;
	c->Up = dl.Col_Header_Table[col].Up;
	c->Up->Down = c;
	dl.Col_Header_Table[col].Up = c;
	c->Down = &dl.Col_Header_Table[col];
	
	dl.Row_Header_Table[row].Children_Number++;
	c->Right = &dl.Row_Header_Table[row];
	c->Left = dl.Row_Header_Table[row].Left;
	c->Left->Right = c;
	dl.Row_Header_Table[row].Left = c;
}

void DL_Load(DancingLink & dl, std::string & filename) {
	std::ifstream filein(filename.c_str());
	try {
		if (!filein)
			throw std::runtime_error("No such file.");
	}
	catch (std::exception &ex) {
		std::cout << "File Error" << std::endl;
	}
	int Row, Col, N;
	int row_No, col_No;
	filein >> Row;
	filein >> Col;
	filein >> N;
	//std::cout << "total : " << Row << " rows " << Col << " cols " << std::endl;
	DL_Init(dl, Row, Col);
	for (int i = 1; i <= N; i++)
	{
		filein >> row_No >> col_No;
		Cell_Insert(dl, row_No, col_No);
	}
	return;
}

int Select_Next_Column_Sorting(DancingLink & dl) {
	Cell * c = &dl.DL_Header;
	int m = INT_MAX;
	int location = dl.DL_Header.Right->Col;
	for (Cell *i = c->Right; i != c; i = i->Right) {
		if (i->Children_Number <= m) location = i->Col;
	}
	return location;
}

int Select_Next_Column_Simply(DancingLink & dl) {
	// std::cout << "dl.DL_Header.Right->Col : " << dl.DL_Header.Right->Col << std::endl;
	return dl.DL_Header.Right->Col;
}

void Select_All_Rows_Cols(DancingLink & dl, int & target_row, std::set<int> & row_set, std::set<int> & col_set) {
	Cell *r = &dl.Row_Header_Table[target_row];
	col_set.insert(r->Right->Col);
	row_set.insert(target_row);
	
	for (Cell *i = r->Right; i != r; i = i->Right) {
		col_set.insert(i->Col);
		for (Cell *j = dl.Col_Header_Table[i->Col].Down; j != &dl.Col_Header_Table[i->Col]; j = j->Down) {
			row_set.insert(j->Row);
		}
	}
}

void Remove_Single_Row(DancingLink & dl,const int & row) {
	Cell *r = &dl.Row_Header_Table[row];
	UD_remove(*r);
	r->Children_Number = 0;
	// UD_remove(f);
	for (Cell *j = r->Right; j != r; j = j->Right) {
		UD_remove(*j);
		dl.Col_Header_Table[j->Col].Children_Number--;
	}
}

void Remove_Single_Col(DancingLink & dl,const int & col) {
	Cell *c = &dl.Col_Header_Table[col];
	LR_remove(*c);
	c->Children_Number = 0;
	// UD_remove(f);
	for (Cell *j = c->Down; j != c; j = j->Down) {
		LR_remove(*j);
		dl.Row_Header_Table[j->Row].Children_Number--;
	}
}

void Remove_Rows_Cols(DancingLink & dl, std::set<int> & row_set, std::set<int> & col_set) {
	std::set<int>::iterator it;
	for (it = row_set.begin(); it != row_set.end(); it++)
		Remove_Single_Row(dl, *it);
	for (it = col_set.begin(); it != col_set.end(); it++){
		 LR_remove(dl.Col_Header_Table[*it]);
		 }
		//Remove_Single_Col(dl,*it);
	/*
	std::cout << "Col_header : " << std::endl;
	for (Cell *it = dl.DL_Header.Right; it != &dl.DL_Header; it = it->Right)
		std::cout << (*it).Col << " ";
	std::cout << "\nRow_Header : " << std::endl;
	for (Cell *it = dl.DL_Header.Down; it != &dl.DL_Header; it = it->Down)
		std::cout << (*it).Row << std::endl;
	std::cout << std::endl;
	*/
}

void Recover_Single_Row(DancingLink & dl, const int & row) {
	Cell * f = &dl.Row_Header_Table[row];
	UD_recover(*f);
	// UD_recover(*f);
	// dl.Col_Header_Table[f->Col].Children_Number++;
	for(Cell * i = f->Right; i != f; i = i->Right){
		UD_recover(*i);
		dl.Col_Header_Table[i->Col].Children_Number++;
	}
}

void Recover_Rows_Cols(DancingLink & dl, std::set<int> & row_set, std::set<int>  & col_set) {
	std::set<int>::iterator it;
	for (it = row_set.begin(); it != row_set.end(); it++)
		Recover_Single_Row(dl, *it);
	for (it = col_set.begin(); it != col_set.end(); it++)
		LR_recover(dl.Col_Header_Table[*it]);
}

bool DancingLinksSolver(DancingLink & dl, std::vector<int> & result_vec, int (*Select_Next_Column_Method)(DancingLink&)) {
	// print(dl);
	if (dl.DL_Header.Right == &dl.DL_Header)
	{
		Print_Result(result_vec);
		return true;
	}
	int c = (*Select_Next_Column_Method)(dl);
	// std::cout << "current column : " << c << std::endl;
	Cell *col = &dl.Col_Header_Table[c];
	for (Cell * j = col->Down; j != col; j = j->Down) {
		std::set<int> row_set;
		std::set<int> col_set;
		// std::cout << "j->Row : \t" << j->Row << std::endl;
		result_vec.push_back(j->Row);
		Select_All_Rows_Cols(dl, j->Row, row_set, col_set);
		
		/*
		std::cout << "removed_row : ";
		for (std::set<int>::iterator it = row_set.begin(); it != row_set.end(); it++)
			std::cout << *it << " ";
		std::cout << "\nremoved_col : ";
		for (std::set<int>::iterator it = col_set.begin(); it != col_set.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl << std::endl;
		*/

		Remove_Rows_Cols(dl, row_set, col_set);
		if (DancingLinksSolver(dl, result_vec, Select_Next_Column_Method))
			return true;
		result_vec.pop_back();
		Recover_Rows_Cols(dl, row_set, col_set);
	}
	return false;
}
/*
int main()
{
	std::string Temp_File_Name = "in";
	int k = 3;
	int i = 1;
	while (i <= 3) {
		std::cout << "***************************************" << std::endl;
		std::cout << "Enter file " << i << std::endl;
		std::vector<int> result_vec;
		std::string filename = Temp_File_Name + std::to_string(i) + ".txt";
		DancingLink dl;
		DL_Load(dl, filename);
		if (DancingLinksSolver(dl, result_vec, Select_Next_Column_Simply))
			std::cout << "File " << i << " done....." << std::endl;
		else
			std::cout << "File " << i << " has no solution..." << std::endl;
		std::cout << "***************************************" << std::endl;
		i++;
	}
	system("pause");
	return 0;
}
*/