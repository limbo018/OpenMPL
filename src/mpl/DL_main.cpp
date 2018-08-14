#include"DL_MPL.h"

int main()
{
	std::string InFileName = "edge.txt";
	std::string ExactCoverProblem = "MPLD_edge.txt";
	MPLD_Solver(InFileName, ExactCoverProblem, false);
	system("pause");
	return 0;
}