#include"DL_MPL.h"

int main()
{
	std::string InFileName = "karate.txt";
	std::string ExactCoverProblem = "MPLD_karate.txt";
	MPLD_Solver(InFileName, ExactCoverProblem, false, 3);
	system("pause");
	return 0;
}