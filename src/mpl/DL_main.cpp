#include"DL_MPL.h"

// Use the API in SimpleMPL.h to create the Component file for dlx solver;
// then call MPLD_Solver()

int main()
{
	std::string InFileName = "karate.txt";
	std::string ExactCoverProblem = "MPLD_karate.txt";
	MPLD_Solver(InFileName, ExactCoverProblem, false, 3);
	system("pause");
	return 0;
}