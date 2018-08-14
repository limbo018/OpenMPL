#include"DL_MPL.h"

// Use the API in SimpleMPL.h to create the Component file for dlx solver;
// then call MPLD_Solver()

int main()
{
	std::string InFileName = "edge.txt";
	std::string ExactCoverProblem = "MPLD_edge.txt";
	MPLD_Solver(InFileName, ExactCoverProblem, false, 2);
	system("pause");
	return 0;
}