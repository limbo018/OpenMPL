#include"DL_MPL.h"

// Use the API in SimpleMPL.h to create the Component file for dlx solver;
// then call MPLD_Solver()


// ===============================================================================
// InFileName : 		a file stores edge list of the conflict graph.
// ExactCoverProblem : 	a file stores the intermediate dancing link structure information 
// result_file :		a file stores the result

int DL_main()
{
	std::string InFileName = "karate.txt";
	std::string ExactCoverProblem = "MPLD_karate.txt";
	std::string result_file = "result.txt";
	MPLD_Solver(InFileName, ExactCoverProblem, false, 3, result_file);
	return 0;
}