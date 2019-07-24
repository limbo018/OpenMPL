#include"DL_MPL.h"
#include<time.h>
#include <assert.h>
#include <cstdlib>
// Use the API in SimpleMPL.h to create the Component file for dlx solver;
// then call MPLD_Solver()


// ===============================================================================
// InFileName : 		a file stores edge list of the conflict graph.
// ExactCoverProblem : 	a file stores the intermediate dancing link structure information 
// result_file :		a file stores the result

int main()
{
	clock_t start, finish;
	std::string InFileName = "2.txt";
	std::string ExactCoverProblem = "MPLD_total_C1.txt";
	std::string result_file = "total_C1_result.txt";
	std::cout<<"Here"<<std::endl;
	start = clock();
	std::cout<<"Ready to solver"<<std::endl;
	MPLD_Solver(InFileName, ExactCoverProblem, false, 3, result_file);
	finish = clock();
	double totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "\ntime : " << totaltime << std::endl;
	system("pause");
	return 0;
}