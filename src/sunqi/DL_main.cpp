#include"DL_MPL.h"


int main()
{
	clock_t start, finish;
	std::string InFileName = "total_C1.txt";
	std::string ExactCoverProblem = "MPLD_total_C1.txt";
	start = clock();
	MPLD_Solver(InFileName, ExactCoverProblem, false, 3);
	finish = clock();
	double totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "\ntime : " << totaltime << std::endl;
	system("pause");
	return 0;
}