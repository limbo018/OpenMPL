/*************************************************************************
    > File Name: main.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri May 22 07:58:01 2015
 ************************************************************************/

#include <iostream>
#include <chrono>
#include <ratio>
#include "SimpleMPL.h"

int main(int argc, char** argv)
{
	SimpleMPL::SimpleMPL mpl;
    mpl.print_welcome();

    // be aware of sizes of every object
    SimpleMPL::mplPrint(SimpleMPL::kDEBUG, "size of Shape = %u bytes\n", sizeof(SimpleMPL::Shape));
    SimpleMPL::mplPrint(SimpleMPL::kDEBUG, "size of Rectangle = %u bytes\n", sizeof(SimpleMPL::Rectangle<int>));
    SimpleMPL::mplPrint(SimpleMPL::kDEBUG, "size of Polygon = %u bytes\n", sizeof(SimpleMPL::Polygon<int>));

  auto start = std::chrono::high_resolution_clock::now();
	mpl.run(argc, argv);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::ratio<1, 1>> duration_s (end - start); 
  mplPrint(SimpleMPL::kINFO, "program takes %g seconds\n", duration_s.count());

	return 0;
}
