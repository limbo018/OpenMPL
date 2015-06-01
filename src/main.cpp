/*************************************************************************
    > File Name: main.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri May 22 07:58:01 2015
 ************************************************************************/

#include <iostream>
#include <boost/timer/timer.hpp>
#include "SimpleMPL.h"
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
	SimpleMPL::SimpleMPL mpl;
    mpl.print_welcome();

	boost::timer::auto_cpu_timer timer;

	mpl.run(argc, argv);

	return 0;
}
