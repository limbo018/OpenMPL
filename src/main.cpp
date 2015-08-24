/*************************************************************************
    > File Name: main.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri May 22 07:58:01 2015
 ************************************************************************/

#include <iostream>
#include <boost/timer/timer.hpp>
#include "SimpleMPL.h"

int main(int argc, char** argv)
{
	SimpleMPL::SimpleMPL mpl;
    mpl.print_welcome();

    char buf[256];
    SimpleMPL::mplSPrint(SimpleMPL::kINFO, buf, "program takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer (buf);

	mpl.run(argc, argv);

	return 0;
}
