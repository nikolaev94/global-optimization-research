
#ifndef COUNTER_PLOT_DRAWER_H
#define COUNTER_PLOT_DRAWER_H

#include <fstream>
#include <vector>

#include <opt_problem.h>

#include "user_input.h"

class ContourPlotDataDumper
{
private:
	OptProblem::OptProblemP target_problem;

	const unsigned int DEFAULT_GRID_SIZE = 60;

public:
	ContourPlotDataDumper(OptProblem::OptProblemP in_target_problem) :
		target_problem(in_target_problem) {}

	void dump(const std::string&);
};

#endif
