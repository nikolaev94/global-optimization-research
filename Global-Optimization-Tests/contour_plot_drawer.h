
#ifndef COUNTER_PLOT_DRAWER_H
#define COUNTER_PLOT_DRAWER_H

#include <vector>

#include <discpp.h>

#include <opt_problem.h>


#include "user_input.h"

class ContourPlotDrawer
{
private:
	OptProblem::OptProblemP target_problem;

public:
	ContourPlotDrawer(OptProblem::OptProblemP in_target_problem) :
		target_problem(in_target_problem) {}

	void draw();
};

#endif
