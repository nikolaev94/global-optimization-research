
#include <cstdlib>

#include <iostream>
#include <list>
#include <memory>
#include <fstream>
#include <string>

#include <solver.h>

#include <discpp.h>

#include "problem_loader.h"
#include "user_input.h"

int main(int argc, char* argv[])
{
	UserParameters user_parameters;

	user_parameters.parse_arguments_from_command_line(argc, argv);

	ProblemLoader loader;

	Solver::problem_list problems;

	loader.load_series(user_parameters, problems);

	std::string prefix = loader.get_output_file_prefix(user_parameters);

	Solver::Input input = { UserParameters::LEFT_BOUND, UserParameters::RIGHT_BOUND,
		user_parameters.get_precision(), user_parameters.get_method_parameter(),
		user_parameters.get_workers_number(), user_parameters.get_dimension(),
		user_parameters.get_solving_method() };

	Solver::Output output(input);

	Solver slv(input, problems);

	slv.run_solver(output);

	output.dump_results_to_file(prefix + "results.csv");

	output.dump_error_metrics_by_trials_to_file(prefix + "errors.csv");

	output.dump_solved_problem_portion_by_trials_to_file(prefix + "portions.csv");

	output.dump_method_trials_to_file(prefix + "trials.txt");

	/*
	double xray[50], yray[50], zmat[50][50];
	int n = 50, i, j;
	double  fpi = 3.14159 / 180.0, step, x, y;
	double  zlev;
	Dislin g;

	step = 360.0 / (n - 1);

	for (i = 0; i < n; i++)
	{
		xray[i] = i * step;
		yray[i] = i * step;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			x = xray[i] * fpi;
			y = yray[j] * fpi;
			zmat[i][j] = 2 * sin(x) * sin(y);
		}
	}

	g.scrmod("revers");
	g.setpag("da4p");
	g.metafl("cons");
	g.disini();
	g.complx();
	g.pagera();

	g.titlin("Contour Plot", 1);
	g.titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)", 3);

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.intax();
	g.axspos(450, 2670);
	g.graf(0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0);

	g.height(30);
	for (i = 0; i < 9; i++)
	{
		zlev = -2.0 + i * 0.5;
		g.setclr((i + 1) * 25);
		if (i == 4)
			g.labels("none", "contur");
		else
			g.labels("float", "contur");

		g.contur(xray, n, yray, n, (double *)zmat, zlev);
	}

	g.height(50);
	g.color("fore");
	g.title();
	g.disfin();
	return 0;*/
	return 0;
}
