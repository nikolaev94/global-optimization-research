
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
		user_parameters.get_solving_method(),
		user_parameters.do_use_neighbour_nodes_optimization() };

	Solver::Output output(input);

	Solver slv(input, problems);

	slv.run_solver(output);

	output.dump_results_to_file(prefix + "results.csv");

	output.dump_error_metrics_by_trials_to_file(prefix + "errors.csv");

	output.dump_solved_problem_portion_by_trials_to_file(prefix + "portions.csv");

	output.dump_method_trials_to_file(prefix + "trials.txt");

	return 0;
}
