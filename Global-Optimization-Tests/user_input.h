#ifndef USER_INPUT_H
#define USER_INPUT_H

#include <iostream>
#include <string>

#include <solver.h>

enum FunctionClass
{
	GKLS_SIMPLE, GKLS_HARD, VAGRIS
};

struct UserParameters
{
private:
	unsigned int dimensions;
	unsigned int num_workers;
	unsigned int series_size;
	unsigned int problem_no;
	unsigned int num_constraints;
	double method_parameter;
	double precision;
	FunctionClass function_class;
	Solver::SolvingMethod solving_method;
	bool use_neighbour_nodes_optimization;
	bool single_problem_mode;
	bool use_sleep_mode;

	void show_usage();

public:
	static double LEFT_BOUND;

	static double RIGHT_BOUND;

	UserParameters() : dimensions(2), num_workers(2), method_parameter(5.0),
		num_constraints(0), precision(0.1), series_size(10), problem_no(0),
		function_class(FunctionClass::GKLS_SIMPLE), solving_method(Solver::SEQUENTIAL),
		use_neighbour_nodes_optimization(true),single_problem_mode(false),
		use_sleep_mode(false) {}

	void parse_arguments_from_command_line(int argc, char* argv[]);
	unsigned int get_num_constrains() const;
	unsigned int get_dimension() const;
	unsigned int get_series_size() const;
	unsigned int get_selected_problem_no() const;
	FunctionClass get_function_class() const;
	unsigned int get_workers_number() const;
	double get_method_parameter() const;
	Solver::SolvingMethod get_solving_method() const;
	double get_precision() const;
	bool do_use_neighbour_nodes_optimization() const;
	bool is_single_problem_mode() const;
	bool do_use_sleep_mode() const;
};

#endif
