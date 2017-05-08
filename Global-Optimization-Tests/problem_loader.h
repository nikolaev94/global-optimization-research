
#ifndef PROBLEM_LOADER_H
#define PROBLEM_LOADER_H

#include <vector>
#include <memory>

#include <grishagin_function.h>
#include <grishagin_problem.h>

//#include <gkls_function.h>
//#include <gkls_problem.h>

#include <gkls\gkls_function.hpp>
#include <GKLSProblem.h>



#include <solver.h>

#include "user_input.h"

class ProblemLoader
{
private:
	//void set_GKLS_Simple_parameters(unsigned dimension);
	//void set_GKLS_Hard_parameters(unsigned dimension);

	//int set_GKLS_class_parameters(const UserParameters& user_parameter);

	//void print_GKLS_class_parameters();

	gkls::GKLSClass map_GKLS_function_class(const FunctionClass&);

	void load_GKLS_series(const UserParameters&, Solver::problem_list&);

	void load_vagris_series(const UserParameters&, Solver::problem_list&);
public:
	ProblemLoader();
	~ProblemLoader();

	std::string get_output_file_prefix(const UserParameters&);

	void load_series(const UserParameters&, Solver::problem_list&);
};

#endif
