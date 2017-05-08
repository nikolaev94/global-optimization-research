
#ifndef PROBLEM_LOADER_H
#define PROBLEM_LOADER_H

#include <vector>
#include <memory>

#include <grishagin_function.h>
#include <grishagin_problem.h>

#include <gkls\gkls_function.hpp>
#include <GKLSProblem.h>

#include <solver.h>

#include "user_input.h"

class ProblemLoader
{
private:
	gkls::GKLSClass map_GKLS_function_class(const FunctionClass&);

	void load_GKLS_series(const UserParameters&, Solver::problem_list&);

	void load_vagris_series(const UserParameters&, Solver::problem_list&);
public:

	std::string get_output_file_prefix(const UserParameters&);

	void load_series(const UserParameters&, Solver::problem_list&);
};

#endif
