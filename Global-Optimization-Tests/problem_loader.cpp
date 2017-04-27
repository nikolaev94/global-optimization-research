
#include "problem_loader.h"


ProblemLoader::ProblemLoader() {}

ProblemLoader::~ProblemLoader() {}


void ProblemLoader::set_GKLS_Hard_parameters(unsigned dimension)
{
	switch (dimension)
	{
	case 2:
		GklsFunction::GKLS_global_dist = 0.9;
		GklsFunction::GKLS_global_radius = 0.1;
		break;
	case 3:
		GklsFunction::GKLS_global_dist = 0.9;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	case 4:
		GklsFunction::GKLS_global_dist = 0.9;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	case 5:
		GklsFunction::GKLS_global_dist = 0.66;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	default:
		break;
	}
}


void ProblemLoader::set_GKLS_Simple_parameters(unsigned dimension)
{
	switch (dimension)
	{
	case 2:
		GklsFunction::GKLS_global_dist = 0.9;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	case 3:
		GklsFunction::GKLS_global_dist = 0.66;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	case 4:
		GklsFunction::GKLS_global_dist = 0.66;
		GklsFunction::GKLS_global_radius = 0.2;
		break;
	case 5:
		GklsFunction::GKLS_global_dist = 0.66;
		GklsFunction::GKLS_global_radius = 0.3;
		break;
	default:
		break;
	}
}

int ProblemLoader::set_GKLS_class_parameters(const UserParameters& user_parameter)
{
	GklsFunction::GKLS_dim = user_parameter.get_dimensions();
	GklsFunction::GKLS_set_default();

	switch (user_parameter.get_function_class())
	{
	case FunctionClass::GKLS_SIMPLE:
		set_GKLS_Simple_parameters(user_parameter.get_dimensions());
		break;
	case FunctionClass::GKLS_HARD:
		set_GKLS_Hard_parameters(user_parameter.get_dimensions());
		break;
	default:
		break;
	}

	return GklsFunction::GKLS_parameters_check();
}


void ProblemLoader::load_GKLS_series(const UserParameters& user_parameters,
	Solver::problem_list& problems)
{
	if (set_GKLS_class_parameters(user_parameters))
	{
		std::cerr << "Failed to set GKSL parameters" << std::endl;

		exit(EXIT_FAILURE);
	}


	for (std::size_t i = 1; i <= 100; i++)
	{
		GklsFunction::GklsFunctionPtr gkls_function(new GklsFunction(GklsFunction::D2, i));
		problems.push_back(OptProblem::OptProblemPtr(new GklsProblem(gkls_function)));
	}
}


void ProblemLoader::load_vagris_series(Solver::problem_list& problems)
{
	for (size_t i = 1; i <= 100; i++)
	{
		GrishaginFunction::GrishaginFunctionPtr vagris_function(new GrishaginFunction(i));
		problems.push_back(OptProblem::OptProblemPtr(new GrishaginProblem(vagris_function)));
	}
}


void ProblemLoader::print_GKLS_class_parameters()
{
	std::cout << "Dim: " << GklsFunction::GKLS_dim
		<< " Global dist: " << GklsFunction::GKLS_global_dist
		<< " Global radius: " << GklsFunction::GKLS_global_radius << std::endl;
}

std::string ProblemLoader::get_output_file_prefix(const UserParameters& user_parameters)
{
	std::string output_file_prefix;
	switch (user_parameters.get_function_class())
	{
	case FunctionClass::GKLS_HARD:
		output_file_prefix += "gkls_hard_";
		break;
	case FunctionClass::GKLS_SIMPLE:
		output_file_prefix += "gkls_simple_";
		break;
	case FunctionClass::VAGRIS:
		output_file_prefix += "vagris_";
		break;
	default:
		break;
	}

	output_file_prefix += std::string("dim")
		+ std::to_string(user_parameters.get_dimensions()) + "_";

	switch (user_parameters.get_solving_method())
	{
	case Solver::SolvingMethod::DYNAMIC:
		output_file_prefix += "dynamic_";
		break;
	case Solver::SolvingMethod::SIMULTANEOUS:
		output_file_prefix += "simult_";
		break;
	case Solver::SolvingMethod::SEQUENTIAL:
		output_file_prefix += "sequential_";
		break;
	default:
		break;
	}

	return output_file_prefix;
}


void ProblemLoader::load_series(const UserParameters& user_parameters,
	Solver::problem_list& problems)
{
	switch (user_parameters.get_function_class())
	{
	case FunctionClass::GKLS_SIMPLE:
	case FunctionClass::GKLS_HARD:
		load_GKLS_series(user_parameters, problems);
		print_GKLS_class_parameters();
		break;
	case FunctionClass::VAGRIS:
		load_vagris_series(problems);
		break;
	default:
		break;
	}
}

