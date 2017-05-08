
#include "problem_loader.h"


gkls::GKLSClass ProblemLoader::map_GKLS_function_class(const FunctionClass& function_class)
{
	switch (function_class)
	{
	case GKLS_SIMPLE:
		return gkls::GKLSClass::Simple;
	case GKLS_HARD:
		return gkls::GKLSClass::Hard;
	default:
		return gkls::GKLSClass::Simple;
	}
}


void ProblemLoader::load_GKLS_series(const UserParameters& user_parameters,
	Solver::problem_list& problems)
{

	auto series_size = user_parameters.get_series_size();

	auto selected_problem_no = user_parameters.get_selected_problem_no();

	auto gkls_function_class =
		map_GKLS_function_class(user_parameters.get_function_class());

	auto num_constrains = user_parameters.get_num_constrains();

	if (!selected_problem_no)
	{
		for (unsigned int no = 1; no <= series_size; no++)
		{
			GKLSProblem::GKLSFunctionPtr objective(new gkls::GKLSFunction());

			objective->SetFunctionClass(
				gkls_function_class, user_parameters.get_dimension());

			objective->SetType(GKLSProblem::DEFAULT_GKLS_FUNCTION_TYPE);

			objective->SetFunctionNumber(no);

			std::vector<GKLSProblem::GKLSFunctionPtr> constrains;

			for (unsigned int constraint_no = 0; constraint_no < num_constrains;)
			{
				GKLSProblem::GKLSFunctionPtr constraint_function(new gkls::GKLSFunction());

				constraint_function->SetFunctionClass(gkls_function_class,
					user_parameters.get_dimension());

				constraint_function->SetType(GKLSProblem::DEFAULT_GKLS_FUNCTION_TYPE);

				++constraint_no;

				constraint_function->SetFunctionNumber(no + constraint_no);

				constrains.push_back(constraint_function);
			}

			problems.push_back(OptProblem::OptProblemPtr(new GKLSProblem(objective,
				constrains)));
		}
	}
	else
	{
		GKLSProblem::GKLSFunctionPtr objective(new gkls::GKLSFunction());

		objective->SetFunctionClass(
			gkls_function_class, user_parameters.get_dimension());

		objective->SetType(GKLSProblem::DEFAULT_GKLS_FUNCTION_TYPE);

		objective->SetFunctionNumber(selected_problem_no);


		std::vector<GKLSProblem::GKLSFunctionPtr> constrains;


		for (unsigned int constraint_no = 0; constraint_no < num_constrains;)
		{
			GKLSProblem::GKLSFunctionPtr constraint_function(new gkls::GKLSFunction());

			constraint_function->SetFunctionClass(gkls_function_class,
				user_parameters.get_dimension());

			constraint_function->SetType(GKLSProblem::DEFAULT_GKLS_FUNCTION_TYPE);

			++constraint_no;

			constraint_function->SetFunctionNumber(selected_problem_no + constraint_no);

			constrains.push_back(constraint_function);
		}

		problems.push_back(OptProblem::OptProblemPtr(new GKLSProblem(objective, constrains)));
	}
}


void ProblemLoader::load_vagris_series(const UserParameters& user_parameters,
	Solver::problem_list& problems)
{
	auto series_size = user_parameters.get_series_size();

	auto selected_problem_no = user_parameters.get_selected_problem_no();

	if (!selected_problem_no)
	{
		for (unsigned int no = 1; no <= series_size; no++)
		{
			GrishaginFunction::GrishaginFunctionPtr vagris_function(
				new GrishaginFunction(no));
			problems.push_back(OptProblem::OptProblemPtr(
				new GrishaginProblem(vagris_function)));
		}
	}
	else
	{
		GrishaginFunction::GrishaginFunctionPtr vagris_function(
			new GrishaginFunction(selected_problem_no));
		problems.push_back(OptProblem::OptProblemPtr(
			new GrishaginProblem(vagris_function)));
	}
}


//void ProblemLoader::print_GKLS_class_parameters()
//{
//	std::cout << "Dim: " << GklsFunction::GKLS_dim
//		<< " Global dist: " << GklsFunction::GKLS_global_dist
//		<< " Global radius: " << GklsFunction::GKLS_global_radius << std::endl;
//}


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
		+ std::to_string(user_parameters.get_dimension()) + "_";

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
		break;
	case FunctionClass::VAGRIS:
		load_vagris_series(user_parameters, problems);
		break;
	default:
		break;
	}
}

