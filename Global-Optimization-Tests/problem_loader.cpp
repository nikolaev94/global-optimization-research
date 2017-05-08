
#include "problem_loader.h"


ProblemLoader::ProblemLoader() {}

ProblemLoader::~ProblemLoader() {}


//void ProblemLoader::set_GKLS_Hard_parameters(unsigned dimension)
//{
//	switch (dimension)
//	{
//	case 2:
//		GklsFunction::GKLS_global_dist = 0.9;
//		GklsFunction::GKLS_global_radius = 0.1;
//		break;
//	case 3:
//		GklsFunction::GKLS_global_dist = 0.9;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	case 4:
//		GklsFunction::GKLS_global_dist = 0.9;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	case 5:
//		GklsFunction::GKLS_global_dist = 0.66;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	default:
//		break;
//	}
//}


//void ProblemLoader::set_GKLS_Simple_parameters(unsigned dimension)
//{
//	switch (dimension)
//	{
//	case 2:
//		GklsFunction::GKLS_global_dist = 0.9;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	case 3:
//		GklsFunction::GKLS_global_dist = 0.66;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	case 4:
//		GklsFunction::GKLS_global_dist = 0.66;
//		GklsFunction::GKLS_global_radius = 0.2;
//		break;
//	case 5:
//		GklsFunction::GKLS_global_dist = 0.66;
//		GklsFunction::GKLS_global_radius = 0.3;
//		break;
//	default:
//		break;
//	}
//}

//int ProblemLoader::set_GKLS_class_parameters(const UserParameters& user_parameter)
//{
//	GklsFunction::GKLS_dim = user_parameter.get_dimensions();
//	GklsFunction::GKLS_set_default();
//
//	switch (user_parameter.get_function_class())
//	{
//	case FunctionClass::GKLS_SIMPLE:
//		set_GKLS_Simple_parameters(user_parameter.get_dimensions());
//		break;
//	case FunctionClass::GKLS_HARD:
//		set_GKLS_Hard_parameters(user_parameter.get_dimensions());
//		break;
//	default:
//		break;
//	}
//
//	return GklsFunction::GKLS_parameters_check();
//}


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

	if (!selected_problem_no)
	{
		for (unsigned int no = 1; no <= series_size; no++)
		{
			GKLSProblem::GKLSFunctionPtr objective(new gkls::GKLSFunction());

			objective->SetFunctionClass(
				gkls_function_class, user_parameters.get_dimension());

			objective->SetType(gkls::GKLSFuncionType::TD2);

			objective->SetFunctionNumber(no);


			problems.push_back(OptProblem::OptProblemPtr(new GKLSProblem(objective)));

			/*GklsFunction::GklsFunctionPtr gkls_function(
				new GklsFunction(GklsFunction::D2, no));
			problems.push_back(OptProblem::OptProblemPtr(new GklsProblem(gkls_function)));*/
		}
	}
	else
	{
		GKLSProblem::GKLSFunctionPtr objective(new gkls::GKLSFunction());

		objective->SetFunctionClass(
			gkls_function_class, user_parameters.get_dimension());

		objective->SetType(gkls::GKLSFuncionType::TD);

		objective->SetFunctionNumber(1);


		std::vector<GKLSProblem::GKLSFunctionPtr> my_list;




		GKLSProblem::GKLSFunctionPtr constraint(new gkls::GKLSFunction());

		constraint->SetFunctionClass(gkls_function_class, user_parameters.get_dimension());

		constraint->SetType(gkls::GKLSFuncionType::TD);

		constraint->SetFunctionNumber(2);

		my_list.push_back(constraint);


		problems.push_back(OptProblem::OptProblemPtr(new GKLSProblem(objective, my_list)));

		/*GklsFunction::GklsFunctionPtr gkls_function(
			new GklsFunction(GklsFunction::D2, selected_problem_no));
		problems.push_back(OptProblem::OptProblemPtr(new GklsProblem(gkls_function)));*/
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

