
#include "user_input.h"


double UserParameters::LEFT_BOUND = 0.0f;

double UserParameters::RIGHT_BOUND = 1.0f;

void UserParameters::show_usage()
{
	std::cout << "Usage: Global-Optimization-Tests.exe [options]" << std::endl;
	std::cout << "-f --function-class <GKLS_Simple|GKLS_Hard|Vagris>" << std::endl;
	std::cout << "-s --strategy <Dynamic|Simultaneous|Sequential>" << std::endl;
	std::cout << "-d --dim <integer>: GKLS problems dimension" << std::endl;
	std::cout << "-p --parameter <float>: Method initial parameter" << std::endl;
	std::cout << "-e --eps <float>: Method precison" << std::endl;
	std::cout << "-w --workers <integer>: Number of workers" << std::endl;
	std::cout << "-n <integer>: Size of problem series" << std::endl;
	std::cout << "-i <integer>: Solve problem in demo mode" << std::endl;
	std::cout << "-c <integer>: Number of problem constraints" << std::endl;
	std::cout << "--use-all-nodes: Use all trials to calculate lower lip constant" << std::endl;
	std::cout << "" << std::endl;

	exit(EXIT_SUCCESS);
}

void UserParameters::parse_arguments_from_command_line(int argc, char* argv[])
{
	for (int i = 1; i < argc; i++)
	{
		std::string option(argv[i]);
		if (option == "-h" || option == "--help")
		{
			show_usage();
		}
		else if (option == "-f" || option == "--function-class")
		{
			if (i + 1 < argc)
			{
				std::string function_class_name(argv[++i]);

				if (function_class_name == "GKLS" || function_class_name == "GKLS_Simple")
				{
					this->function_class = FunctionClass::GKLS_SIMPLE;
				}
				else if (function_class_name == "GKLS_Hard")
				{
					this->function_class = FunctionClass::GKLS_HARD;
				}
				else if (function_class_name == "Vagris")
				{
					this->function_class = FunctionClass::VAGRIS;
				}
				else
				{
					std::cerr << option << ": argument "
						<< function_class_name << " is invalid." << std::endl;
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				std::cerr << option << ": option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-s" || option == "--strategy")
		{
			if (i + 1 < argc)
			{
				std::string solving_method_name(argv[++i]);

				if (solving_method_name == "dynamic" || solving_method_name == "Dynamic"
					|| solving_method_name == "dyn")
				{
					this->solving_method = Solver::DYNAMIC;
				}
				else if (solving_method_name == "simultaneous"
					|| solving_method_name == "Simultaneous"
					|| solving_method_name == "sim")
				{
					this->solving_method = Solver::SIMULTANEOUS;
				}
				else if (solving_method_name == "sequential"
					|| solving_method_name == "Sequential"
					|| solving_method_name == "seq")
				{
					this->solving_method = Solver::SEQUENTIAL;
				}
				else
				{
					std::cerr << option << ": argument "
						<< solving_method_name << " is invalid." << std::endl;
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-d" || option == "--dim")
		{
			if (i + 1 < argc)
			{
				int dim = atoi(argv[++i]);

				if (dim < 2)
				{
					std::cerr << option << ": argument "
						<< dim << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->dimensions = static_cast<unsigned int> (dim);
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-p" || option == "--parameter")
		{
			if (i + 1 < argc)
			{
				double method_parameter = atof(argv[++i]);
				if (method_parameter < 1.0 + DBL_EPSILON)
				{
					std::cerr << option << ": argument "
						<< method_parameter << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->method_parameter = method_parameter;
				}
			}
			else
			{
				std::cerr << option
					<< " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-e" || option == "--eps")
		{
			if (i + 1 < argc)
			{
				double eps = atof(argv[++i]);
				if (eps < DBL_EPSILON)
				{
					std::cerr << option
						<< ": argument " << eps << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->precision = eps;
				}
			}
			else
			{
				std::cerr << option
					<< " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-w" || option == "--workers")
		{
			if (i + 1 < argc)
			{
				int num_workers = atoi(argv[++i]);
				if (num_workers < 1)
				{
					std::cerr << option << ": argument "
						<< num_workers << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->num_workers = static_cast<unsigned> (num_workers);
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-n")
		{
			if (i + 1 < argc)
			{
				int series_size = atoi(argv[++i]);

				if (series_size < 1)
				{
					std::cerr << option << ": argument "
						<< series_size << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->series_size = static_cast<unsigned> (series_size);
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-i")
		{
			if (i + 1 < argc)
			{
				int problem_no = atoi(argv[++i]);

				if (problem_no < 1)
				{
					std::cerr << option << ": argument "
						<< problem_no << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->problem_no = static_cast<unsigned> (problem_no);

					this->single_problem_mode = true;
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "-c")
		{
			if (i + 1 < argc)
			{
				int num_constraints = atoi(argv[++i]);

				if (num_constraints < 1)
				{
					std::cerr << option << ": argument "
						<< problem_no << " is invalid" << std::endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					this->num_constraints = static_cast<unsigned> (num_constraints);
				}
			}
			else
			{
				std::cerr << option << " option requires one argument." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (option == "--use-all-nodes")
		{
			this->use_neighbour_nodes_optimization = false;
		}
		else if (option == "--dump-nodes")
		{
			this->dump_trial_nodes = true;
		}
		else
		{
			std::cerr << option << ": unknown option." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	if (function_class == FunctionClass::VAGRIS)
	{
		dimensions = 2;
	}
}

unsigned int UserParameters::get_dimension() const
{
	return this->dimensions;
}

unsigned int UserParameters::get_series_size() const
{
	return this->series_size;
}


unsigned int UserParameters::get_selected_problem_no() const
{
	return this->problem_no;
}


FunctionClass UserParameters::get_function_class() const
{
	return this->function_class;
}

unsigned int UserParameters::get_workers_number() const
{
	return this->num_workers;
}

double UserParameters::get_method_parameter() const
{
	return this->method_parameter;
}

Solver::SolvingMethod UserParameters::get_solving_method() const
{
	return this->solving_method;
}

double UserParameters::get_precision() const
{
	return this->precision;
}


unsigned int UserParameters::get_num_constrains() const
{
	return this->num_constraints;
}


bool UserParameters::do_use_neighbour_nodes_optimization() const
{
	return this->use_neighbour_nodes_optimization;
}


bool UserParameters::do_dump_trial_nodes() const
{
	return this->dump_trial_nodes;
}


bool UserParameters::is_single_problem_mode() const
{
	return this->single_problem_mode;
}
