
#include "Solver.h"

#include <iostream>
#include <fstream>

Solver::Solver() // : scheduler(1)
{
}

Solver::Solver(const Input& _input, const problem_list& _problems)
	: scheduler(tbb::task_scheduler_init::deferred), input(_input), problems(_problems)
{
	init_tbb(input.num_threads);
}

Solver::Solver(const Input& _input) : scheduler(tbb::task_scheduler_init::deferred), input(_input)
{
	init_tbb(input.num_threads);
}


Solver::~Solver()
{
}


void Solver::Output::dump_results_to_file(const std::string& filename)
{
	std::ofstream ofs(filename);

	ofs << "xmin" << ';' << "zmin" << ';' << "error" << ';'
		<< "Problem trials" << ';' << "Iterations" << ';'
		<< "Total trials" << ';' << std::endl;

	for (const auto& res : this->results)
	{
		ofs << res.xmin << ';' << res.zmin << ';' << res.error << ';'
			<< res.trials << ';' << res.iterations << ';'
			<< res.total_trials << ';' << std::endl;
	}

	ofs << std::endl << "Precise" << ';' << "Parameter" << ';'
		<< "Number of workers" << ';' << "Total trials" << ';'
		<< "Total iterations" << ';' << "Elapsed time" << ';' << std::endl;

	ofs << this->original_input.method_eps << ';' << this->original_input.method_param << ';'
		<< this->original_input.num_threads << ';'
		<< MethodData::global_trials_count << ';'
		<< MethodData::global_iterations_count << ';' << this->elapsed_time << ';' << std::endl;

	ofs.close();
}


void Solver::Output::dump_error_metrics_by_trials_to_file(const std::string& filename)
{
	std::ofstream ofs(filename);

	/*
	ofs << "Number of trials" << ';' << "Average error" << ';';

	for (const auto& avg : this->metrics.average_error_by_trials)
	{
		ofs << avg.first << ';' << avg.second << ';';
	}

	ofs << ';' << "Number of trials" << ';' << "Max error" << ';' << std::endl;

	for (const auto& maximum : this->metrics.max_error_by_trials)
	{
		ofs << ';' << maximum.first << '; ' << maximum.second << '; ' << std::endl;
	}
	*/

	ofs << "Num trials" << ';' << "Average error" << ';' << "Max error" << std::endl;

	for (const auto& error: this->metrics.errors_by_trials)
	{
		ofs << error.first << ';' << error.second.first << ';' << error.second.second << std::endl;
	}

	ofs.close();
}


void Solver::Output::dump_solved_problem_portion_by_trials_to_file(const std::string& filename)
{
	std::ofstream ofs(filename);

	ofs << "Num trials" << ';' << "Solved problem" << std::endl;

	for (const auto& portion: this->metrics.solved_problems_portion_by_trials)
	{
		ofs << portion.first << ';' << portion.second << std::endl;
	}

	ofs.close();
}


void Solver::init_tbb(int number_threads)
{
	this->scheduler.initialize(number_threads);
}


void Solver::Trial::perform_trial(problem_iterator problem)
{
	OptProblem* opt_problem = problem->get();
	for (size_t i = 0; i < opt_problem->getConstraintsNumber(); i++)
	{
		double func_value = opt_problem->getContraintValue(i, this->x);
		if (func_value > 0)
		{
			this->z = func_value;
			this->nu = i + 1;
			this->admissible = false;
			return;
		}
	}
	this->z = opt_problem->getObjectiveValue(this->x);
	this->nu = opt_problem->getConstraintsNumber() + 1;
	this->admissible = true;
}

bool Solver::TrialSubset::use_neighbour_nodes = true;

double Solver::TrialSubset::get_subset_max_difference()
{
	if (TrialSubset::use_neighbour_nodes && this->subset.size() > this->SUBSET_SIZE_LIMIT)
	{
		return calc_max_difference_neighbours();
	}
	else
	{
		return calc_max_difference_all_trials();
	}
}


double Solver::TrialSubset::calc_max_difference_all_trials()
{
	double maxDiff = -DBL_MAX;
	for (auto i = this->subset.begin(); i != --this->subset.end(); i++)
	{
		for (auto j = std::next(i); j != this->subset.end(); j++)
		{
			if (fabs(j->z - i->z) / (j->x - i->x) > maxDiff)
			{
				maxDiff = fabs(j->z - i->z) / (j->x - i->x);
			}

		}
	}
	if (maxDiff < DBL_EPSILON)
	{
		maxDiff = 1.0;
	}
	return maxDiff;
}


double Solver::TrialSubset::calc_max_difference_neighbours()
{
	double maxDiff = -DBL_MAX;
	for (auto i = this->subset.begin(); i != --this->subset.end(); i++)
	{
		auto j = std::next(i);
		if (fabs(j->z - i->z) / (j->x - i->x) > maxDiff)
		{
			maxDiff = fabs(j->z - i->z) / (j->x - i->x);
		}
	}
	if (maxDiff < DBL_EPSILON)
	{
		maxDiff = 1.0;
	}
	return maxDiff;
}


void Solver::TrialSubset::calc_subset_lower_lip_const()
{
	if (this->subset.size() < 2) {
		this->lip_const = 1.0;
	} else {
		this->lip_const = this->get_subset_max_difference();
	}
}


void Solver::TrialSubset::calc_subset_min_estimate()
{
	double zmin = DBL_MAX;
	for (const auto& trial : this->subset)
	{
		if (trial.z < zmin)
		{
			zmin = trial.z;
		}
	}
	this->min_estimator = zmin;
}


/*
void Solver::get_trial_subsets_by_index(const std::set<Trial>& trials, trial_subsets& subsets)
{
	
	for each (auto index in subsets)
	{
		for each (auto trial in subsets[index.first])
		{
			std::cout << index.first << ' ' << trial.x << ' ' << trial.z << ' ' << trial.nu << std::endl;
		}
	}

	std::for_each(trials.begin(), trials.end(), [&subsets](const Trial& trial) {
		std::cout << trial.x << ' ' << trial.z << ' ' << trial.nu << std::endl;
		subsets[trial.nu].insert(trial);
	});

	for each (auto index in subsets)
	{
		for each (auto trial in subsets[index.first].subset)
		{
			std::cout << index.first << ' ' << trial.x << ' ' << trial.z << ' ' << trial.nu << std::endl;
		}
	}
}
*/


double Solver::FunctionStatsInfo::initial_parameter = 10.0;


void Solver::FunctionStatsInfo::update()
{
	this->calc_count++;
	if (this->calc_count > this->CALC_COUNT_LIMIT) {
		this->parameter = this->DEFAULT_PARAMETER;
	}
}

/*
size_t Solver::MethodDataHasher::operator()(const MethodData& method_data) const
{
	return method_data.problem_number;
}
*/


size_t Solver::ProblemIteratorHasher::operator() (const problem_iterator& iterator) const
{
	return (size_t)&(*iterator);
}


Solver::MethodData::MethodData(problem_iterator _problem) : problem(_problem)
{
	this->trials.emplace(this->input.left, 0.0, 0);
	this->trials.emplace(this->input.right, 0.0, 0);

	Trial mid_node((this->input.left + this->input.right) / 2.0);
	mid_node.perform_trial(this->problem);
	this->trials.insert(mid_node);

	this->sln_estimator.error = fabs(input.right - input.left);

	this->init_function_stats();

	this->init_trial_subsets();

	this->update_trial_subsets();
}

/*
bool Solver::MethodData::operator==(const MethodData& cmp) const
{
	return this->problem_number == cmp.problem_number;
}
*/

/*
size_t Solver::MethodData::hash() const
{
	return this->problem_number;
}
*/

void Solver::MethodData::init_function_stats()
{
	OptProblem* opt_problem = this->problem->get();
	for (size_t i = 0; i <= opt_problem->getConstraintsNumber() + 1; i++)
	{
		this->function_stats.emplace(i, this->input.method_param);
	}
}


void Solver::MethodData::init_trial_subsets()
{
	std::for_each(trials.begin(), trials.end(), [&](const Trial& trial) {
		subsets[trial.nu].insert(trial);
	});
}


void Solver::MethodData::update_trial_subsets()
{
	this->update_method_parameters();

	this->calc_lower_lip_const();

	this->calc_min_estimators();
}

void Solver::MethodData::calc_lower_lip_const()
{
	for (auto i = subsets.begin(); i != --subsets.end(); i++)
	{
		i->second.min_estimator = 0.0;
	}

	auto max_nu = subsets.rbegin();

	max_nu->second.calc_subset_min_estimate();
}


void Solver::MethodData::calc_min_estimators()
{
	for (auto& subset_index: subsets)
	{
		subset_index.second.calc_subset_lower_lip_const();
	}
}


void Solver::MethodData::update_method_parameters()
{
	for (auto& subset_index : subsets)
	{
		subset_index.second.method_parameter = function_stats.at(subset_index.first).parameter;
	}
}


void Solver::MethodData::add_trial(const Trial& another_trial)
{
	this->trials.insert(another_trial);

	this->trials_count++;

	this->update_stats(another_trial);

	this->update_solution(another_trial);

	this->add_trial_to_subset(another_trial);
}


void Solver::MethodData::add_trial_to_subset(const Trial& another_trial)
{
	this->subsets[another_trial.nu].insert(another_trial);
}


void Solver::MethodData::update_stats(const Trial& another_trial)
{
	for (auto& function_index : function_stats)
	{
		if (function_index.first <= another_trial.nu)
		{
			function_index.second.update();
		}
	}
}


void Solver::MethodData::update_solution(const Trial& another_trial)
{
	if (another_trial.admissible)
	{
		OptProblem* opt_problem = this->problem->get();
		double error = opt_problem->getReferenceMinError(another_trial.x);

		if (error < this->sln_estimator.error)
		{
			this->sln_estimator.error = error;
			this->sln_estimator.xmin = another_trial.x;
			this->sln_estimator.zmin = another_trial.z;

			if (error < this->input.method_eps)
			{
				// this->total_trials_count = MethodData::global_trials_count;

				this->total_trials_count = MethodData::global_trials_count;
				this->total_iterations_count = MethodData::global_iterations_count;
				this->method_finished = true;
			}
		}
	}
}


void Solver::MethodData::construct_segment_set(std::multiset<Interval>& segments) const
{
	for (auto i = trials.begin(); i != --trials.end(); i++)
	{
		auto next = std::next(i);

		Interval curr_segment(problem, *i, *next);

		curr_segment.calc_characteristic(subsets);

		segments.insert(curr_segment);
	}
}


bool Solver::MethodData::is_finished() const
{
	return this->method_finished;
}


Solver::Input Solver::MethodData::input;

// unsigned Solver::MethodData::global_trials_count = 0;

unsigned Solver::MethodData::global_iterations_count = 0;

tbb::atomic<unsigned> Solver::MethodData::global_trials_count;

void Solver::MethodData::set_method_input(const Input& user_input)
{
	Solver::MethodData::input = user_input;
}


unsigned Solver::MethodData::get_num_trials() const
{
	return this->trials_count;
}


double Solver::MethodData::get_error_value() const
{
	return this->sln_estimator.error;
}


/*
void Solver::IntervalContainer::split_interval()
{

}
*/

void Solver::Interval::calc_characteristic(const trial_subsets& subsets)
{
	double delta = this->get_interval_length();
	if (this->left_node.nu == this->right_node.nu)
	{
		double min_estimate_nu = subsets.at(this->left_node.nu).min_estimator;

		this->greater_nu_method_param = subsets.at(this->left_node.nu).method_parameter;
		this->greater_nu_lipconst = subsets.at(this->left_node.nu).lip_const;

		this->charact = delta + pow(this->right_node.z - this->left_node.z, 2.0)
			/ (pow(this->greater_nu_method_param, 2.0) * pow(this->greater_nu_lipconst, 2.0) * delta)
			- 2.0 * (this->right_node.x + this->left_node.x - 2.0 * min_estimate_nu)
			/ (this->greater_nu_method_param * this->greater_nu_lipconst);
	}
	else if (this->right_node.nu > this->left_node.nu)
	{
		double min_estimate_nu = subsets.at(this->right_node.nu).min_estimator;

		this->greater_nu_lipconst = subsets.at(this->right_node.nu).lip_const;
		this->greater_nu_method_param = subsets.at(this->right_node.nu).method_parameter;

		this->charact = 2.0 * delta - 4.0 * (this->right_node.z - min_estimate_nu)
			/ (this->greater_nu_method_param * this->greater_nu_lipconst);
	}
	else
	{
		double min_estimate_nu = subsets.at(this->left_node.nu).min_estimator;

		this->greater_nu_lipconst = subsets.at(this->left_node.nu).lip_const;
		this->greater_nu_method_param = subsets.at(this->left_node.nu).method_parameter;

		this->charact = 2.0 * delta - 4.0 * (this->left_node.z - min_estimate_nu)
			/ (this->greater_nu_method_param * this->greater_nu_lipconst);
	}
}


double Solver::Interval::get_interval_length()
{
	OptProblem* opt_problem = problem->get();
	unsigned problem_dim = opt_problem->getDimention();
	double delta = 0.0;
	if (problem_dim > 1)
	{
		if (problem_dim == 2)
		{
			delta = sqrt(fabs(right_node.x - left_node.x));
		}
		else
		{
			delta = pow(fabs(right_node.x - left_node.x), 1.0 / problem_dim);
		}
	}
	else
	{
		delta = fabs(right_node.x - left_node.x);
	}
	return delta;
}


double Solver::sgn(double arg)
{
	if (arg < DBL_EPSILON && arg > -(DBL_EPSILON))
	{
		return 0.0;
	}
	if (arg > 0.0)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}


double Solver::Interval::get_new_point() const
{
	double x = 0.0;
	if (right_node.nu != left_node.nu)
	{
		x = (right_node.x + left_node.x) / 2.0;
	}
	else
	{
		OptProblem* opt_problem = problem->get();
		unsigned problem_dim = opt_problem->getDimention();
		if (problem_dim > 1)
		{
			x = (right_node.x + left_node.x) / 2.0 - sgn(right_node.z - left_node.z)
				* pow(fabs(right_node.z - left_node.z) / greater_nu_lipconst, 1.0 * problem_dim)
				/ (2.0 * greater_nu_method_param);
		}
		else
		{
			x = (right_node.x + left_node.x) / 2.0
				- (right_node.z - left_node.z) / (2.0 * greater_nu_lipconst * greater_nu_method_param);
		}
	}
	return x;
}


Solver::Trial Solver::Interval::create_new_trial() const
{
	Trial another_trial(this->get_new_point());
	another_trial.perform_trial(this->problem);
	return another_trial;
}

tbb::mutex Solver::SimultaneousMethodDataContainer::mutex;

void Solver::SimultaneousMethodDataContainer::construct_and_merge_sets(std::multiset<Interval>& merged_set)
{
	for (const auto& method_data_index : this->problem_series_container)
	{
		auto& problem_data = method_data_index.second;
		if (problem_data.is_finished())
		{
			continue;
		}

		std::multiset<Interval> problem_segments;
		problem_data.construct_segment_set(problem_segments);
		merged_set.insert(problem_segments.begin(), problem_segments.end());
	}
}


void Solver::SimultaneousMethodDataContainer::add_new_trial(problem_iterator target_problem, Trial new_trial)
{
	auto target_data = this->problem_series_container.find(ProblemIterator(target_problem));

	SimultaneousMethodDataContainer::mutex.lock();
	target_data->second.add_trial(new_trial);
	SimultaneousMethodDataContainer::mutex.unlock();
}


void Solver::SimultaneousMethodDataContainer::complete_iteration()
{
	MethodData::global_iterations_count++;

	for (auto& method_data_index : this->problem_series_container)
	{
		auto& problem_data = method_data_index.second;
		// std::cout << problem_data.trials_count << std::endl;
		if (!problem_data.is_finished())
		{
			problem_data.update_trial_subsets();
			// problem_data.iterations_count++;
		}
	}
}


void Solver::SimultaneousMethodDataContainer::parallel_perform_iteration(const std::multiset<Interval>& segments_to_consider)
{
	tbb::parallel_do(segments_to_consider.begin(), segments_to_consider.end(),
		[this](const Interval& segment) -> void
	{
		Trial another_trial = segment.create_new_trial();

		/*
		* Updating global trial counter
		*/
		++MethodData::global_trials_count;

		this->add_new_trial(segment.problem, another_trial);
	});

	/*
	for (auto segment_iter = all_segments.begin(); segment_iter != segments_end; ++segment_iter)
	{
		Trial another_trial = segment_iter->create_new_trial();

		/*
		* Updating global trial counter
		*
		MethodData::global_trials_count++;

		problems_method_data.add_new_trial(segment_iter->problem, another_trial);
	}*/
}


void Solver::MethodDataContainer::update_metrics(MetricsContainer& metrics)
{
	update_errors(metrics.errors_by_trials);

	update_portion(metrics.solved_problems_portion_by_trials);
}


void Solver::MethodDataContainer::update_portion(portion_vector& solved_problems_portion)
{
	unsigned num_problems_finished = 0;

	for (const auto& method_data_index : this->problem_series_container)
	{
		if (method_data_index.second.is_finished())
		{
			++num_problems_finished;
		}
	}

	double portion = 1.0 * num_problems_finished / this->problem_series_container.size();
	solved_problems_portion.emplace_back(MethodData::global_trials_count, portion);
}

void Solver::MethodDataContainer::update_errors(errors_vector& errors_by_trials)
{
	double average_error = 0.0;
	double max_error = -(DBL_MAX);

	// unsigned num_trials = 0;
	for (const auto& method_data_index : this->problem_series_container)
	{
		// num_trials += method_data_index.second.get_num_trials();

		double method_error = method_data_index.second.get_error_value();

		average_error += method_error;
		if (method_error > max_error)
		{
			max_error = method_error;
		}
	}
	average_error /= this->problem_series_container.size();
	errors_by_trials.push_back(std::make_pair(MethodData::global_trials_count, std::make_pair(average_error, max_error)));
	// errors_by_trials.push_back(std::make_pair(num_trials, std::make_pair(average_error, max_error)));

	// metrics.errors_by_trials.push_back(std::make_pair(num_trials, std::make_pair(average_error, max_error)));
}


bool Solver::MethodDataContainer::is_all_finished()
{
	bool done = true;
	for (const auto& method_data_index : this->problem_series_container)
	{
		done &= method_data_index.second.is_finished();
	}

	return done;
}


void Solver::SimultaneousMethodDataContainer::add_problem(problem_iterator problem)
{
	this->problem_series_container.emplace(problem, problem);
}


Solver::MetricsContainer::MetricsContainer(Input input)
{
	double dist = fabs(input.right - input.left);
	this->errors_by_trials.push_back(std::make_pair(0, std::make_pair(dist, dist)));
	// this->errors_by_trials.emplace_back(0, dist, dist);

	/*
	this->average_error_by_trials.emplace_back(0, dist);
	this->max_error_by_trials.emplace_back(0, dist);
	*/
	

	// this->max_error_by_trials.emplace_back(0, 0.5 * distance);
}


void Solver::MethodDataContainer::dump_solving_results(std::list<ProblemSolvingResult>& results)
{
	for (const auto& method_data_index : this->problem_series_container)
	{
		results.emplace_back(method_data_index.second);
	}
}


void Solver::DynamicMethodDataContainer::enqueue_problem(problem_iterator problem)
{
	this->problem_series_container.emplace(problem, problem);


	this->problem_queue.emplace(problem);

	//this->problem_queue.push(problem);
}


/*
bool Solver::DynamicMethodDataContainer::queue_is_empty()
{
	return this->problem_queue.empty();
}
*/


void Solver::DynamicMethodDataContainer::init_workers(unsigned int num_threads)
{
	for (size_t i = 0; i < num_threads; i++)
	{
		take_problem_from_queue();
	}
}


void Solver::DynamicMethodDataContainer::take_problem_from_queue()
{
	if (!this->problem_queue.empty())
	{
		auto another_problem_to_solve = problem_queue.front();
		problem_queue.pop();
		//auto target_problem_data = this->problem_series_container.find(another_problem_to_solve);
		auto target_problem_data = this->problem_series_container.find(ProblemIterator(another_problem_to_solve));

		active_solving_problems.push_front(std::ref(target_problem_data->second));
		// active_solving_problems.push_front(target_problem_data->second);
	}
}

void Solver::DynamicMethodDataContainer::perform_iteration()
{
	/*
	cilk::reducer<cilk::op_add<unsigned>> trials_sum(0);

	cilk_for(auto solving_problem_reference : this->active_solving_problems)
	{
		auto& solving_problem = solving_problem_reference.get();

		std::multiset<Interval> interval_set;
		solving_problem.construct_segment_set(interval_set);
		auto best_segment = interval_set.begin();
		Trial another_trial = best_segment->create_new_trial();

		/*
		* Updating global trial counter
		
		MethodData::global_trials_count++;

		solving_problem.add_trial(another_trial);
	}
	*/

	for (auto solving_problem_reference : this->active_solving_problems)
	{
		auto& solving_problem = solving_problem_reference.get();

		std::multiset<Interval> interval_set;
		solving_problem.construct_segment_set(interval_set);
		auto best_segment = interval_set.begin();
		Trial another_trial = best_segment->create_new_trial();

		/*
		 * Updating global trial counter
		 */
		MethodData::global_trials_count++;

		solving_problem.add_trial(another_trial);
	}
}


void Solver::DynamicMethodDataContainer::parallel_perform_iteration()
{
	tbb::parallel_do(active_solving_problems.begin(), active_solving_problems.end(),
		[](std::reference_wrapper<MethodData> solving_problem_reference) -> void
	{
		auto& solving_problem = solving_problem_reference.get();

		std::multiset<Interval> interval_set;
		solving_problem.construct_segment_set(interval_set);
		auto best_segment = interval_set.begin();
		Trial another_trial = best_segment->create_new_trial();

		/*
		* Updating global trial counter
		*/
		++MethodData::global_trials_count;
		//MethodData::global_trials_count++;

		solving_problem.add_trial(another_trial);
	});
}

void Solver::DynamicMethodDataContainer::complete_iteration()
{
	MethodData::global_iterations_count++;

	/*
	 * note: there is no increment in the loop construct
	 */
	for (auto ref_iterator = active_solving_problems.begin(); ref_iterator != active_solving_problems.end();)
	{
		auto& solving_problem = (*ref_iterator).get();
		if (solving_problem.is_finished())
		{
			/*
			 * removing problem from active list
			 */
			ref_iterator = active_solving_problems.erase(ref_iterator);
			/*
			 * take another one
			 */
			take_problem_from_queue();
		}
		else
		{
			solving_problem.update_trial_subsets();

			++ref_iterator;
		}
	}
}

Solver::problem_list::const_iterator Solver::ProblemIterator::problem_list_begin;


bool Solver::ProblemIterator::operator<(const ProblemIterator& cmp) const
{
	return (std::distance(problem_list_begin, this->problem_ptr)
		< std::distance(problem_list_begin, cmp.problem_ptr));
}

void Solver::run_simultaneous_search(Output& out)
{
	Solver::MethodData::set_method_input(input);
	tbb::tick_count start_time, finish_time;

	start_time = tbb::tick_count::now();
	SimultaneousMethodDataContainer problems_method_data;
	for (auto i = problems.begin(); i != problems.end(); ++i)
	{
		problems_method_data.add_problem(i);
	}

	do
	{
		std::multiset<Interval> all_segments;
		problems_method_data.construct_and_merge_sets(all_segments);

		auto segments_end = all_segments.begin();
		/*
		 * Need to check if end is reached
		 */
		if (all_segments.size() < input.num_threads)
		{
			segments_end = all_segments.end();
		}
		else
		{
			std::advance(segments_end, input.num_threads);
		}

		problems_method_data.parallel_perform_iteration(std::multiset<Interval>(all_segments.begin(), segments_end));

		/*for (auto segment_iter = all_segments.begin(); segment_iter != segments_end; ++segment_iter)
		{
			Trial another_trial = segment_iter->create_new_trial();

			/*
			 * Updating global trial counter
			 *
			MethodData::global_trials_count++;

			problems_method_data.add_new_trial(segment_iter->problem, another_trial);
		}*/

		problems_method_data.update_metrics(out.metrics);

		problems_method_data.complete_iteration();
	} while (!problems_method_data.is_all_finished());

	finish_time = tbb::tick_count::now();
	out.elapsed_time = (finish_time - start_time).seconds();
	problems_method_data.dump_solving_results(out.results);
}


void Solver::run_dynamic_search(Output& out)
{
	Solver::MethodData::set_method_input(this->input);
	tbb::tick_count start_time, finish_time;

	start_time = tbb::tick_count::now();
	DynamicMethodDataContainer problems_method_data;

	for (auto i = this->problems.begin(); i != this->problems.end(); ++i)
	{
		problems_method_data.enqueue_problem(i);
	}

	problems_method_data.init_workers(input.num_threads);

	do
	{
		problems_method_data.parallel_perform_iteration();
		// problems_method_data.perform_iteration();

		problems_method_data.update_metrics(out.metrics);

		problems_method_data.complete_iteration();
	} while (!problems_method_data.is_all_finished());

	finish_time = tbb::tick_count::now();
	out.elapsed_time = (finish_time - start_time).seconds();

	problems_method_data.dump_solving_results(out.results);
}


void Solver::run_sequential_search(Output& out)
{
	Solver::MethodData::set_method_input(this->input);
	tbb::tick_count start_time, finish_time;

	start_time = tbb::tick_count::now();

}

void Solver::run_solver(Output& out) {
	Solver::ProblemIterator::problem_list_begin = this->problems.begin();

	switch (this->input.solving_method)
	{
	case SIMULTANEOUS:
		run_simultaneous_search(out);
		break;
	case DYNAMIC:
		run_dynamic_search(out);
		break;
	case SEQUENTIAL:
		run_sequential_search(out);
		break;
	default:
		break;
	}
}

