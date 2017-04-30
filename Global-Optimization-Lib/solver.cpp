
#include "Solver.h"

#include <iostream>
#include <fstream>

Solver::Solver(const Input& _input, const problem_list& _problems)
	: scheduler(tbb::task_scheduler_init::deferred), input(_input), problems(_problems)
{
	init_tbb(input.num_threads);
}

Solver::~Solver() {}


void Solver::Output::dump_results_to_file(const std::string& filename)
{
	std::ofstream ofstream(filename);

	ofstream << "xmin" << ';' << "zmin" << ';' << "error" << ';'
		<< "Problem trials" << ';' << "Iterations" << ';'
		<< "Total trials" << ';' << "Time" << ';' << std::endl;

	for (const auto& res : results)
	{
		ofstream << res.xmin << ';' << res.zmin << ';' << res.error << ';'
			<< res.trials << ';' << res.iterations << ';'
			<< res.total_trials << ';' << res.elapsed_time
			<< ';' << std::endl;
	}

	ofstream << std::endl << "Precise" << ';' << "Parameter" << ';'
		<< "Number of workers" << ';' << "Total trials" << ';'
		<< "Total iterations" << ';' << "Elapsed time" << ';' << std::endl;

	ofstream << this->original_input.method_eps << ';' << this->original_input.method_param << ';'
		<< this->original_input.num_threads << ';'
		<< MethodData::global_trials_count << ';'
		<< MethodData::global_iterations_count << ';' << this->elapsed_time << ';' << std::endl;

	ofstream.close();
}


void Solver::Output::dump_error_metrics_by_trials_to_file(const std::string& filename)
{
	std::ofstream ofstream(filename);

	ofstream << "Num trials" << ';' << "Average error" << ';' << "Max error" << std::endl;

	for (const auto& error: this->metrics.errors_by_trials)
	{
		ofstream << error.first << ';' << error.second.first <<
			';' << error.second.second << std::endl;
	}

	ofstream.close();
}


void Solver::Output::dump_solved_problem_portion_by_trials_to_file(const std::string& filename)
{
	std::ofstream ofstream(filename);

	ofstream << "Num trials" << ';' << "Solved problem" << std::endl;

	for (const auto& portion: this->metrics.solved_problems_portion_by_trials)
	{
		ofstream << portion.first << ';' << portion.second << std::endl;
	}

	ofstream.close();
}


void Solver::Output::add_problem_solving_result(const MethodData &solved_problem_data)
{
	results.emplace_back(solved_problem_data);
}


void Solver::init_tbb(int number_threads)
{
	this->scheduler.initialize(number_threads);
}


void Solver::Trial::perform_trial(problem_iterator problem)
{
	//OptProblem* opt_problem = problem->get();

	for (size_t i = 0; i < (*problem)->getConstraintsNumber(); i++)
	{
		double func_value = (*problem)->getContraintValue(i, x);
		if (func_value > 0)
		{
			this->z = func_value;
			this->nu = i + 1;
			this->admissible = false;
			return;
		}
	}
	this->z = (*problem)->getObjectiveValue(this->x);
	this->nu = (*problem)->getConstraintsNumber() + 1;
	this->admissible = true;
}

bool Solver::TrialSubset::USE_NEIGHBOUR_NODES = true; ///false; // 

void Solver::TrialSubset::calc_subset_min_estimate()
{
	double zmin = DBL_MAX;

	for (const auto& trial : subset)
	{
		if (trial.z < zmin)
		{
			zmin = trial.z;
		}
	}
	min_estimator = zmin;
}


double Solver::TrialSubset::get_subset_max_difference()
{
	if (TrialSubset::USE_NEIGHBOUR_NODES && subset.size() > SUBSET_SIZE_LIMIT)
	{
		return calc_max_difference_neighbours_1();
	}
	else
	{
		return calc_max_difference_all_trials();
	}
}


double Solver::TrialSubset::calc_max_difference_all_trials()
{
	double maxDiff = 0.0;

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


double Solver::TrialSubset::calc_max_difference_neighbours_1()
{
	double maxDiff = 0.0;
	for (auto subset = this->subset.begin(); subset != --this->subset.end(); subset++)
	{
		auto j = std::next(subset);
		if (fabs(j->z - subset->z) / (j->x - subset->x) > maxDiff)
		{
			maxDiff = fabs(j->z - subset->z) / (j->x - subset->x);
		}
	}
	if (maxDiff < DBL_EPSILON)
	{
		maxDiff = 1.0;
	}
	return maxDiff;
}


double Solver::TrialSubset::calc_relative_difference_between_nodes(const Trial& lhs_node,
	const Trial& rhs_node)
{
	return fabs(rhs_node.z - lhs_node.z) / fabs(rhs_node.x - lhs_node.x);
}


double Solver::TrialSubset::calc_max_difference_between_neighbours(const Trial& another_trial)
{
	auto position = subset.find(another_trial);

	double maxDiff = 0.0;

	if (position != subset.begin())
	{
		double difference = calc_relative_difference_between_nodes(
			*std::prev(position), *position);

		if (difference > maxDiff)
		{
			maxDiff = difference;
		}
	}

	if (std::next(position) != subset.end())
	{
		double difference = calc_relative_difference_between_nodes(
			*position, *std::next(position));

		if (difference > maxDiff)
		{
			maxDiff = difference;
		}
	}

	return maxDiff;
}


double Solver::TrialSubset::calc_max_difference_between_all_trials(const Trial& another_trial)
{
	auto position = subset.find(another_trial);

	double maxDiff = 0.0;

	for (auto trial = subset.begin(); trial != subset.end(); ++trial)
	{
		if (trial != position)
		{
			double difference = calc_relative_difference_between_nodes(
				*trial, *position);
			if (difference > maxDiff)
			{
				maxDiff = difference;
			}
		}
	}

	return maxDiff;
}


void Solver::TrialSubset::calc_subset_lower_lip_const()
{
	if (this->subset.size() < 2)
	{
		this->lip_const = 1.0;
	}
	else
	{
		this->lip_const = this->get_subset_max_difference();
	}
}


double Solver::TrialSubset::get_subset_lip_const_lower_estimation(const Trial& another_trial)
{
	if (subset.size() < 2)
	{
		return 1.0;
	}
	else
	{
		double maxDiff = 0.0;

		if (TrialSubset::USE_NEIGHBOUR_NODES)
		{
			if (subset.size() > TrialSubset::SUBSET_SIZE_LIMIT)
			{
				maxDiff = calc_max_difference_between_neighbours(another_trial);
			}
			else
			{
				maxDiff = calc_max_difference_between_all_trials(another_trial);
			}
		}
		else
		{
			maxDiff = calc_max_difference_between_all_trials(another_trial);
		}

		if (maxDiff < DBL_EPSILON)
		{
			return 1.0;
		}

		return maxDiff;
	}
}


bool Solver::TrialSubset::update_subset_lip_const_lower_estimation(const Trial& another_trial)
{
	double new_lip_const =
		get_subset_lip_const_lower_estimation(another_trial);

	if (new_lip_const > lip_const)
	{
		lip_const = new_lip_const;
		return true;
	}

	return false;
}


bool Solver::TrialSubset::check_and_set_zero()
{
	auto prev_min_estimator = min_estimator;

	min_estimator = 0.0;

	if (prev_min_estimator)
	{
		return true;
	}
	return false;
}


bool Solver::TrialSubset::update_subset_min_estimate(const Trial& another_trial)
{
	double zmin = min_estimator;

	if (!zmin)
	{
		zmin = DBL_MAX;
	}

	if (another_trial.z < zmin)
	{
		min_estimator = another_trial.z;
		return true;
	}
	return false;
}


bool Solver::TrialSubset::update_subset_method_parameter(double in_parameter)
{
	auto prev_parameter = method_parameter;

	method_parameter = in_parameter;

	if (prev_parameter != method_parameter)
	{
		return true;
	}
	return false;
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
	calc_count++;

	if (calc_count > CALC_COUNT_LIMIT)
	{
		parameter = DEFAULT_PARAMETER;
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


Solver::MethodData::MethodData(problem_iterator _problem) : problem(_problem)
{
	starting_stamp = tbb::tick_count::now();

	trials.emplace(this->input.left, 0.0, 0);
	trials.emplace(this->input.right, 0.0, 0);

	Trial mid_node((this->input.left + this->input.right) / 2.0);

	mid_node.perform_trial(this->problem);

	trials.insert(mid_node);

	sln_estimator.set_initial_error(fabs(input.right - input.left));

	init_function_stats();

	init_trial_subsets();

	update_trial_subsets(mid_node);

	init_segment_set();
}


void Solver::MethodData::init_function_stats()
{
	for (std::size_t constraint_no = 0;
		constraint_no <= (*problem)->getConstraintsNumber() + 1; constraint_no++)
	{
		function_stats.emplace(constraint_no, input.method_param);
	}
}


void Solver::MethodData::init_trial_subsets()
{
	std::for_each(trials.begin(), trials.end(), [&](const Trial& trial) {
		subsets[trial.nu].insert(trial);
	});
}


void Solver::MethodData::init_segment_set()
{
	for (auto left_node = trials.begin(); left_node != --trials.end(); left_node++)
	{
		auto next = std::next(left_node);

		Interval curr_segment(problem, *left_node, *next);

		curr_segment.calc_characteristic(subsets);

		segment_set.push_back(curr_segment);
	}
}


void Solver::MethodData::add_trial_to_subset(const Trial& another_trial)
{
	subsets[another_trial.nu].insert(another_trial);
}


bool Solver::MethodData::update_trial_subsets(const Trial& another_trial)
{
	bool subset_parameter_changed = false;

	subset_parameter_changed |= update_method_parameters(another_trial);

	subset_parameter_changed |= update_lip_const_lower_estimation(another_trial);

	subset_parameter_changed |= update_min_estimators(another_trial);

	return subset_parameter_changed;
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

bool Solver::MethodData::update_method_parameters(const Trial& another_trial)
{
	update_stats(another_trial);

	return subsets[another_trial.nu].update_subset_method_parameter(
		function_stats[another_trial.nu].parameter);
}


bool Solver::MethodData::update_lip_const_lower_estimation(const Trial& another_trial)
{
	return subsets[another_trial.nu].
		update_subset_lip_const_lower_estimation(another_trial);
}


bool Solver::MethodData::update_min_estimators(const Trial& another_trial)
{
	auto max_nu = subsets.rbegin()->first;

	if (another_trial.nu == max_nu)
	{
		bool min_estimate_changed = false;

		for (auto subset = subsets.begin(); subset != --subsets.end(); subset++)
		{
			min_estimate_changed |= subset->second.check_and_set_zero();
		}


		min_estimate_changed |=
			subsets[another_trial.nu].update_subset_min_estimate(another_trial);

		return min_estimate_changed;
	}

	return false;
}


void Solver::MethodData::update_solution(const Trial& another_trial)
{
	if (another_trial.admissible)
	{
		double error = (*problem)->getReferenceMinError(another_trial.x);

		if (error < sln_estimator.error)
		{
			this->sln_estimator.error = error;
			this->sln_estimator.xmin = another_trial.x;
			this->sln_estimator.zmin = another_trial.z;

			if (error < input.method_eps)
			{
				calc_elapsed_time();
				this->total_trials_count = MethodData::global_trials_count;
				this->total_iterations_count = MethodData::global_iterations_count;
				this->method_finished = true;
			}
		}
	}
}


void Solver::MethodData::update_interval_charateristic(Interval& interval_to_process)
{
	interval_to_process.calc_characteristic(subsets);
}


Solver::Trial Solver::MethodData::get_new_trial()
{
	auto best_interval_iterator = segment_set.begin();

	return best_interval_iterator->create_new_trial();
}


void Solver::MethodData::update_interval_charateristics()
{
	for (auto& segment : segment_set)
	{
		segment.calc_characteristic(subsets);
	}
}

void Solver::MethodData::perform_iteration(MethodData &problem_data)
{
	problem_data.sort_segment_set();

	Trial another_trial = problem_data.get_new_trial();

	problem_data.add_new_trial(another_trial);

	problem_data.split_best_interval(another_trial);

	if (problem_data.do_update_interval_charateristics)
	{
		problem_data.update_interval_charateristics();
	}
}


void Solver::MethodData::update_trial_subsets()
{
	this->update_method_parameters();

	this->update_lower_lip_const();

	this->calc_min_estimators();
}


void Solver::MethodData::update_lower_lip_const()
{
	/*
	for (auto subset = subsets.begin(); subset != --subsets.end(); subset++)
	{
	subset->second.min_estimator = 0.0;
	}

	auto max_nu = subsets.rbegin();

	max_nu->second.calc_subset_lower_lip_const();
	*/
	
}


void Solver::MethodData::update_method_parameters()
{
	for (auto& subset_index : subsets)
	{
		subset_index.second.method_parameter =
			function_stats.at(subset_index.first).parameter;
	}
}



// probably only update is needed!
void Solver::MethodData::calc_min_estimators()
{
	for (auto& subset_index: subsets)
	{
		subset_index.second.calc_subset_min_estimate();

		// subset_index.second.calc_subset_lower_lip_const();
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


void Solver::MethodData::add_new_trial(const Trial& another_trial)
{
	if (method_finished)
	{
		return;
	}

	trials.insert(another_trial);

	trials_count++;

	update_solution(another_trial);

	add_trial_to_subset(another_trial);

	do_update_interval_charateristics = update_trial_subsets(another_trial);
}



bool Solver::MethodData::is_finished() const
{
	return this->method_finished;
}


void Solver::MethodData::merge_segment_set_into(std::list<Interval>& out_segment_set)
{
	segment_set.sort();

	out_segment_set.merge(segment_set);
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
	return trials_count;
}


double Solver::MethodData::get_error_value() const
{
	return sln_estimator.error;
}


void Solver::MethodData::split_best_interval(const Trial &new_trial)
{
	auto best_interval_iterator = segment_set.begin();

	split_interval(std::make_pair(*best_interval_iterator, new_trial));

	segment_set.erase(best_interval_iterator);
}


void Solver::MethodData::split_interval(const std::pair<Interval, Trial> &new_trial_data)
{
	auto old_interval = new_trial_data.first;

	Interval left_subinterval(problem, old_interval.left_node, new_trial_data.second);
	left_subinterval.calc_characteristic(subsets);

	segment_set.push_back(left_subinterval);

	Interval right_subinterval(problem, new_trial_data.second, old_interval.right_node);
	right_subinterval.calc_characteristic(subsets);

	segment_set.push_back(right_subinterval);
}


void Solver::MethodData::sort_segment_set()
{
	segment_set.sort();
}

void Solver::MethodData::parallel_perform_iteration()
{
	sort_segment_set();

	tbb::concurrent_vector<std::pair<Interval, Trial>> new_trials;

	auto segments_end = segment_set.begin();

	if (segment_set.size() < input.num_threads)
	{
		segments_end = segment_set.end();
	}
	else
	{
		std::advance(segments_end, input.num_threads);
	}

	tbb::parallel_do(segment_set.begin(), segments_end,
		[&new_trials](const Interval &target_segment) {
		Trial new_trial = target_segment.create_new_trial();

		new_trials.push_back(std::make_pair(target_segment, new_trial));
	});

	for (const auto &trial_to_process : new_trials)
	{
		auto target_interval_position = std::find(segment_set.begin(),
			segments_end, trial_to_process.first);

		segment_set.erase(target_interval_position);

		add_new_trial(trial_to_process.second);

		split_interval(trial_to_process);		
	}

	if (do_update_interval_charateristics)
	{
		for (auto &segment : segment_set)
		{
			segment.calc_characteristic(subsets);
		}

		do_update_interval_charateristics = false;
	}
}


void Solver::MethodData::calc_elapsed_time()
{
	elapsed_time_in_seconds = (tbb::tick_count::now() - starting_stamp).seconds();
}



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
	//OptProblem* opt_problem = problem->get();

	unsigned problem_dim = (*problem)->getDimention();
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


double Solver::Interval::get_new_point() const
{
	double x = 0.0;
	if (right_node.nu != left_node.nu)
	{
		x = (right_node.x + left_node.x) / 2.0;
	}
	else
	{
		//OptProblem* opt_problem = problem->get();

		unsigned problem_dim = (*problem)->getDimention();
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
	Trial another_trial(get_new_point());

	another_trial.perform_trial(problem);

	return another_trial;
}


double Solver::Interval::sgn(double arg) const
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


void Solver::SimultaneousMethodDataContainer::sort_segment_set()
{
	all_segments.sort();
}

void Solver::SimultaneousMethodDataContainer::add_new_trial(const std::pair<problem_iterator, Trial>& new_trial_data)
{
	auto target_problem_data = problem_series_container.find(ProblemIterator(new_trial_data.first));

	target_problem_data->second.add_new_trial(new_trial_data.second);
}


void Solver::SimultaneousMethodDataContainer::split_interval(
	const std::pair<Interval, Trial>& new_trial_data)
{
	auto old_interval = new_trial_data.first;

	auto target_problem_data =
		problem_series_container.find(ProblemIterator(old_interval.problem));

	Interval left_subinterval(old_interval.problem,
		old_interval.left_node, new_trial_data.second);

	target_problem_data->second.update_interval_charateristic(left_subinterval);

	all_segments.push_back(left_subinterval);

	Interval right_subinterval(old_interval.problem,
		new_trial_data.second, old_interval.right_node);

	target_problem_data->second.update_interval_charateristic(right_subinterval);

	all_segments.push_back(right_subinterval);
}


void Solver::SimultaneousMethodDataContainer::update_segment_set()
{
	all_segments.erase(std::remove_if(all_segments.begin(), all_segments.end(),
		[this](const Interval &segment) {
		auto target_method_data = problem_series_container.find(segment.problem);

		return target_method_data->second.is_finished();
	}), all_segments.end());


	for (auto& segment : all_segments)
	{
		auto target_method_data = problem_series_container.find(segment.problem);

		if (target_method_data->second.do_update_interval_charateristics)
		{
			target_method_data->second.update_interval_charateristic(segment);
		}	
	}
}


void Solver::SimultaneousMethodDataContainer::parallel_perform_iteration()
{
	sort_segment_set();

	auto segments_end = all_segments.begin();

	if (all_segments.size() < MethodData::input.num_threads)
	{
		segments_end = all_segments.end();
	}
	else
	{
		std::advance(segments_end, MethodData::input.num_threads);
	}

	tbb::concurrent_vector<std::pair<Interval, Trial>> new_trials;

	tbb::parallel_do(all_segments.begin(), segments_end,
		[&new_trials](const Interval &target_segment) {
		Trial new_trial = target_segment.create_new_trial();

		new_trials.push_back(std::make_pair(target_segment, new_trial));
	});

	for (const auto &trial_to_process : new_trials)
	{
		auto target_interval_position = std::find(all_segments.begin(),
			all_segments.end(), trial_to_process.first);

		all_segments.erase(target_interval_position);

		++MethodData::global_trials_count;

		add_new_trial(std::make_pair(trial_to_process.first.problem,
			trial_to_process.second));

		split_interval(trial_to_process);
	}

	++MethodData::global_iterations_count;

	update_segment_set();
}


void Solver::SimultaneousMethodDataContainer::merge_segment_sets()
{
	for (auto& method_data_index : problem_series_container)
	{
		auto& problem_data = method_data_index.second;

		problem_data.merge_segment_set_into(all_segments);
	}
}


void Solver::SimultaneousMethodDataContainer::add_problem(problem_iterator problem)
{
	problem_series_container.emplace(problem, problem);
}




void Solver::MethodDataContainer::update_metrics(MetricsContainer& metrics)
{
	update_errors(metrics.errors_by_trials);

	update_portion(metrics.solved_problems_portion_by_trials);
}


void Solver::MethodDataContainer::update_portion(portion_vector& solved_problems_portion)
{
	unsigned num_problems_finished = 0;

	for (const auto& method_data_index : problem_series_container)
	{
		if (method_data_index.second.is_finished())
		{
			++num_problems_finished;
		}
	}

	double portion = 1.0 * num_problems_finished / problem_series_container.size();
	solved_problems_portion.emplace_back(MethodData::global_trials_count, portion);
}

void Solver::MethodDataContainer::update_errors(errors_vector& errors_by_trials)
{
	double average_error = 0.0;
	double max_error = -(DBL_MAX);

	// unsigned num_trials = 0;
	for (const auto& method_data_index : problem_series_container)
	{
		// num_trials += method_data_index.second.get_num_trials();

		double method_error = method_data_index.second.get_error_value();

		average_error += method_error;
		if (method_error > max_error)
		{
			max_error = method_error;
		}
	}
	average_error /= problem_series_container.size();
	errors_by_trials.push_back(
		std::make_pair(MethodData::global_trials_count,
		std::make_pair(average_error, max_error)));
}


bool Solver::MethodDataContainer::is_all_finished()
{
	bool done = true;
	for (const auto& method_data_index : problem_series_container)
	{
		done &= method_data_index.second.is_finished();
	}

	return done;
}


Solver::MetricsContainer::MetricsContainer(Input input)
{
	double dist = fabs(input.right - input.left);
	errors_by_trials.push_back(std::make_pair(0, std::make_pair(dist, dist)));

	// this->max_error_by_trials.emplace_back(0, 0.5 * distance);
}



void Solver::MethodDataContainer::add_problem(problem_iterator problem)
{
	problem_series_container.emplace(problem, problem);
}


void Solver::MethodDataContainer::dump_solving_results(
	std::list<ProblemSolvingResult>& results)
{
	for (const auto& method_data_index : problem_series_container)
	{
		results.emplace_back(method_data_index.second);
	}
}


void Solver::DynamicMethodDataContainer::enqueue_problems(
	const problem_list& problems_to_solve)
{
	for (auto iterator = problems_to_solve.cbegin();
		iterator != problems_to_solve.cend(); ++iterator)
	{
		problem_queue.push(iterator);
	}
}

void Solver::DynamicMethodDataContainer::init_workers(unsigned int num_threads)
{
	for (unsigned int i = 0; i < num_threads; i++)
	{
		take_problem_from_queue();
	}
}


void Solver::DynamicMethodDataContainer::take_problem_from_queue()
{
	if (!problem_queue.empty())
	{
		auto another_problem_to_solve = problem_queue.front();

		problem_queue.pop();

		auto target_problem_data =
			problem_series_container.find(ProblemIterator(another_problem_to_solve));

		active_solving_problems.push_front(std::ref(target_problem_data->second));
	}
}

void Solver::DynamicMethodDataContainer::perform_iteration()
{
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
		++MethodData::global_trials_count;
		MethodData::perform_iteration(solving_problem_reference);
	});

	++MethodData::global_iterations_count;

	for (auto ref_iterator = active_solving_problems.begin();
		ref_iterator != active_solving_problems.end();)
	{
		auto& solving_problem = (*ref_iterator).get();

		if (solving_problem.is_finished())
		{
			ref_iterator = active_solving_problems.erase(ref_iterator);

			take_problem_from_queue();
		}
		else
		{
			++ref_iterator;
		}
	}



}

void Solver::DynamicMethodDataContainer::complete_iteration()
{
	MethodData::global_iterations_count++;

	/*
	 * note: there is no increment in the loop construct
	 */
	for (auto ref_iterator = active_solving_problems.begin();
		ref_iterator != active_solving_problems.end();)
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
	tbb::tick_count start_time;
	tbb::tick_count	finish_time;

	start_time = tbb::tick_count::now();

	SimultaneousMethodDataContainer problems_method_data(problems);

	do
	{
		problems_method_data.parallel_perform_iteration();

		problems_method_data.update_metrics(out.metrics);
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

	DynamicMethodDataContainer problems_method_data(problems);

	do
	{
		problems_method_data.parallel_perform_iteration();

		problems_method_data.update_metrics(out.metrics);
	} while (!problems_method_data.is_all_finished());

	finish_time = tbb::tick_count::now();
	out.elapsed_time = (finish_time - start_time).seconds();

	problems_method_data.dump_solving_results(out.results);
}


void Solver::run_sequential_search(Output& out)
{
	tbb::tick_count start_time;
	tbb::tick_count finish_time;

	start_time = tbb::tick_count::now();

	for (auto problem = this->problems.begin();
		problem != this->problems.end(); ++problem)
	{
		MethodData method_data(problem);
		do
		{
			method_data.parallel_perform_iteration();
		} while (!method_data.is_finished());
		out.add_problem_solving_result(method_data);
	}

	finish_time = tbb::tick_count::now();
	out.elapsed_time = (finish_time - start_time).seconds();
}

void Solver::run_solver(Output& out) {
	Solver::ProblemIterator::problem_list_begin = problems.begin();

	Solver::MethodData::set_method_input(input);

	switch (input.solving_method)
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

