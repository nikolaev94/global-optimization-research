
#ifndef SOLVER_H
#define SOLVER_H

#include "opt_problem.h"

#include <chrono>
#include <queue>
#include <map>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>


#include <tbb/concurrent_vector.h>
#include <tbb/parallel_do.h>
#include <tbb/parallel_sort.h>
#include <tbb/mutex.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>


//#include <mpi.h>



class Solver
{
private:
	struct Trial;
	struct MethodData;

public:
	typedef std::list< std::shared_ptr<OptProblem> > problem_list;
	typedef problem_list::const_iterator problem_iterator;

	struct ProblemSolvingResult;
	struct MetricsContainer;

	enum SolvingMethod
	{
		SEQUENTIAL, DYNAMIC, SIMULTANEOUS
	};

	struct Input
	{
		double left, right;
		double method_eps, method_param;
		unsigned int num_threads;
		unsigned int problems_dimension;
		SolvingMethod solving_method;

		Input() :
			left(0.0), right(1.0), method_eps(1.e-3), method_param(3.0),
			num_threads(1), problems_dimension(2),
			solving_method(SolvingMethod::SEQUENTIAL) {}

		Input(double in_left, double in_right, double in_method_eps,
			double in_method_param, unsigned int in_num_threads,
			unsigned int in_problems_dimension = 2,
			SolvingMethod in_solving_method = SolvingMethod::SEQUENTIAL) :
			left(in_left), right(in_right), method_eps(in_method_eps),
			method_param(in_method_param), num_threads(in_num_threads),
			problems_dimension(in_problems_dimension), solving_method(in_solving_method) {}
	};

	struct ProblemSolvingResult
	{
		problem_iterator solved_problem;

		double xmin;
		double zmin;
		double error;
		bool correctness;
		unsigned trials;
		double elapsed_time;

		/*
		 * Total number of algorithm iterations completed before the problem was resolved
		 */
		unsigned iterations;
		/*
		 * Total number of all problems trials performed until the problem was resolved
		 */
		unsigned total_trials;
		
		ProblemSolvingResult() : correctness(true) {}
		ProblemSolvingResult(const MethodData& method_data)
			: solved_problem(method_data.problem)
		{
			// this->iterations = method_data.iterations_count;

			this->xmin = method_data.sln_estimator.xmin;
			this->zmin = method_data.sln_estimator.zmin;
			this->error = method_data.sln_estimator.error;
			this->trials = method_data.trials_count;
			this->iterations = method_data.total_iterations_count;
			this->total_trials = method_data.total_trials_count;
			this->elapsed_time = method_data.elapsed_time_in_seconds;
			this->correctness = true;
		}
	};

	

	typedef std::vector <std::pair<unsigned, std::pair<double, double>>> errors_vector;
	typedef std::vector <std::pair<unsigned, double>> portion_vector;

	struct MetricsContainer
	{
		errors_vector errors_by_trials;
		portion_vector solved_problems_portion_by_trials;

		MetricsContainer() {}
		MetricsContainer(Input);
	};

	struct Output
	{
		double elapsed_time;
		const Input& original_input;
		std::list<ProblemSolvingResult> results;
		MetricsContainer metrics;

		Output(const Input& input) : original_input(input),
			metrics(input), elapsed_time(0.0) {}

		void add_problem_solving_result(const MethodData&);
		void dump_results_to_file(const std::string&);
		void dump_error_metrics_by_trials_to_file(const std::string&);
		void dump_solved_problem_portion_by_trials_to_file(const std::string&);
	};

private:
	


	struct Trial
	{
		double x, z;
		unsigned nu;
		bool admissible;
		
		Trial() {}
		Trial(double _x) : x(_x) {}
		Trial(const Trial& src) : x(src.x), z(src.z), nu(src.nu),
			admissible(src.admissible) {}
		
		Trial(double _x, double _z, unsigned _nu) : x(_x), z(_z), nu(_nu) {}

		void perform_trial(problem_iterator problem);
		bool operator< (const Trial& comp) const { return this->x < comp.x; }
	};

	struct TrialSubset
	{
		/*
		 * Use only neighbour trials to calculate lower lip const
		 */
		static bool USE_NEIGHBOUR_NODES;
		
		std::set<Trial> subset;
		double lip_const;
		double min_estimator;
		const std::size_t SUBSET_SIZE_LIMIT = 10;

		double method_parameter;

		TrialSubset(): lip_const(1.0), min_estimator(0.0) {}

		/*TrialSubset(const std::set<Trial>& in_subset) : subset (in_subset),
			lip_const(1.0), min_estimator(0.0) {}*/

		void insert(const Trial& trial) { this->subset.insert(trial); }

		bool check_and_set_zero();
		bool update_subset_min_estimate(const Trial&);

		bool update_subset_lip_const_lower_estimation(const Trial&);

		bool update_subset_method_parameter(double in_parameter);

		void calc_subset_min_estimate();
	private:
		unsigned int dimension;

		double get_subset_lip_const_lower_estimation(const Trial&);

		double calc_max_difference_between_neighbours(const Trial&);
		double calc_max_difference_between_all_trials(const Trial&);
		double calc_relative_difference_between_nodes(const Trial&, const Trial&);
	};

	typedef std::map< unsigned, TrialSubset > trial_subsets;
	
	struct Interval
	{
		problem_iterator problem;
		Trial left_node, right_node;
		double charact;

		double greater_nu_lipconst;
		double greater_nu_method_param;

		void calc_characteristic(const trial_subsets&);

		Trial create_new_trial() const;

		Interval() = delete;

		Interval(problem_iterator in_problem) : problem(in_problem) {}
		Interval(problem_iterator in_problem, const Trial& in_left_node,
			const Trial& in_right_node)
			: problem(in_problem), left_node(in_left_node), right_node(in_right_node) {}

		Interval(const Interval& src) : problem(src.problem), left_node(src.left_node), right_node(src.right_node),
			charact(src.charact), greater_nu_lipconst(src.greater_nu_lipconst), greater_nu_method_param(src.greater_nu_method_param) {}

		// By default, compare intervals by charateristics
		bool operator< (const Interval& interval) const
		{
			return this->charact > interval.charact;
		}

		bool operator== (const Interval& interval) const
		{
			return this->right_node.x == interval.right_node.x;
		}

	private:
		double sgn(double arg) const;
		double get_new_point() const;
		double get_interval_length();
	};

	struct IntervalNodeComparator
	{
		bool operator() (const Interval& lhs, const Interval& rhs) const
		{
			return lhs.right_node < rhs.right_node;
		}
	};

	struct IntervalCharateristicCompartor
	{
		bool operator() (const Interval& lhs, const Interval& rhs) const
		{
			return lhs.charact > rhs.charact;
		}
	};

	struct FunctionStatsInfo
	{
		static double initial_parameter;

		const double DEFAULT_PARAMETER = 2.0;
		const unsigned CALC_COUNT_LIMIT = 20;

		unsigned calc_count;
		double parameter;

		void update();
		FunctionStatsInfo() : calc_count(0), parameter(initial_parameter) {}
		FunctionStatsInfo(double in_parameter) : calc_count(0), parameter(in_parameter) {}
	};

	typedef std::map<unsigned, FunctionStatsInfo> function_stats_t;

	struct SolutionEstimator
	{
		double xmin;
		double zmin;
		double error;

		void set_initial_error(double in_error)
		{
			error = in_error;
		}

		SolutionEstimator() : xmin(0.0), zmin(DBL_MAX), error(DBL_MAX) {}
	};

	struct MethodData
	{
		static Input input;
		static unsigned global_iterations_count;

		static unsigned global_trials_count;

		static tbb::mutex mu;
		//static tbb::atomic<unsigned int> global_trials_count;

		static void set_method_input(const Input&);

		static void perform_iteration(MethodData&);


		//Trial get_best_interval();

		//void perform_iteration();

		//bool operator==(const MethodData&) const;

		unsigned trials_count = 0;
		unsigned total_trials_count = 0;
		unsigned total_iterations_count = 0;

		double elapsed_time_in_seconds = 0.0;

		SolutionEstimator sln_estimator;
		problem_iterator problem;

		bool do_update_interval_charateristics = false;

		void update_interval_charateristic(Interval&);

		void update_interval_charateristics();

		void merge_segment_set_into(std::list<Interval>&);

		void parallel_perform_iteration();
		
		void dump_solving_result(std::list<ProblemSolvingResult>& results);

		void add_new_trial(const Trial&);
		
		bool is_finished() const;

		unsigned get_num_trials() const;

		double get_error_value() const;

		void sort_segment_set();

		MethodData(problem_iterator);

	private:
		tbb::tick_count starting_stamp;
		tbb::tick_count finishing_stamp;

		bool method_finished = false;

		std::set<Trial> trials;

		std::map<unsigned, TrialSubset> subsets;

		std::list<Interval> segment_set;

		std::map<unsigned, FunctionStatsInfo> function_stats;
		
		void init_function_stats();
		void init_trial_subsets();
		void init_segment_set();

		void add_trial_to_subset(const Trial&);

		void update_solution(const Trial&);

		bool update_trial_subsets(const Trial&);
		void update_stats(const Trial&);
		bool update_method_parameters(const Trial&);
		bool update_min_estimators(const Trial&);
		bool update_lip_const_lower_estimation(const Trial&);

		Trial get_new_trial();

		void split_best_interval(const Trial&);

		void split_interval(const std::pair<Interval, Trial>&);

		void calc_elapsed_time();
	};

	struct MethodDataKey
	{
		problem_iterator problem;

		bool operator==(const MethodDataKey& cmp) const;

		MethodDataKey() = delete;

		MethodDataKey(problem_iterator in_problem) : problem(in_problem) {}
	};


	struct MethodDataKeyHasher
	{
		std::size_t operator() (const MethodDataKey&) const;
	};


	struct MethodDataContainer
	{
		bool is_all_finished();

		void dump_solving_results(std::list<ProblemSolvingResult>&);

		void update_metrics(MetricsContainer&);

		MethodDataContainer(const problem_list& in_problems)
		{
			for (auto problem = in_problems.begin();
				problem != in_problems.end(); ++problem)
			{
				add_problem(problem);
			}
		}

	private:
		void add_problem(problem_iterator);

	protected:
		std::unordered_map<MethodDataKey, MethodData, MethodDataKeyHasher> problem_series_container;

		std::vector<problem_iterator> problem_order_in_series;

		void update_errors(errors_vector&);
		void update_portion(portion_vector&);
	};


	struct SimultaneousMethodDataContainer : public MethodDataContainer
	{	
		void add_problem(problem_iterator);

		void parallel_perform_iteration();

		SimultaneousMethodDataContainer(const problem_list &in_problems)
			: MethodDataContainer(in_problems)
		{
			merge_segment_sets();
		}

	private:
		std::list<Interval> all_segments;

		//std::vector<Interval> all_segments;

		void merge_segment_sets();

		void sort_segment_set();

		void add_new_trial(const std::pair<problem_iterator, Trial>&);

		void split_interval(const std::pair<Interval, Trial>&);

		void update_segment_set();
	};

	struct DynamicMethodDataContainer : public MethodDataContainer
	{
		void parallel_perform_iteration();

		void complete_iteration();

		DynamicMethodDataContainer(const problem_list& in_problems)
			: MethodDataContainer(in_problems)
		{
			enqueue_problems(in_problems);

			init_workers(MethodData::input.num_threads);
		}
	private:

		void init_workers(unsigned);

		void enqueue_problems(const problem_list&);

		void take_problem_from_queue();


		//std::list<problem_iterator> active_solving_problems1;

		std::list<std::reference_wrapper<MethodData>> active_solving_problems;

		std::queue<problem_iterator> problem_queue;
	};


	Input input;
	problem_list problems;

	tbb::task_scheduler_init scheduler;
	void init_tbb(int number_threads = 1);

	typedef std::set<Interval>::const_iterator target_interval;

	bool check_all_problem(const MethodData&);

	void run_simultaneous_search(Output&);
	void run_dynamic_search(Output&);
	void run_sequential_search(Output&);

public:

	Solver() = delete;
	Solver(const Input&, const problem_list&);
	~Solver();

	void run_solver(Output&);
};

#endif
