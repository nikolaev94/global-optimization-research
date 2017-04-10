
#ifndef SOLVER_H
#define SOLVER_H

#include "all.h"
#include "opt_problem.h"

#include <queue>
#include <map>
#include <set>
#include <unordered_map>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_do.h>
#include <tbb/mutex.h>
#include <tbb/spin_mutex.h>
#include <tbb/tick_count.h>


/*
#include <cilk\cilk.h>
#include <cilk\cilk_stub.h>
#include <cilk\reducer.h>
#include <cilk\reducer_opadd.h>
*/

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
		SolvingMethod solving_method;

		Input() :
			left(0.0), right(1.0), method_eps(1.e-3), method_param(3.0),
			num_threads(1), solving_method(SolvingMethod::SEQUENTIAL) {}

		Input(double _left, double _right, double _method_eps,
			double _method_param, unsigned int _num_threads,
			SolvingMethod _solving_method = SolvingMethod::SEQUENTIAL) :
			left(_left), right(_right), method_eps(_method_eps), method_param(_method_param),
			num_threads(_num_threads), solving_method(_solving_method) {}
	};

	struct ProblemSolvingResult
	{
		problem_iterator solved_problem;

		double xmin, zmin;
		double error;
		bool correctness;
		unsigned trials;
		/*
		 * Total number of algorithm iterations completed before the problem was resolved
		 */
		unsigned iterations;
		/*
		 * Total number of all problems trials performed until the problem was resolved
		 */
		unsigned total_trials;
		
		ProblemSolvingResult() : correctness(true) {}
		ProblemSolvingResult(const MethodData& method_data) : solved_problem(method_data.problem)
		{
			// this->iterations = method_data.iterations_count;
			this->xmin = method_data.sln_estimator.xmin;
			this->zmin = method_data.sln_estimator.zmin;
			this->error = method_data.sln_estimator.error;
			this->trials = method_data.trials_count;
			this->iterations = method_data.total_iterations_count;
			this->total_trials = method_data.total_trials_count;
			this->correctness = true;
		}
	};

	typedef std::vector <std::pair<unsigned, std::pair<double, double>>> errors_vector;
	typedef std::vector <std::pair<unsigned, double>> portion_vector;

	struct MetricsContainer
	{
		/*
		std::vector<std::pair<unsigned, double>> average_error_by_trials;
		std::vector<std::pair<unsigned, double>> max_error_by_trials;
		*/

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

		Output(const Input& input) : original_input(input), metrics(input), elapsed_time(0.0) {}
		void dump_results_to_file(const std::string&);
		void dump_error_metrics_by_trials_to_file(const std::string&);
		void dump_solved_problem_portion_by_trials_to_file(const std::string&);
	};

private:
	struct ProblemIterator
	{
		static problem_iterator problem_list_begin;
		problem_iterator problem_ptr;

		bool operator<(const ProblemIterator& cmp) const;

		ProblemIterator(problem_iterator iterator) : problem_ptr(iterator) {}
	};

	

	struct Trial
	{
		double x, z;
		unsigned nu;
		bool admissible;
		
		Trial() {}
		Trial(double _x) : x(_x) {}
		Trial(const Trial& src) : x(src.x), z(src.z), nu(src.nu), admissible(src.admissible) {}
		
		Trial(double _x, double _z, unsigned _nu) : x(_x), z(_z), nu(_nu) {}

		void perform_trial(problem_iterator problem);
		bool operator< (const Trial& comp) const { return this->x < comp.x; }
	};

	struct TrialSubset
	{
		/*
		 * Use only neighbour trials to calculate lower lip const
		 * This optimization can be used only for one-function problems
		 */
		static bool use_neighbour_nodes;
		std::set<Trial> subset;
		double lip_const;
		double min_estimator;
		const size_t SUBSET_SIZE_LIMIT = 20;

		double method_parameter;

		TrialSubset() : lip_const(1.0), min_estimator(0.0) {}
		TrialSubset(const std::set<Trial>& _subset) : subset (_subset), lip_const(1.0), min_estimator(0.0) {}
		void insert(const Trial& trial) { this->subset.insert(trial); }

		void calc_subset_lower_lip_const();
		void calc_subset_min_estimate();
	private:
		double get_subset_max_difference();

		double calc_max_difference_neighbours();
		double calc_max_difference_all_trials();
	};


	typedef std::map< unsigned, TrialSubset > trial_subsets;
	

	struct Interval
	{
		problem_iterator problem;
		Trial left_node, right_node;
		double charact;

		double greater_nu_lipconst;
		double greater_nu_method_param;

		Interval() {}
		Interval(problem_iterator _problem) : problem(_problem) {}
		Interval(problem_iterator _problem, const Trial& _left_node, const Trial& _right_node)
			: problem(_problem), left_node(_left_node), right_node(_right_node) {}

		// By default, compare intervals by chars
		bool operator< (const Interval& interval) const { return this->charact > interval.charact; }

		void calc_characteristic(const trial_subsets&);
		double get_new_point() const;
		Trial create_new_trial() const;
	private:
		double get_interval_length();
	};

	typedef std::multiset<Interval>::iterator interval_iterator;

	struct IntervalNodeComparator {
		bool operator() (const Interval& lhs, const Interval& rhs) const
		{
			return lhs.right_node < rhs.right_node;
		}
	};

	/*

	struct IntervalContainer
	{
	void split_interval(const Interval&, const Trial&);

	private:
	std::set<Interval> segments;
	//std::set<Interval, IntervalNodeComparator> segments;
	// std::multiset<Trial> segment_set;
	};

	*/

	struct FunctionStatsInfo
	{
		static double initial_parameter;

		const double DEFAULT_PARAMETER = 2.0;
		const unsigned CALC_COUNT_LIMIT = 20;

		unsigned calc_count;
		double parameter;

		void update();
		FunctionStatsInfo() : calc_count(0), parameter(initial_parameter) {}
		FunctionStatsInfo(double _parameter) : calc_count(0), parameter(_parameter) {}
	};

	typedef std::map<unsigned, FunctionStatsInfo> function_stats_t;

	struct SolutionEstimator
	{
		double xmin;
		double zmin;
		double error;

		SolutionEstimator() : xmin(0.0), zmin(DBL_MAX), error(DBL_MAX) {}
	};


	struct MethodData
	{
		static Input input;
		static void set_method_input(const Input&);
		static unsigned global_iterations_count;
		static tbb::atomic<unsigned> global_trials_count;

		// unsigned iterations_count = 0;
		unsigned trials_count = 0;
		unsigned total_trials_count = 0;
		unsigned total_iterations_count = 0;
		SolutionEstimator sln_estimator;
		problem_iterator problem;

		void add_trial(const Trial&);
		void update_trial_subsets();
		
		void construct_segment_set(std::multiset<Interval>&) const;
		bool is_finished() const;

		unsigned get_num_trials() const;
		double get_error_value() const;

		MethodData(problem_iterator);

		// bool operator==(const MethodData&) const;

	private:
		bool method_finished = false;
		std::set<Trial> trials;
		trial_subsets subsets;


		//std::multiset<>



		// std::multiset<Interval>;

		function_stats_t function_stats;
		
		void init_function_stats();
		void init_trial_subsets();

		void add_trial_to_subset(const Trial&);
		void update_stats(const Trial&);
		void update_solution(const Trial&);

		void update_method_parameters();
		void calc_lower_lip_const();
		void calc_min_estimators();
	};

	struct ProblemIteratorHasher
	{
		size_t operator() (const problem_iterator&) const;
	};

	struct MethodDataContainer
	{
		bool is_all_finished();
		void update_metrics(MetricsContainer&);
		void dump_solving_results(std::list<ProblemSolvingResult>&);
	protected:
		std::map<ProblemIterator, MethodData> problem_series_container;

		void update_errors(errors_vector&);
		void update_portion(portion_vector&);

		// std::unordered_map<problem_iterator, MethodData, ProblemIteratorHasher> problem_series_container;
	};


	struct SimultaneousMethodDataContainer : public MethodDataContainer
	{	
		static tbb::mutex mutex;
		void construct_and_merge_sets(std::multiset<Interval>&);
		void add_new_trial(problem_iterator, Trial);
		void complete_iteration();
		void add_problem(problem_iterator);
		void parallel_perform_iteration(const std::multiset<Interval>&);

		SimultaneousMethodDataContainer() {}
	};

	struct DynamicMethodDataContainer : public MethodDataContainer
	{
		// bool queue_is_empty();

		void enqueue_problem(problem_iterator);
		void init_workers(unsigned);
		void perform_iteration();
		void parallel_perform_iteration();
		void complete_iteration();

		DynamicMethodDataContainer() {}
	private:
		void take_problem_from_queue();

		std::list<std::reference_wrapper<MethodData>> active_solving_problems;
		std::queue<problem_iterator> problem_queue;
		// std::list<MethodData> active_solving_problems;
	};


	/*class DynamicMethodIterationExecutor
	{
	private:
		std::list<std::reference_wrapper<MethodData>> method_data;
	public:
		DynamicMethodIterationExecutor(std::list<std::reference_wrapper<MethodData>>&);
		void operator()(const tbb::blocked_range<unsigned>& range);
	};*/

	

	Input input;
	problem_list problems;

	tbb::task_scheduler_init scheduler;
	void init_tbb(int number_threads = 1);

	typedef std::set<Interval>::const_iterator target_interval;

	static double sgn(double arg);

	bool check_all_problem(const MethodData&);

	void run_simultaneous_search(Output&);
	void run_dynamic_search(Output&);
	void run_sequential_search(Output&);

public:

	Solver();
	Solver(const Input&);
	Solver(const Input&, const problem_list&);
	~Solver();

	void run_solver(Output&);
};

#endif
