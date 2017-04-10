
/*
#include "all.h"
#include "infrastructure.h"

void generateShekelFunctionsToFile(const char* filename) {
	std::uniform_real_distribution<double> un(1.0, 3.0); // ki
	std::uniform_real_distribution<double> un2(0.0, 10.0); // ai, ci
	std::random_device rd;
	std::default_random_engine re(rd());
	auto rnd = std::bind(un, re);
	auto rnd2 = std::bind(un2, re);
	FILE* file = fopen(filename, "w");
	fprintf(file, "%d\n", 30);
	for (int i = 0; i < MAX_FUNCSNUM; i++) {
		for (int j = 0; j < 10; j++) {
			fprintf(file, "%f ", rnd());
		}
		for (int j = 0; j < 20; j++) {
			fprintf(file, "%f ", rnd2());
		}
		fputc('\n', file);
	}
	fclose(file);
}

void generateHillFunctionsToFile(const char* filename) {
	std::uniform_real_distribution<double> un(-1.0, 1.0);
	std::random_device rd;
	std::default_random_engine re(rd());
	auto rnd = std::bind(un, re);
	FILE* file = fopen(filename, "w");
	fprintf(file, "%d\n", 28);
	for (int i = 0; i < MAX_FUNCSNUM; i++) {
		for (int j = 0; j < 28; j++) {
			fprintf(file, "%f ", rnd());
		}
		fputc('\n', file);
	}
	fclose(file);
}

void generateIntegralFunctionsToFile(const char* filename) {
	std::uniform_real_distribution<double> un(-1.0, 1.0); // ai, bi
	std::uniform_real_distribution<double> un2; // alpha
	std::random_device rd;
	std::default_random_engine re(rd());
	auto rnd = std::bind(un, re);
	auto rnd2 = std::bind(un2, re);
	FILE* file = fopen(filename, "w");
	fprintf(file, "%d\n", 29);
	for (int i = 0; i < MAX_FUNCSNUM; i++) {
		for (int j = 0; j < 28; j++) {
			fprintf(file, "%f ", rnd());
		}
		fprintf(file, "%f\n", rnd2());
	}
	fclose(file);
}

void readData(const char* filename,
	int problemsNum, int funcsNum, double delta,
	std::vector<std::unique_ptr<Function>>* series) {
	FILE* file = fopen(filename, "r");
	int lineSize = 0;
	fscanf(file, "%d", &lineSize);
	double* buf = new double[problemsNum * funcsNum * lineSize];
	for (int i = 0; i < problemsNum * funcsNum * lineSize; i++) {
		fscanf(file, "%lf", &buf[i]);
	}
	fclose(file);
	series->clear();
	switch (lineSize) {
	case 28:
		for (int i = 0; i < problemsNum; i++) {
			std::unique_ptr<Function> temp(new Hill(funcsNum));
			series->push_back(std::move(temp));
		}
		break;
	case 29:
		for (int i = 0; i < problemsNum; i++) {
			std::unique_ptr<Function> temp(new Integral(funcsNum));
			series->push_back(std::move(temp));
		}
		break;
	case 30:
		for (int i = 0; i < problemsNum; i++) {
			std::unique_ptr<Function> temp(new CShekel(funcsNum));
			series->push_back(std::move(temp));
		}
		break;
	default:
		break;
	}
	int j = 0;
	for (size_t i = 0; i < series->size(); i++) {
		double* slice = new double[lineSize * funcsNum];
		std::copy(&buf[j], &buf[j + lineSize * funcsNum],
			stdext::checked_array_iterator<double*>(slice, lineSize * funcsNum));
		series->at(i)->setData(slice, delta);
		delete[] slice;
		j += lineSize * funcsNum;
	}
	delete[] buf;
}

void solveProblem(Function* problem, InputParameters inputArgs, OutputResult& result,
	unsigned* stats, ZedGraph::PointPairList^ pointList, double times[]) {
	OutputResult enumResult;
	double startTime = omp_get_wtime();
	enumerativeSearch(problem, inputArgs, enumResult);
	times[0] = omp_get_wtime() - startTime;
	startTime = omp_get_wtime();
	globalSearch(problem, inputArgs, enumResult.xmin, result, stats, pointList);
	times[1] = omp_get_wtime() - startTime;
}

void updateSolvedCounts(std::vector<unsigned>& counts, unsigned count) {
	bool thresholdIsReached = false;
	for (size_t i = 0; i < counts.size(); i++) {
		if (thresholdIsReached) {
			counts[i]++;
			continue;
		}
		unsigned threshold = i * 10;
		if (count <= threshold) {
			counts[i]++;
			thresholdIsReached = true;
		}
	}
}

void solveSeries(std::vector<std::unique_ptr<Function>>* series,
	InputParameters input, const char* outputFilename,
	std::vector<unsigned>& solvedCounts,
	std::vector<double>& sumErrors,
	std::vector<double>& maxErrors) {
	solvedCounts.assign(SOLVED_COUNTS_SIZE, 0);
	if (input.solvingMethod == SeriesSolvingMethod::SEQUENTIAL) {
		OutputResult* enumResult = new OutputResult[series->size()];
		OutputResult* methodResult = new OutputResult[series->size()];
		double* enumElapsedTimePerTask = new double[series->size()];
		double* methodElapsedTimePerTask = new double[series->size()];

		for (size_t taskIndex = 0; taskIndex < series->size(); taskIndex++) {
			Function* task = series->at(taskIndex).get();
			double startTime = omp_get_wtime();
			enumerativeSearch(task, input, enumResult[taskIndex]);
			enumElapsedTimePerTask[taskIndex] = omp_get_wtime() - startTime;
			startTime = omp_get_wtime();
			globalSearch(task, input, enumResult[taskIndex].xmin, methodResult[taskIndex]);
			methodElapsedTimePerTask[taskIndex] = omp_get_wtime() - startTime;
			updateSolvedCounts(solvedCounts, methodResult[taskIndex].itersCount);
		}

		FILE* file = fopen(outputFilename, "w");
		fprintf(file, "№ задачи;Число итераций;x*;z*;Время работы;Число итераций;x*;z*;Время работы\n");
		for (size_t taskIndex = 0; taskIndex < series->size(); taskIndex++) {
			fprintf(file, "%d;%d;%f;%f;%f;", taskIndex + 1, enumResult[taskIndex].itersCount,
				enumResult[taskIndex].xmin, enumResult[taskIndex].zmin, enumElapsedTimePerTask[taskIndex]);
			fprintf(file, "%d;%f;%f;%f\n", methodResult[taskIndex].itersCount, methodResult[taskIndex].xmin,
				methodResult[taskIndex].zmin, methodElapsedTimePerTask[taskIndex]);
		}
		fclose(file);

		delete[] enumResult;
		delete[] methodResult;
		delete[] enumElapsedTimePerTask;
		delete[] methodElapsedTimePerTask;

	} else {
		unsigned blocksNum = series->size() / input.blockSize;
		unsigned totalTrialsNum = 0;
		unsigned totalItersLowerBound = 0;

		for (size_t blockIndex = 0; blockIndex < blocksNum; blockIndex++) {
			unsigned** stats = new unsigned*[input.blockSize];
			Function** block = new Function*[input.blockSize];
			OutputResult* enumResult = new OutputResult[input.blockSize];
			double* slnEstimators = new double[input.blockSize];
			OutputResult* methodResult = new OutputResult[input.blockSize];

			for (size_t problemIndex = 0; problemIndex < input.blockSize; problemIndex++) {
				block[problemIndex] = series->at(blockIndex * input.blockSize + problemIndex).get();
				enumerativeSearch(block[problemIndex], input, enumResult[problemIndex]);
				slnEstimators[problemIndex] = enumResult[problemIndex].xmin;
			}

			for (size_t problemIndex = 0; problemIndex < input.blockSize; problemIndex++) {
				stats[problemIndex] = new unsigned[block[problemIndex]->getFuncsNum()];
			}

			switch (input.solvingMethod) {
			case SeriesSolvingMethod::SIMULTANEOUS:
				globalSearchSimultaneous(block, input, slnEstimators, methodResult, stats, totalTrialsNum, totalItersLowerBound,
					sumErrors, maxErrors);
				break;
			case SeriesSolvingMethod::DYNAMIC:
				globalSearchDynamic(block, input, slnEstimators, methodResult, stats, totalTrialsNum, totalItersLowerBound,
					sumErrors, maxErrors);
				break;
			default:
				break;
			}

			FILE* file = fopen(outputFilename, "w");
			fprintf(file, "№ задачи;x*;z*;Число итераций;x*;z*;Число итераций;Общее число итераций;Число подсчетов функций\n");
			for (size_t problemIndex = 0; problemIndex < input.blockSize; problemIndex++) {
				fprintf(file, "%u;%f;%f;%d;", blockIndex * input.blockSize + problemIndex + 1,
					enumResult[problemIndex].xmin, enumResult[problemIndex].zmin, enumResult[problemIndex].itersCount);
				fprintf(file, "%f;%f;%d;%u;", methodResult[problemIndex].xmin, methodResult[problemIndex].zmin,
					methodResult[problemIndex].itersCount, methodResult[problemIndex].totalItersCount);
				fprintf(file, "%u", stats[problemIndex][0]);
				for (int funcIndex = 1; funcIndex < sizeof(stats[problemIndex]) / sizeof(stats[problemIndex][0]) + 1; funcIndex++) {
					fprintf(file, ",%u", stats[problemIndex][funcIndex]);
				}
				fprintf(file, "\n");
				updateSolvedCounts(solvedCounts, methodResult[problemIndex].totalItersCount);
			}
			fclose(file);

			for (size_t problemIndex = 0; problemIndex < input.blockSize; problemIndex++) {
				delete[] stats[problemIndex];
			}
			delete[] stats;
			delete[] slnEstimators;
			delete[] methodResult;
			delete[] enumResult;
			delete[] block;
		}
	}
}
*/
