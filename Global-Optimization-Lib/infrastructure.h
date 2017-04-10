
/*
#ifndef INFRASTRUCTURE_H
#define INFRASTRUCTURE_H

#include "function.h"
#include "shekel.h"
#include "hill.h"
#include "integral.h"
#include "method.h"

enum GraphAbsciss {
	ITERATIONS,
	TRIALS
};

enum GraphType {
	SOLVED_TASKS_BY_ITERATIONS,
	MEAN_ERROR_AMOUNT_BY_ITERATIONS,
	MAX_ERROR_AMOUNT_BY_ITERATIONS
};

const unsigned SOLVED_COUNTS_SIZE = 1000;

void generateShekelFunctionsToFile(const char* filename);

void generateHillFunctionsToFile(const char* filename);

void generateIntegralFunctionsToFile(const char* filename);

void readData(const char* filename, int problemsNum, int funcsNum, double delta, std::vector<std::unique_ptr<Function>>* series);

void solveSeries(std::vector<std::unique_ptr<Function>>* series,
	InputParameters input, const char* outputFilename,
	std::vector<unsigned>& solvedCounts,
	std::vector<double>& sumErrors,
	std::vector<double>& maxErrors);

void solveProblem(Function* task, InputParameters inputArgs, OutputResult& result,
	unsigned* stats, ZedGraph::PointPairList^ pointList, double times[]);
	

#endif
*/
