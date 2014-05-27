/* movingPeaks.cpp 
 * Copyright (C) 2014 Pavel Grunt 
 * This is free software; you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License. 
 * 
 * LEM in Moving Peaks Benchmark enviroment
 * 
 */  


#include "Lem.h"
#include "movpeaks.h"
#include <iostream>
#include <sstream>
#include <utility>
#include <cstring>


#ifdef SVM
#ifndef CLSVAL
#define CLSVAL 1
#endif
#endif
#ifdef COMBINED
#ifndef CLSVAL
#define CLSVAL 2
#endif
#endif

#ifndef CLSVAL
#define CLSVAL 0
#endif

namespace {
	int dimensions = 5;
	int peaks = 10, change_freq = 5000;
	std::pair<double, double> userrange({0, 100});
	Specimen specimen;
	unsigned groups = 2;
	double maxGroupSize = 0.05, byLearning = 0.7, byEvolution = 0.2;
	unsigned populationSize = 100, generations = 5000, progressProbe = 4, mutationProbe = 800;
	unsigned minimize = 0;
	unsigned classifierType = CLSVAL;

	
	Specimen createSpecimen(std::pair<double, double> & range, unsigned dimensions) {
		Specimen s(dimensions);
		for (unsigned i = 0; i < dimensions; ++i)
		{
			s[i] = range;
		}
		return s;
	}

	double fitness(std::vector<double> &v) {
		return eval_movpeaks(&v[0]);
	}
	double dummyFitness(std::vector<double> &v) {
		return dummy_eval(&v[0]);
	}

	void helpMsg() {
		generations = 0;
		std::cout << "LEM - Moving Peaks Benchmark \n" ;
		std::cout << "Pocet vrcholu    \t-peaks <unsigned>\n";
		std::cout << "Frekvence zmeny  \t-freq <unsigned>\n";
		std::cout << "Velikost populace\t-p <unsigned>\n";
		std::cout << "Pocet generaci   \t-g <unsigned>\n";
		std::cout << "Pocet dimenzi    \t-d <unsigned>\n";
		std::cout << "Spodni hranice   \t-lb <double>\n";
		std::cout << "Horni hranice    \t-ub <double>\n";
		std::cout << "Pocet skupin     \t-grs <unsigned>\n";
		std::cout << "Velikost skupiny\t-mg (0,0; 1,0)\n";
		std::cout << "=====================================\n";
	}

// COMMAND LINE
	void parseCommandLine(int argc, char const *argv[]) {
		std::string input;
		if (argc == 2) helpMsg();
		for (int argsi = 0; argsi < (argc - 1); ++argsi) {
			if (0 == strcmp(argv[argsi],"-p")) {
				input = argv[++argsi];
				std::stringstream(input) >> populationSize;
			} else if (0 == strcmp(argv[argsi],"-h")) {
				helpMsg();
			} else if (0 == strcmp(argv[argsi],"-lb")) {
				input = argv[++argsi];
				std::stringstream(input) >> userrange.first;
			} else if (0 == strcmp(argv[argsi],"-ub")) {
				input = argv[++argsi];
				std::stringstream(input) >> userrange.second;
			} else if (0 == strcmp(argv[argsi],"-pp")) {
				input = argv[++argsi];
				std::stringstream(input) >> progressProbe;
			} else if (0 == strcmp(argv[argsi],"-mp")) {
				input = argv[++argsi];
				std::stringstream(input) >> mutationProbe;
			} else if (0 == strcmp(argv[argsi],"-mg")) {
				input = argv[++argsi];
				std::stringstream(input) >> maxGroupSize;
			} else if (0 == strcmp(argv[argsi],"-min")) {
				input = argv[++argsi];
				std::stringstream(input) >> minimize;
			} else if (0 == strcmp(argv[argsi],"-l")) {
				input = argv[++argsi];
				std::stringstream(input) >> byLearning;
			} else if (0 == strcmp(argv[argsi],"-e")) {
				input = argv[++argsi];
				std::stringstream(input) >> byEvolution;
			} else if (0 == strcmp(argv[argsi],"-d")) {
				input = argv[++argsi];
				std::stringstream(input) >> dimensions;
			} else if (0 == strcmp(argv[argsi],"-freq")) {
				input = argv[++argsi];
				std::stringstream(input) >> change_freq;
			} else if (0 == strcmp(argv[argsi],"-peaks")) {
				input = argv[++argsi];
				std::stringstream(input) >> peaks;
			} else if (0 == strcmp(argv[argsi],"-g")) {
				input = argv[++argsi];
				std::stringstream(input) >> generations;
			} else if (0 == strcmp(argv[argsi],"-grs")) {
				input = argv[++argsi];
				std::stringstream(input) >> groups;
			} 
		}
	}
}

#include "cmath"
int main(int argc, char const *argv[])
{
	parseCommandLine(argc, argv);

	init_peaks(dimensions, peaks, change_freq);

	specimen = createSpecimen(userrange, dimensions);
	Lem l(fitness, specimen, populationSize, byLearning, byEvolution, groups, maxGroupSize, minimize, classifierType);
	l.setParams(progressProbe, mutationProbe);

	for (unsigned i = 1; i < generations;++i) {
		l.step();
		l.evaluatePopulation(l.getPopulation(), dummyFitness);
		l.sortPopulation(l.getPopulation());
	}



	// Statistic s = l.getStatistic();
	// std::cout << "fitness" << '\n';
	// for (unsigned i = 0; i < s.bests.size(); ++i) {
	// 	std::cout << s.bests[i] << '\n';
	// }


	// std::cout << "OFFLINE ERROR\n";
	std::cout << get_offline_error() << '\n';
	return 0;
}