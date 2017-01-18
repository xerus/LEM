/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * movingSin.cpp -- LEM for a simple dynamic optimization problem
 * 
 ******************************************************************************
 *
 * Author: Pavel Grunt (xgrunt01@stud.fit.vutbr.cz)
 *        
 * Year: 2014
 * ----------------------------------------------------------------------------
 *
 * LEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LEM. If not, see <http://www.gnu.org/licenses/>.
 *
\*****************************************************************************/

#include "Lem.h"
#include <iostream>
#include <sstream>
#include <utility>
#include <cstring>
#include <cmath>

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
	const double PI = 3.141592653589793238462643383279502884;
	int dimensions = 8;
	int M = 10, T = 50; // T in {50,100,200}
	double k = 0.5; // k in {0.5, 4}
	std::pair<double, double> userrange({0, 80*PI});
	Specimen specimen;
	unsigned groups = 2;
	double maxGroupSize = 0.05, byLearning = 0.7, byEvolution = 0.2;
	unsigned populationSize = 100, generations = 1000, progressProbe = 4, mutationProbe = 80;
	unsigned minimize = 0;

	int G = 0;
	unsigned classifierType = CLSVAL;
	
	Specimen createSpecimen(std::pair<double, double> & range, unsigned dimensions) {
		Specimen s(dimensions);
		for (unsigned i = 0; i < dimensions; ++i)
		{
			s[i] = range;
		}
		return s;
	}

	double peakfunc(const double x, const int g) {
		double point = x - floor(g / T) * k * PI;
		if (point > 0 && point < PI) {
			return sin(point);
		}
		return 0;
	}

	double fitness(std::vector<double> &v) {
		double mean = 0.0;
		for (auto x : v) {
			mean += peakfunc(x, G);
		}		
		return mean / v.size();
	}

	void helpMsg() {
		generations = 0;
		std::cout << "LEM - Pohybliva sinusovka \n";
		std::cout << "Pocet opakovani  \t-M <unsigned>\n";
		std::cout << "Perioda zmeny    \t-T <unsigned>\n";
		std::cout << "Koeficient zmeny \t-k <double>\n";
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
			} else if (0 == strcmp(argv[argsi],"-g")) {
				input = argv[++argsi];
				std::stringstream(input) >> generations;
			} else if (0 == strcmp(argv[argsi],"-grs")) {
				input = argv[++argsi];
				std::stringstream(input) >> groups;
			} else if (0 == strcmp(argv[argsi],"-T")) {
				input = argv[++argsi];
				std::stringstream(input) >> T;
			} else if (0 == strcmp(argv[argsi],"-M")) {
				input = argv[++argsi];
				std::stringstream(input) >> M;
			} else if (0 == strcmp(argv[argsi],"-k")) {
				input = argv[++argsi];
				std::stringstream(input) >> k;
			} 
		}
	}
}


int main(int argc, char const *argv[])
{
	parseCommandLine(argc, argv);
	specimen = createSpecimen(userrange, dimensions);

#ifdef VERBOSE
	std::vector<double> stats(generations, 0);
#endif
	double EcSum = 0;

	for (int i = 0; i < M; ++i) {
		G = 0;
		Lem l(fitness, specimen, populationSize, byLearning, byEvolution, groups, maxGroupSize, minimize, classifierType);
		l.setParams(progressProbe, mutationProbe);
		for (int e = generations; G <= e; ) {
#ifdef VERBOSE			
			stats[G] += l.getPopulation().back().fitness;
#endif
			EcSum += l.getPopulation().back().fitness; // best of generation
			++G;
			l.evaluatePopulation(l.getPopulation(), fitness);
			l.sortPopulation(l.getPopulation());
			l.step();
		}
	}

#ifdef VERBOSE
	for (auto Fbg : stats) {
		std::cout << Fbg / M << '\n';
	}
#endif

	EcSum /= M;
	EcSum /= generations;

	std::cout << "Collective mean fitness\n" << EcSum << '\n';
	return 0;
}