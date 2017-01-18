/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * main.cpp -- LEM for function optimization
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
#include <iostream>
#include <sstream>
#include <memory>


#include <vector>
#include <utility>
#include <algorithm>

#include <random>
#include <numeric>

#include <cmath>
#include <cstring>

#include "Lem.h"
#include "Statistic.h"

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

// AdaBoost
#ifndef CLSVAL
#define CLSVAL 0
#endif

namespace {
	std::pair<double, double> sphererange({-600, 600});
	std::pair<double, double> rosenbrockrange({-2.0, 10});
	std::pair<double, double> griewangkrange({-5.12, 5.11});
	std::pair<double, double> ackleyrange({-32.768, 32.768});

	const double PI = 3.141592653589793;

	std::pair<double, double> userrange({-600, 600});
	ObjectiveFuncType fitnessFunc = nullptr;
	unsigned groups = 2;
	double maxGroupSize = 0.05, byLearning = 0.7, byEvolution = 0.2;
	unsigned populationSize = 100, dimensions = 10, generations = 100, progressProbe = 4, mutationProbe = 80;
	unsigned minimize = 2;
	unsigned classifierType = CLSVAL;


	Specimen createSpecimen(std::pair<double, double> & range, unsigned dimensions) {
		Specimen s(dimensions);
		for (unsigned i = 0; i < dimensions; ++i)
		{
			s[i] = range;
		}
		return s;
	}

	Specimen specimen;

	static std::random_device rd;
	static std::default_random_engine r(rd());
	double rand_range(double lowerBound, double upperBound) {
		std::uniform_real_distribution<double> d(lowerBound, upperBound);
		return d(r);
	}
	double normal_distribution(double mean, double stddev) {
		std::normal_distribution<double> distribution(mean, stddev);
		return distribution(r);
	}


	// sphere
	double sphere(std::vector<double> & rep) {
		return std::inner_product(rep.begin(), rep.end(),
								  rep.begin(), 0.0);
	}

	// rosenbrock
	double rosenbrock(std::vector<double> & rep) {
		double total = 0.0;
		for (unsigned i = 0, e = rep.size() - 1; i < e; ++i) {
			double rep2 = std::pow(rep[i], 2);
			total += 100.0 * std::pow(rep[i+1] - rep2, 2) + std::pow(rep[i] - 1.0, 2);
		}
		return total;
	}	

	// griewangk
	double griewangk(std::vector<double> & rep) {
		double cosprod = 1;
		double total = 0;
		for (unsigned i = 0, e = rep.size(); i < e; i++) {
			cosprod *= cos(rep[i] / sqrt(i+1));
			total += std::pow(rep[i], 2) / 4000.0;
		}
		total = 1 + total - cosprod;
		return total;
	}

	double ackley(std::vector<double> & rep) {
		double sqsum = 0;
		double cossum = 0;
		double total;
		for (unsigned i = 0, e = rep.size(); i < e; i++) {
			sqsum += std::pow(rep[i], 2);
			cossum += cos(2*PI*rep[i]);
		}
		sqsum = sqsum / rep.size();
		cossum = cossum / rep.size();
		total = 20 + exp(1) - 20*exp(-0.2 * sqrt(sqsum)) - exp(cossum);
		return total;
	}

	double rastrigin(std::vector<double> & rep) {
		double cossum = 0;
		for (unsigned i = 0, e = rep.size(); i < e; i++) {
			cossum += std::pow(rep[i], 2) - 10 * cos(2*PI*rep[i]);
		}
		return 10 * rep.size() + cossum;
	}

	double schwefel(std::vector<double> & rep) {
		return std::accumulate(rep.begin(), rep.end(), 418.9828872724338 * rep.size(),
			 [](double part, double val) {return part - val*sin(sqrt(fabs(val)));});
	}

	double step(std::vector<double> & rep) {
		return std::accumulate(rep.begin(), rep.end(), 0,
			 [](double part, double val) {return part + floor(val);});
	}
	double funcvalue(double fitness) {
		return (1.0 / fitness) - 1.;
	}

	void helpMsg() {
		generations = 0;
		std::cout << "LEM - hledani minima funkci \n" << "Funkce \t-f {ackley, griewangk, rastrigin, rosenbrock, schwefel, sphere}\n";
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
			} else if (0 == strcmp(argv[argsi],"-f")) {
				input = argv[++argsi];
				if ("sphere" == input) {
					fitnessFunc = sphere;
				} else if ("rosenbrock" == input) {
					fitnessFunc = rosenbrock;	
				} else if ("ackley" == input) {
					fitnessFunc = ackley;	
				} else if ("griewangk" == input) {
					fitnessFunc = griewangk;	
				} else if ("schwefel" == input) {
					fitnessFunc = schwefel;
				} else if ("step" == input) {
					fitnessFunc = step;
				} else if ("rastrigin" == input) {
					fitnessFunc = rastrigin;
				}
			}
		}
		if (fitnessFunc == nullptr) {
			fitnessFunc = sphere;
		}
	}

	inline double functionValue(double v) {
		if (minimize == 2) return (1.0 / v) - 1; // inverseFitness
		if (minimize == 3) return 10000000 - v;
		return v;
	}
}

int main(int argc, char const *argv[])
{
	parseCommandLine(argc, argv);
	specimen =  createSpecimen(userrange, dimensions);
	Lem l(fitnessFunc, specimen, populationSize, byLearning, byEvolution, groups, maxGroupSize, minimize, classifierType);


	l.run(generations, progressProbe, mutationProbe);
	Statistic s = l.getStatistic();

	std::cout << "fitness" << ';' << "value" << "\n";
	for (unsigned i = 0; i < s.bests.size(); ++i) {
		// double av = s.averages[i];
		// std::cout << r << ';' <<  funcvalue(r) << ';' << av << ';' << funcvalue(av) << '\n';
		std::cout << functionValue(s.bests[i]) << ';' <<  functionValue(s.averages[i]) << ";\t";
		// for (auto & pp : s.groupsStatistics[i]) {
		// 	std::cout << pp.first << ';' << pp.second.first << ';' << pp.second.second << ";\t";
		// }
		std::cout << '\n';
	}
	// std::cout << "SIZE" << s.representations.second.size();
	
	// std::cout << "# BEST INDIVIDUALS \n";
	// for (auto & ch : s.representations.second) {
	// 	std::cout << "# [" ;
	// 	for (double chi : ch)
	// 	{
	// 		std::cout << chi << ", ";
	// 	}
	// 	std::cout << "]\t";
	// 	std::cout << fitnessFunc(ch) << ';' << funcvalue(fitnessFunc(ch)) << '\n';
	// }
	// std::vector<double> vv = {420.968746,420.968746};
	// std::cout << fitnessFunc(vv) << '\n';
	// std::cout << "#BEST FITNESS; FUNCTION VALUE\n";
	// std::cout << '#' << s.representations.first << ';' << funcvalue(s.representations.first) << '\n';

	// std::cout << s.bests.size() << " " << s.groupsStatistics.size() << '\n';

	return 0;
}