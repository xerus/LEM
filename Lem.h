/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * Lem.h -- Definition of LEM class
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

#ifndef LEM_H
#define LEM_H

#include "Statistic.h"
#include "Chromosome.h"
#include <vector>
#include <utility>

using ObjectiveFuncType = double (*) (std::vector<double> &);
using Specimen = std::vector<std::pair<double, double>>;

using CHROMOSOME_REP = std::vector<double>;
using ChromosomePtr = Chromosome<CHROMOSOME_REP>;
using Population = std::vector<ChromosomePtr>;
using Group = std::vector<CHROMOSOME_REP>;
using ChromosomeComparsion = bool (*) (const ChromosomePtr &, const ChromosomePtr &);
using FitnessWheelType = unsigned (*) (const double & , Population & );
using VoidFuncType = void (*) ();
using RankSurvType = void (*) (Population &, Population &);
using StatisticFuncType = bool (*) (Population & , Statistic & ,  std::pair<double, std::vector<std::vector<double>>> & );

class Lem
{
public:
	//
	// Konstruktor - parametry objektivni funkce, vzor jedince
	// classifierType: 1 - SVM
	//                 2 - Stridani AdaBoost a SVM
	//             jinak - AdaBoost
	Lem(ObjectiveFuncType f, Specimen specimen,
	 	const unsigned popSize = 100,
	 	double byLearning = 0.70, double byEvolution = 0.20,
	 	const unsigned groups = 2, const double maxGroupSize = 0.30,
	 	unsigned minimize = 0, unsigned classifierType = 0);

	void step();
	void run(const unsigned steps, int progressProbe = 4, int mutationProbe = 10);

	std::vector<Chromosome<std::vector<double>>> & getPopulation();
	void sortPopulation(Population &p);
	void evaluatePopulation(Population & p, ObjectiveFuncType ff);
	void evaluatePopulation(Population & p);
	void setParams(const int progressProbe = 4, const int mutationProbe = 10);

	Statistic getStatistic();
private:
    unsigned generationNumber;
    unsigned classifierCount;    
    std::vector<unsigned> indexes;
    std::vector<Group> TrainGroups;

    Population population, children;

    Statistic stats;
    std::pair<double, std::vector<std::vector<double>>> absBest;

    int progressCounter, initProgressProbe;
    int mutationCounter, initMutationProbe;

    Specimen specimen;
    unsigned byLearning, byEvolution, byRandom;
    unsigned groupCombinations, childrenByGroup, groupSize;
    bool progress;
    // ukazatele na funkce
    double (Lem::*fitnessFunc)(std::vector<double>&);
    ObjectiveFuncType objFunc;
    unsigned generateByLearning;
    ChromosomeComparsion compareChromosomes;
    FitnessWheelType fitnessWheel;
    RankSurvType rank_survival;
    StatisticFuncType updateStatistic;

    //
    // funkce pro vytvareni populace
    //
    void generateByLearningAdaBoost();
    void generateByLearningSVM();
    void generateByEvolution();
    void generateRandomly();
    // transformace fitness funkce;
    inline double inverseFitness(std::vector<double> & representation);
    inline double substractFitness(std::vector<double> & representation);
    inline double originalFitness(std::vector<double> & representation);
    // statistiky
    void computeGroupsStatistics(std::vector<GroupStatistic> & gs, Population & p);
    void childrenStatistic(Population & p, Statistic & stat);
    //
    void indexesForGroups(unsigned groups, unsigned popSize);
    void extremaSelection(Population & p);

};


#endif
