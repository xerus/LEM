/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * Lem.cpp -- The implementation of LEM
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
#include "Statistic.h"
#include "Chromosome.h"

#include "SVM.h"
#include "AdaBoost.h"

#define MMIN(x, y) (((x) < (y)) ? (x) : (y))
#define MMAX(x, y) (((x) > (y)) ? (x) : (y))

#define GENE_TYPE double

#include <numeric>
#include <random>
#include <algorithm>

#include <memory>

#include <vector>
#include <utility>
#include <set>
#include <map>
#include <iostream>

namespace {
    std::random_device rd;
    std::default_random_engine r(rd());


    // 
    // Vraci dvojici (spodni hranice, horni hranice) pro zvolenou dimenzi
    //
    std::pair<GENE_TYPE, GENE_TYPE> & dimensionRange(Specimen & specimen, unsigned dim) {
        return specimen[dim];
    }

    //
    // Vraci maximum rozdilu hodnot ze vsech dimenzi
    //
    double getDiameter(std::vector<std::vector<double>> & examples) {
        double diameter = 0.001;
        double min, max;
        unsigned size = (examples.size() > 0) ? examples[0].size() : 0;
        for (unsigned di = 0; di < size; ++di) {
            min = std::numeric_limits<double>::max();
            max = std::numeric_limits<double>::lowest();
            for (auto & ex : examples) {
                min = std::min(min,ex[di]);
                max = std::max(max,ex[di]);
            }
            diameter = std::max(diameter, max - min);
        }
        return diameter;
    }


    template <typename TP>
    TP binom(TP n, TP k) {
        if (k == 1 || k == (n - 1)) return n;
        if (k == 0 || n == k) return 1;
        // k = (k < n - k) ? n - k : k;
        return binom(n - 1, k - 1) + binom(n - 1, k);
    }

    unsigned combinationCount(unsigned groups) {
        return binom(groups, (unsigned) 2);
    }



    //
    // Porovnani fitness jedincu
    //
    bool compareMaximizeChromosomes(const ChromosomePtr &ch1, const ChromosomePtr &ch2) {
        return ch1.fitness < ch2.fitness;
    }

    bool compareMinimizeChromosomes(const ChromosomePtr &ch1, const ChromosomePtr &ch2) {
        return ch1.fitness > ch2.fitness;
    }


    //
    // Generovani nahodnych cisel
    //
    //
    double normal_distribution(double mean, double stddev) {
        std::normal_distribution<double> distribution(mean, stddev);
        return distribution(r);
    }

    double rand_range(double lowerBound, double upperBound) {
        std::uniform_real_distribution<double> d(MMIN(lowerBound, upperBound),MMAX(lowerBound, upperBound));
        return d(r);            
    }

    template <typename TT>
    TT rand_range(TT lowerBound, TT upperBound) {
        std::uniform_int_distribution<TT> d(MMIN(lowerBound, upperBound),MMAX(lowerBound, upperBound));
        return d(r);            
    }

    template <typename TT>
    TT rand_range(std::pair<TT, TT> & bounds) {
        return rand_range(bounds.first, bounds.second);
    }

    int rand_int(int lowerBound, int upperBound) {
        std::uniform_int_distribution<int> d(MMIN(lowerBound, upperBound),MMAX(lowerBound, upperBound));
        return d(r);
    }

    int rand_int(int ub) {
        return rand_int(0, ub);
    }

    //
    // Kopirovani hodnot z jednoho jedince do druheho
    //
    void valueTransfer(const ChromosomePtr & from, ChromosomePtr & to) {
        to.fitness = from.fitness;
        std::copy(from.representation.begin(), from.representation.end(),
            to.representation.begin());
    }

    //
    // Vrati hodnotu pro zvolenou dimenzi nahodne vybraneho jedince z populace
    //
    double randomValueFromPopulation(Population & p, const int dimension) {
        int pos = rand_int(0, p.size() - 1);
        return p[pos].representation[dimension];
    }

    //
    // Vypocte prumernou hodnotu fitness v populaci
    //
    double sumPopulationFitness(Population & p) {
        double total = 0;
        for (auto & ch : p) {
            total += ch.fitness;
        }
        return total;
    }

    //
    // Vypocte prumernou hodnotu fitness v populaci
    //

    template <typename InputIterator1, typename InputIterator2>
    double averagePopulationFitness(InputIterator1 start, InputIterator2 end) {
        unsigned i = 0;
        double sum = 0;
        for (;start != end; ++start, ++i) {
            sum += (*start).fitness;
        }
        return sum / i;
    }

    //
    // Vypocte prumernou hodnotu fitness v populaci
    //
    double averagePopulationFitness(Population & p) {
        return averagePopulationFitness(p.begin(), p.end());
    }

    //
    // Zjisti nejlepsi jedince v populaci
    //
    std::vector<std::vector<double>> bestRepresentations(Population & p) {
        std::vector<std::vector<double>> representations;
        const double bestFitness = p.back().fitness;
        int i = p.size() - 1;
        while (i >= 0 && p[i].fitness == bestFitness) {
            representations.emplace_back((p[i].representation));
            --i;
        }
        return representations;
    }


    //
    // Pridani statistickych udaju o generaci
    // Vraci udaj, zda se zlepsila fitness
    //
    bool updateStatisticMaximize(Population & population, Statistic & stats,  std::pair<double, std::vector<std::vector<double>>> & absBest) {
        bool progress;
        double lastBest = stats.bests.back();
        stats.averages.emplace_back(averagePopulationFitness(population));
        stats.bests.emplace_back(population.back().fitness);

        progress = population.back().fitness > lastBest;
        if (progress && absBest.first < population.back().fitness) {
            absBest.first = population.back().fitness;
            absBest.second = bestRepresentations(population);
        }
        return progress;
    }

    bool updateStatisticMinimize(Population & population, Statistic & stats, std::pair<double, std::vector<std::vector<double>>> & absBest) {
        bool progress;
        double lastBest = stats.bests.back();
        stats.averages.emplace_back(averagePopulationFitness(population));
        stats.bests.emplace_back(population.back().fitness);

        progress = population.back().fitness < lastBest;
        if (progress && absBest.first > population.back().fitness) {
            absBest.first = population.back().fitness;
            absBest.second = bestRepresentations(population);
        }
        return progress;
    }


    // Single point crossover
    void crossover(const ChromosomePtr &ch1, const ChromosomePtr &ch2, ChromosomePtr & ch) {
        int point = rand_int(0, ch1.representation.size());
        for (int i = 0; i < point; ++i)
        {
            ch.representation[i] = ch1.representation[i];
        }
        for (unsigned i = point; i < ch2.representation.size(); ++i)
        {
            ch.representation[i] = ch2.representation[i];
        }
    }
    // uniform crossover
    void uniformCrossover(const ChromosomePtr &ch1, const ChromosomePtr &ch2, ChromosomePtr & ch) {
        for (unsigned i = 0; i < ch.representation.size(); ++i)
        {
            ch.representation[i] = (rand_range(0.,1.) < 0.5) ? ch1.representation[i] : ch2.representation[i];
        }
    }
    // Real value crossover
    void valueCrossover(const ChromosomePtr &ch1, const ChromosomePtr &ch2, ChromosomePtr & ch, Specimen & specimen) {
        for (unsigned i = 0; i < ch1.representation.size(); ++i)
        {
            do {
            double a = rand_range(-0.25, 1.25);
            ch.representation[i] = a * ch1.representation[i] + (1-a) * ch2.representation[i];
            } while (ch.representation[i] < dimensionRange(specimen,i).first 
                   || ch.representation[i] > dimensionRange(specimen,i).second) ;
        }
    }


    void randomChromosome(ChromosomePtr & ch, Specimen & specimen) {
        for (unsigned i = 0, e = specimen.size(); i < e; ++i ) {
            ch.representation[i] = rand_range(dimensionRange(specimen, i));
        }
    }

    // ocekava vzestupne serazene dle fitness population a children
    void rank_survival_Maximize(Population & population, Population & children) {
        unsigned popSize = population.size();
        for (unsigned i = 0, j = popSize - 1; i < popSize; ++i)
        {
            if (population[i].fitness < children[j].fitness) {
                valueTransfer(children[j], population[i]);
                --j;
            }
        }
    }
    void rank_survival_Minimize(Population & population, Population & children) {
        unsigned popSize = population.size();
        for (unsigned i = 0, j = popSize - 1; i < popSize; ++i)
        {
            if (population[i].fitness > children[j].fitness) {
                valueTransfer(children[j], population[i]);
                --j;
            }
        }
    }

    void mutateChromosome(ChromosomePtr & child, Specimen & s) {
        for (unsigned i = 0; i < s.size(); ++i) {
            double hunp = MMIN(1.0, 0.01 * (dimensionRange(s, i).second - dimensionRange(s,i).first));
            double chval = child.representation[i];
            child.representation[i] = normal_distribution(chval, hunp);
            while (child.representation[i] < dimensionRange(s, i).first 
                   || child.representation[i] > dimensionRange(s, i).second) {
                child.representation[i] = normal_distribution(chval, hunp);
            }
        }
    }

    void mutateChromosome(const ChromosomePtr & parent, ChromosomePtr & child, Specimen & s, double probability = .5) {
        for (unsigned i = 0; i < s.size(); ++i) {
            if (rand_range(0., 1.) >= probability) continue;
            double hunp = MMIN(1.0, 0.01 * (dimensionRange(s, i).second - dimensionRange(s, i).first));
            child.representation[i] = normal_distribution(parent.representation[i], hunp);
            while (child.representation[i] < dimensionRange(s, i).first 
                   || child.representation[i] > dimensionRange(s, i).second) {
                child.representation[i] = normal_distribution(parent.representation[i], hunp);
            }
        }
    }
  
    // 
    // Vytvoreni jednoho jedince dle SVM
    //
    void svmCrossover(Population &p, SVM & svm, ChromosomePtr & ch, double diameter, Specimen & s, double probability = .5) {
        bool found = false;
        int cnt = 100;
        if (svm.supportVectorsCount()) {
            const std::vector<double> & svref = svm.getSupportVector(rand_int(svm.supportVectorsCount()-1));
            for (unsigned di = 0; di < svref.size(); ++di) {
                if (rand_range(0., 1.) >= probability) ch.representation[di] = svref[di];
                else ch.representation[di] = rand_range((GENE_TYPE) MMAX(dimensionRange(s, di).first,  svref[di] - diameter), 
                                                   (GENE_TYPE) MMIN(dimensionRange(s, di).second, svref[di] + diameter));
            }
        } else {
            while(!found && cnt-- > 0) {
                mutateChromosome(p[rand_int(p.size()/2, p.size()-1)], ch, s);
                if (svm.classifyOne(ch.representation)) return;
            }
            cnt = 100;
            while (cnt-- > 100) {
                randomChromosome(ch, s);
                if (svm.classifyOne(ch.representation)) return;
            }
        }
    }
  
    // 
    // Vytvoreni jednoho jedince dle intervalu (dle adaBoost)
    // 
    void valueDimensionsCrossover(Population &p, RangesForClass & dimensions, Specimen & s, ChromosomePtr & ch) {
        for (int i = 0, e = s.size(); i < e; ++i)
        {
            if (dimensions.count(i)) {
                const unsigned rangeindex = rand_int(dimensions[i].size() - 1);
                unsigned ri = 0;
                for (auto & si : dimensions[i]) {
                    if (rangeindex == ri) {
                        ch.representation[i] = rand_range((GENE_TYPE) MMAX(dimensionRange(s,i).first, si.first), 
                                                          (GENE_TYPE) MMIN(dimensionRange(s,i).second, si.second));
                        break;
                    }
                    ++ri;
                }
            }
            else
                ch.representation[i] = randomValueFromPopulation(p, i);
        }
    }

    //
    // mutuje vsechny jedince ze pupulace p a vlozi je do populace children
    //
    void mutatePopulation(Population & p, Population & children, Specimen & s) {
        for (unsigned i = 0; i < p.size();) {
            mutateChromosome(p[i], children[i], s, 0.25);
            ++i;
        }
    }

    //
    // Vytvori pocatecni populaci
    //
    void createInitialPopulation(Population & p, Specimen & s) {
        for (unsigned i = 0; i < p.size(); )
        {
            randomChromosome(p[i], s);
            ++i;
        }
    }


    //
    // Pravdepodobnostni kolo
    //
    unsigned fitnessMaximizeWheel(const double & sumFitness, Population & p) {
        double target = rand_range(0.0, sumFitness);
        double tsum = 0;
        for (unsigned i = 0; i < p.size(); i++) {
            if (target <= tsum) {
                return i;
            }
            tsum += p[i].fitness;
        }
        return p.size() - 1;
    }   

    unsigned fitnessMinimizeWheel(const double & sumFitness, Population & p) {
        double target = rand_range(0.0, sumFitness);
        double tsum = 0;
        const double max = p[0].fitness;
        for (unsigned i = 0; i < p.size(); i++) {
            if (target <= tsum) {
                return i;
            }
            tsum += (max - p[i].fitness);
        }
        return p.size() - 1;
    }   
}



//
// Konstruktor
//
Lem::Lem(ObjectiveFuncType fitFunc, Specimen spc,
        unsigned pSize, double byL, double byE,
        const unsigned groups, const double maxGroupSize,
        unsigned minimize, unsigned classifierType) 
{   
    objFunc = fitFunc;
    fitnessFunc = &Lem::originalFitness;
    specimen = spc; 
    unsigned popSize = pSize;
    byL = (maxGroupSize * popSize < 1) ? 0 : std::min(std::max(byL, 0.0), 1.0);
    byE = std::min(std::max(byE, 0.0), 1.0 - byL);

    groupSize = std::min(maxGroupSize * popSize, (double) popSize / groups);
    groupCombinations = combinationCount(groups);

    childrenByGroup = (byL * popSize) / groupCombinations;

    byLearning = childrenByGroup * groupCombinations;
    byEvolution = byE * popSize;
    byRandom = popSize - byEvolution - byLearning;
    
    if (minimize == 1) {
        compareChromosomes = compareMinimizeChromosomes;
        fitnessWheel = fitnessMinimizeWheel;
        rank_survival = rank_survival_Minimize;
        updateStatistic = updateStatisticMinimize;
    }
    else  { // (minimize == 0)
        compareChromosomes = compareMaximizeChromosomes;
        fitnessWheel = fitnessMaximizeWheel;
        rank_survival = rank_survival_Maximize;
        updateStatistic = updateStatisticMaximize;
        if (minimize == 2) {
            fitnessFunc = &Lem::inverseFitness;
        } else if (minimize == 3) {
            fitnessFunc = &Lem::substractFitness;
        }
    }  

    generateByLearning = classifierType;

    population = Population(popSize);
    children = Population(popSize);
    for (unsigned i = 0; i < popSize; ++i) {
        population[i].representation = CHROMOSOME_REP(specimen.size());
        children[i].representation = CHROMOSOME_REP(specimen.size());
    }

    indexesForGroups((childrenByGroup > 0) ? groups : 0, popSize);
    createInitialPopulation(population, specimen);
    
    evaluatePopulation(population); 
    childrenStatistic(population, stats);  
    
    sortPopulation(population);

    TrainGroups.resize((childrenByGroup > 0) ? groups : 0);
    for (unsigned i = 0; i < TrainGroups.size(); ++i) {
        // TrainGroups[i].resize(childrenByGroup);
        TrainGroups[i].resize(groupSize);
        for (unsigned chi = 0; chi < TrainGroups[i].size(); ++chi) {
            TrainGroups[i][chi].resize(specimen.size());
        }
    }

    stats.averages.emplace_back(averagePopulationFitness(population));
    stats.bests.emplace_back(population.back().fitness);

    absBest = {population.back().fitness, bestRepresentations(population)};
    setParams();
}

//
// Jeden krok metody
//
void Lem::step() {
    if (!progress && --progressCounter <= 0) {    
        progressCounter = initProgressProbe;
        if (--mutationCounter > 0) {
            mutatePopulation(population, children, specimen);

            evaluatePopulation(children);
            childrenStatistic(children, stats);    
            sortPopulation(children);
            rank_survival(population, children);        
          } else 
        {
            mutationCounter = initMutationProbe;
            createInitialPopulation(population, specimen);
            evaluatePopulation(population);     
            childrenStatistic(population, stats);  
        }
    } else { // increase diversity
        if (progress)
            progressCounter = initProgressProbe;
        if (generateByLearning == 1 || (generateByLearning == 2 && generationNumber % 2 == 1))
            generateByLearningSVM();
        else 
            generateByLearningAdaBoost();
        generateByEvolution();
        generateRandomly();
        evaluatePopulation(children);
        childrenStatistic(children, stats);        
        sortPopulation(children);
        rank_survival(population, children);
    }
    sortPopulation(population);
    progress = updateStatistic(population, stats, absBest);
    ++generationNumber;
}

//
// Nastaveni nekterych ridicich parametru kroku 
//
void Lem::setParams(const int progressProbe, const int mutationProbe) {
    progressCounter = initProgressProbe = progressProbe;
    mutationCounter = initMutationProbe = mutationProbe;
}


//
// Nastaveni poctu kroku k provedení
//
void Lem::run(const unsigned generations, int progressProbe, int mutationProbe) {
    setParams(progressProbe, mutationProbe);
    const unsigned newcapacity = stats.averages.size() + generations;
    stats.averages.reserve(newcapacity);
    stats.bests.reserve(newcapacity);
    stats.groupsStatistics.reserve(newcapacity);
    for (unsigned i = 0; i < generations; ++i)
    {
        step();
    }
}


//
// Generovani jedincu dle AdaBoost klasifikatoru
//
void Lem::generateByLearningAdaBoost() {
    extremaSelection(population);
    for (unsigned i = 0, chi = 0; i < TrainGroups.size(); ++i) {
        for (unsigned j = i + 1; j < TrainGroups.size(); ++j) {
            AdaBoost a(TrainGroups[i], TrainGroups[j], 2 * specimen.size());
            RangesForClass ranges = a.rangesForClass(true);
            for(unsigned dd = 0; dd < childrenByGroup;) {
                valueDimensionsCrossover(population, ranges, specimen, children[chi]);
                {++chi;++dd;}
            }
        }
    }
}

//
// Generovani jedincu dle SVM klasifikatoru
//
void Lem::generateByLearningSVM() {
    extremaSelection(population);
    for (unsigned i = 0, chi = 0; i < TrainGroups.size(); ++i) {
        for (unsigned j = i + 1; j < TrainGroups.size(); ++j) {
            const double diameter = getDiameter(TrainGroups[i]);
            SVM svm(TrainGroups[i], TrainGroups[j], diameter);
            for(unsigned dd = 0; dd < childrenByGroup;) {
                svmCrossover(population, svm, children[chi], diameter, specimen);
                {++chi;++dd;}
            }
        }
    }
}    

//
// Vygeneruje jedenci evolucnimi operatory
//
void Lem::generateByEvolution() {
    double sumFitness = sumPopulationFitness(population);
    for (unsigned i = byLearning, mod = 0, e = byLearning + byEvolution; i < e; ++mod)
    {
        unsigned p1 = fitnessWheel(sumFitness, population);
        unsigned p2 = fitnessWheel(sumFitness, population);
        switch (mod % 3) {
        case 0:
            uniformCrossover(population[p1], population[p2], children[i]);
            break;
        case 1:
            valueCrossover(population[p1], population[p2], children[i], specimen);
            break;
        default:
            mutateChromosome(population[p1], children[i], specimen);
        }
        ++i;
    }       
}

//
// Nahodne vygeneruje jedince
//
void Lem::generateRandomly() {
    for (unsigned i = population.size() - byRandom; i < population.size(); ++i) {
        randomChromosome(children[i], specimen);
    }
}


// 
// Prace s populaci
//
void Lem::sortPopulation(Population &p) {
    std::sort(p.begin(), p.end(), compareChromosomes);
}

void Lem::evaluatePopulation(Population & p, ObjectiveFuncType ff) {
    for (auto & ch : p) {
        ch.fitness = ff(ch.representation);
    }
}


void Lem::evaluatePopulation(Population & p) {
    for (auto & ch : p) {
        ch.fitness = (*this.*fitnessFunc)(ch.representation);
    }
}


void Lem::indexesForGroups(unsigned groups, unsigned popSize) {
    indexes.resize(groups);
    unsigned up, down;

    down = 0;
    up = popSize - groupSize;

    for (unsigned i = 0; i < groups; ++i) {
        if (i % 2 == 0) {
            indexes[i]=(up);
            up -= groupSize;
        } else {
            indexes[i]=(down);
            down += groupSize;
        }
    }

    std::sort(indexes.begin(), indexes.end());
}

void Lem::extremaSelection(Population & p) {
    for (unsigned i = 0; i < indexes.size(); ++i) {     
        for (unsigned ai = 0, j = indexes[i], e = indexes[i] + groupSize; j < e; ++j, ++ai) {
            std::copy(p[j].representation.begin(), p[j].representation.end(),
                TrainGroups[TrainGroups.size() - i - 1][ai].begin());
        }           
    }
}

//
// Prevedeni minimalizace na maximalizaci
// f'(X) = 1 / (1 + f(X))
//
inline double Lem::inverseFitness(std::vector<double> & representation) {
    return 1.0 / (1 + objFunc(representation));
}


//
// Prevedeni minimalizace na maximalizaci
// f'(X) = 10000000 - f(X)
//
inline double Lem::substractFitness(std::vector<double> & representation) {
    return 10000000 - objFunc(representation);
}

//
// f'(X) =  f(X)
//
inline double Lem::originalFitness(std::vector<double> & representation) {
    return objFunc(representation);
}


//
// Vraci udaje o behu -- nejlepsi dosazene hodnoty
//
Statistic Lem::getStatistic() {
    stats.representations = absBest;
    return stats;
}

//
// Vypocet statistic pro jednotlive kombinace skupin
//
void Lem::computeGroupsStatistics(std::vector<GroupStatistic> & gs, Population & p) {
    double groupAvrg;
    for (unsigned i = 0; i < groupCombinations; ++i) {
        groupAvrg = averagePopulationFitness(p.begin() + i * childrenByGroup,
                                             p.begin() + (i+1) * childrenByGroup);
        gs[i] = {groupAvrg, {(*std::min_element(p.begin() + i * childrenByGroup,
                                                p.begin() + (i+1) * childrenByGroup, 
                                                compareChromosomes)).fitness,
                             (*std::max_element(p.begin() + i * childrenByGroup,
                                                p.begin() + (i+1) * childrenByGroup, 
                                                compareChromosomes)).fitness}};
    }
    // by Evolution
    if (byEvolution > 0) {
        groupAvrg = averagePopulationFitness(p.begin() + byLearning,
                                             p.begin() + byLearning + byEvolution);
        gs[groupCombinations] = 
                {groupAvrg, {(*std::min_element(p.begin() + byLearning,
                                                p.begin() + byLearning + byEvolution, 
                                                compareChromosomes)).fitness,
                             (*std::max_element(p.begin() + byLearning,
                                                p.begin() + byLearning + byEvolution, 
                                                compareChromosomes)).fitness}}; 
    }
    // by Random
    if (byRandom > 0) {         
        groupAvrg = averagePopulationFitness(p.begin() + byLearning + byEvolution,
                                             p.end());
        gs[groupCombinations+1] = 
                {groupAvrg, {(*std::min_element(p.begin() + byLearning + byEvolution,
                                                p.end(), 
                                                compareChromosomes)).fitness,
                             (*std::max_element(p.begin() + byLearning + byEvolution,
                                                p.end(), 
                                                compareChromosomes)).fitness}};     
    }
}

//
// Zada nejlepsi a prumernou hodnotu fitness v populaci a pro jednotlive kombinace skupin
//
void Lem::childrenStatistic(Population & p, Statistic & stats) {
    std::vector<GroupStatistic> gs(groupCombinations+2);
    computeGroupsStatistics(gs, p);
    stats.groupsStatistics.emplace_back(gs);
}




//
// Vraci referenci na populaci
//
std::vector<Chromosome<std::vector<double>>> & Lem::getPopulation() {
    return population;
}

