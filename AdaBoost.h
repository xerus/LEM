/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * AdaBoost.h -- Implementation of AdaBoost Classifier
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
#ifndef ADABOOST_H
#define ADABOOST_H

#include <limits>
#include <cmath>

#include <utility>
#include <vector>
#include <set>
#include <map>

#include <memory>
#include <algorithm>

#include "DecisionStump.h"

using RangesForClass = std::map<unsigned, std::set<std::pair<double, double>>>;

class AdaBoost
{
public:
	template <typename T>
	AdaBoost(std::vector<std::pair<std::vector<T>, bool>> & data, int maxClassifiers=10);
	template <typename T>
	AdaBoost(std::vector<std::vector<T>> & dataonly, std::vector<bool> & labelonly, int maxClassifiers=10);
	template <typename T>
	AdaBoost(std::vector<std::vector<T>> & positives, std::vector<std::vector<T>> & negatives, int maxClassifiers=10);
	template <typename T>
	std::vector<bool> classify(std::vector<T> & data);
	template <typename T>
	bool classifyOne(T & data);	
	std::set<int> usedDimensions();
	RangesForClass rangesForClass(bool mainclass);
private:
	template <typename T>
	void selfInitialize(std::vector<std::vector<T>> & dataonly, std::vector<bool> & labelonly, int maxClassifiers);
	std::vector<double> alphas;
	std::vector<DecisionStump> stumps;
};

template <typename T>
void AdaBoost::selfInitialize(std::vector<std::vector<T>> & dataonly, std::vector<bool> & labelonly, int maxClassifiers) {
	alphas.reserve(maxClassifiers);
	std::vector<double> d(dataonly.size(), 1.0 / dataonly.size());
	stumps = std::vector<DecisionStump>(maxClassifiers);
	unsigned dimensions =  dataonly[0].size();

	std::vector<std::set<T>> dimThresholds(dimensions);
	std::vector<std::set<T>> meanvalSets(dimensions);
	for (unsigned ri = 0, e = dataonly.size(); ri < e; ++ri)
	{
		for (unsigned i = 0; i < dimensions; ++i)
		{
			dimThresholds[i].insert(dataonly[ri][i]);
		}
	}

	T setval;
	for (unsigned i = 0; i < dimensions; ++i) {
		bool setptr = false;
		if (dimThresholds[i].size() > 1) {
			for (auto el : dimThresholds[i]) {
				if(setptr) {
					T meanval = (setval + el) / 2.0;
					meanvalSets[i].insert(meanval);
					setval = el;
				} else {
					setptr = true;
					setval = el;
				}
			}
		} else {
			for (auto el : dimThresholds[i]) {
				meanvalSets[i].insert(el);
			}
		}
	}

	std::vector<bool> classified(dataonly.size());
	std::vector<bool> tclassified(dataonly.size());
	for (int classifier = 0; classifier < maxClassifiers; ++classifier)
	{
		double err = std::numeric_limits<double>::infinity();
		for (unsigned dim = 0; dim < dimensions; ++dim)
		{	
			for (bool mainclass : { false, true }) 
			{
				for (auto threshold : meanvalSets[dim]) {
					double oldthrs = threshold;
					int olddim = dim;
					bool oldcls = mainclass;
					if (stumps[classifier].prepared) {
						oldthrs = stumps[classifier].threshold;
						olddim = stumps[classifier].dimension;
						oldcls = stumps[classifier].mainclass;
						stumps[classifier].setup(threshold, dim, mainclass);
					} else {
						stumps[classifier] = DecisionStump(threshold, dim, mainclass);						
					}

					stumps[classifier].classify(dataonly, classified);
					double tmp_err = 0;
					for (unsigned i = 0; i < dataonly.size(); ++i)
					{
						if (labelonly[i] != classified[i]) tmp_err += d[i];
					}
					if (tmp_err < err) {
						err = tmp_err;
						std::copy(classified.begin(), classified.end(), tclassified.begin());
					} else {
						stumps[classifier].setup(oldthrs, olddim, oldcls);
					}
				}
			}	
		}
		if (err >= 0.5) break; // not sufficient

		// best classifier was found
		if (err == 0) {
			alphas.emplace_back(1);
			break;
		} else {
			alphas.emplace_back(0.5 * log10((1 - err) / err));
		}

		double sumD = 0;
		for (unsigned i = 0; i < dataonly.size(); ++i)
		{
			if (labelonly[i] == tclassified[i]) d[i] = d[i] * ((1 - err) / err);
			sumD += d[i];
		}
		// normalize
		for (auto & di : d) {
			di /= sumD;
		}
	}
}

template <typename T>
AdaBoost::AdaBoost(std::vector<std::vector<T>> & dataonly, std::vector<bool> & labelonly, int maxClassifiers) {
	selfInitialize(dataonly, labelonly, maxClassifiers);
}

template <typename T>
AdaBoost::AdaBoost(std::vector<std::pair<std::vector<T>, bool>> & traindata, int maxClassifiers) 
{
	std::vector<std::vector<T>> dataonly(traindata.size());
	std::vector<bool> labelonly(traindata.size());
	for (unsigned i = 0; i < traindata.size(); ++i)
	{
		dataonly[i] = traindata[i].first;
		labelonly[i] = traindata[i].second;
	}
	selfInitialize(dataonly, labelonly, maxClassifiers);
}


template <typename T>
AdaBoost::AdaBoost(std::vector<std::vector<T>> & positives, std::vector<std::vector<T>> & negatives, int maxClassifiers)
{
	std::vector<std::vector<T>> dataonly(positives.size() + negatives.size());
	std::vector<bool> labels(positives.size() + negatives.size(), true);
	
	std::copy(positives.begin(), positives.end(), dataonly.begin());
	std::copy(negatives.begin(), negatives.end(), dataonly.begin() + positives.size());

	std::fill(labels.begin() + positives.size(), labels.end(), false);

	selfInitialize(dataonly, labels, maxClassifiers);
}


template <typename T>
std::vector<bool> AdaBoost::classify(std::vector<T> & data) {
	std::vector<bool> allclassified(data.size(), false);
	for (unsigned i = 0; i < data.size(); ++i) {
		allclassified[i] = classifyOne(data[i]);
	}
	return allclassified;
}

template <typename T>
bool AdaBoost::classifyOne(T & data) {
	double sum = 0;
	for (unsigned ai = 0; ai < alphas.size(); ++ai) {
		sum += (stumps[ai].classifyOne(data)) ? alphas[ai] : -alphas[ai];
	}
	return sum >= 0;
}


#endif
