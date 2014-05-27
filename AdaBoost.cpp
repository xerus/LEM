/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * AdaBoost.cpp -- AdaBoost implementation
 * 
 ******************************************************************************
 *
 * Author: Pavel Grunt (xgrunt01@stud.fit.vutbr.cz)
 *        
 * Year: 2014
 * ----------------------------------------------------------------------------
 *
 * THIS SOFTWARE IS NOT COPYRIGHTED
 *
 * This source code is offered for use in the public domain.
 * You may use, modify or distribute it freely.
 *
 * This source code is distributed in the hope that it will be useful but
 * WITHOUT ANY WARRANTY.  ALL WARRANTIES, EXPRESS OR IMPLIED ARE HEREBY
 * DISCLAIMED.  This includes but is not limited to warranties of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
\*****************************************************************************/
#include "AdaBoost.h"
#include <set>
#include <map>
#include <utility>

#include <algorithm> 
#include <limits>

#include <iostream>

namespace {
	const double absmin = std::numeric_limits<double>::lowest();
	const double absmax = std::numeric_limits<double>::max();

	template <typename V>
	bool hasKey(std::map<unsigned, V> & m, unsigned key) {
		return m.count(key);
	}

	bool inRange(const std::pair<double, double> & range, double number) {
		return number >= range.first && number <= range.second;
	}

	std::pair<double, double> rangeForStump(const double stumpthreshold, 
		const bool stumpclass, const bool mainclass, const double min = absmin, const double max = absmax) 
	{
		if (stumpclass ^ mainclass) {
			return std::pair<double, double>{min, stumpthreshold};
		} else {			
			return std::pair<double, double>{stumpthreshold, max};
		}
	}

	bool intersection(const std::pair<double, double> & fst, const std::pair<double, double> & snd,
		std::pair<double, double> & result) 
	{
		if ((fst.first == snd.first) || (fst.second == snd.second)) {
			return false;
		} else {
			result.first = std::max(fst.first, snd.first);
			result.second = std::min(fst.second, snd.second);
		}
		return result.first < result.second;
	}

	bool isLowerThan(std::pair<double, double> & newSet, std::pair<double, double> & oldSet) 
	{
		if (newSet.first == oldSet.first) {
			return newSet.second < oldSet.second;
		} else if (newSet.second == oldSet.second) {
			return newSet.first > oldSet.first;
		} else {
			return (newSet.second - newSet.first) < (oldSet.second - oldSet.first);
		}
	}

	std::set<std::pair<double, double>> mergeRanges(std::set<std::pair<double, double>> & ranges) {
		std::set<std::pair<double, double>> merged;
		for (auto & ri : ranges) {
			std::pair<double, double> minSet(ri);
			for (auto & rj : ranges) {
				if (ri == rj) continue;
				std::pair<double, double> tmpSet;
				if (intersection(ri, rj, tmpSet) && isLowerThan(tmpSet, minSet)) {
					minSet = tmpSet;
				}
			}
			merged.insert(minSet);
		}
		return merged;
	}

	void printRanges(RangesForClass & r) {
		std::cerr << "Range size: " << r.size() << '\n';
		for (auto & di : r) {
			std::cerr << "dimension " << di.first << ":\t";
			for (auto & ri : di.second) {
				std::cerr << '(' << ri.first << ',' << ri.second << ")\t";
			}
			std::cerr << '\n';
		}
	}
}


std::set<int> AdaBoost::usedDimensions() {
	std::set<int> dimensions;
	for (unsigned i = 0; i < alphas.size(); ++i)
	{
		dimensions.insert(stumps[i].dimension);
	}
	return dimensions;
}

RangesForClass AdaBoost::rangesForClass(bool mainclass) {
	RangesForClass ranges;
	for (unsigned i = 0; i < alphas.size(); ++i) {
		unsigned key = stumps[i].dimension;
		std::pair<double, double> range = rangeForStump(stumps[i].threshold, stumps[i].mainclass, mainclass);
		if (hasKey(ranges, key))
			ranges[key].insert(range);
		else
			ranges[key] = {range};
	}

	for (auto & ri : ranges) {
		ri.second = mergeRanges(ri.second);
	}
	return ranges;
}

