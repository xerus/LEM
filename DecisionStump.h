/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * DecisionStump.h -- Implementation of the Decision Stump Classifier
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
#ifndef DECISION_STUMP_H
#define DECISION_STUMP_H 

#include <vector>
class DecisionStump
{
friend class AdaBoost;	
public:
	DecisionStump() : prepared(false) {};
	template <typename T>
	DecisionStump(T threshold, int dimension, bool mainclass) :
		threshold(threshold), dimension(dimension), mainclass(mainclass), prepared(true) {}
	template <typename T>
	std::vector<bool> classify(std::vector<T> & data);
	template<class T> 
	void classify(std::vector<T> & data, std::vector<bool> & classified);
	template <typename T>
	bool classifyOne(T & data);
private:
	void setup(double threshold, int dimension, bool mainclass);
	double threshold;
	int dimension;
	bool mainclass;
	bool prepared;
};

template<class T> 
std::vector<bool> DecisionStump::classify(std::vector<T> & data) {
	std::vector<bool> classified(data.size());
	if (mainclass) {
		for (unsigned i = 0; i < data.size(); ++i)
		{
			classified[i] = data[i][dimension] >= threshold;
		}
	} else {
		for (unsigned i = 0; i < data.size(); ++i)
		{
			classified[i] = data[i][dimension] < threshold;
		}
	}
	return classified;
}

template<class T> 
void DecisionStump::classify(std::vector<T> & data, std::vector<bool> & classified) {
	if (mainclass) {
		for (unsigned i = 0; i < data.size(); ++i)
		{
			classified[i] = data[i][dimension] >= threshold;
		}
	} else {
		for (unsigned i = 0; i < data.size(); ++i)
		{
			classified[i] = data[i][dimension] < threshold;
		}
	}
}


template <typename T>
bool DecisionStump::classifyOne(T & data) {
	return (mainclass) ? data[dimension] >= threshold : data[dimension] < threshold;
}

#endif