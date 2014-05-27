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