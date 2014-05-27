/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * DecisionStump.cpp -- Implementation of the Decision Stump Classifier
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
#include "DecisionStump.h"

void DecisionStump::setup(double threshold, int dimension, bool mainclass) {
	this->threshold = threshold;
	this->dimension = dimension;
	this->mainclass = mainclass;
	prepared = true;
}


