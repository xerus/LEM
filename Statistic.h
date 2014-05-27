/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * Statistic.h -- Definition of the structure for recording statistics
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
#ifndef STATISTIC_H
#define STATISTIC_H

#include <utility>
#include <vector>

using GroupStatistic = std::pair<double, std::pair<double, double>>;

struct Statistic
{
	std::vector<double> averages;
	std::vector<double> bests;
	std::pair<double, std::vector<std::vector<double>>> representations;
	std::vector<std::vector<GroupStatistic>> groupsStatistics;
};

#endif
