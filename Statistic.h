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
