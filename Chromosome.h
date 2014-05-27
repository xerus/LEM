/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * Chromosome.h -- Definition of the Chromosome (Individual)
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

#ifndef CHROMOSOME_H
#define CHROMOSOME_H

template <typename T>
struct Chromosome
{
	T representation;
	double fitness;
};
#endif