/*****************************************************************************\
 *
 * Learnable Evolution Model for Optimization (LEM)
 *
 * SVM.h -- Definition of the SVM Classifier class
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
#ifndef SVM_H
#define SVM_H

#include <vector>
#include <set>

class SVM
{
public:
    SVM(std::vector<std::vector<double>> & positives, std::vector<std::vector<double>> & negatives, double rbfconst = -123.45, double C=10);
    bool classifyOne(const std::vector<double> & x);
    const std::vector<double> & getSupportVector(unsigned i, bool positives = true);
    unsigned supportVectorsCount(bool positives = true);
private:
    double SMO_Error(const std::vector<double> & x1, const int y1);
    double SMO_Error(const unsigned int index);
    bool SMO_ViolatesKKT(const double alpha, const int y, const double E);
    bool SMO_OptimizationStep(double & alpha1, const std::vector<double> & x1, const int y1, double E1, 
                          double & alpha2, const std::vector<double> & x2, const int y2, double E2);
    void SMO_Learning();
    bool SMO_OptimizeOne(double & alpha1, const std::vector<double> & x1, const int y1, double E1);
    bool SMO_OptimizeOne(const unsigned int index);
    const std::vector<double> & getExample(const unsigned i);

    void computeErrors();
    double calculateOutput(const std::vector<double> & x1);

    std::vector<std::vector<double>> & positives;
    std::vector<std::vector<double>> & negatives;
    std::vector<double> negAlphas;
    std::vector<double> posAlphas;
    std::vector<double> errors;
    std::set<unsigned> posSkip;
    std::set<unsigned> negSkip;
    double rbfconst, b, C;
};



#endif