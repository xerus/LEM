#
# LEM - Makefile
# Pavel Grunt, xgrunt01@stud.fit.vutbr.cz
# 2014
#

#
# eva.fit.vutbr.cz
# CXX=g++47
#
# merlin.fit.vutbr.cz
# CXX=g++-4.7
#
CXXFLAGS=-std=c++11 -O3
OUTPUT=lem

all: lem $(OUTPUT)AdaBoost $(OUTPUT)SVM $(OUTPUT)Combined movingSinA movingPeaksA

movingPeaksA: movingPeaksAdaBoost movingPeaksSVM movingPeaksCombined

movingSinA: movingSinAdaBoost movingSinSVM movingSinCombined

movingSinAdaBoost: Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DADABOOST -o $@ Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 

movingSinSVM: Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DSVM -o $@ Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 

movingSinCombined: Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DCOMBINED -o $@ Lem.o movingSin.cpp SVM.o DecisionStump.o AdaBoost.o 

movingPeaksAdaBoost: Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o
	$(CXX) $(CXXFLAGS) -DADABOOST -o $@ Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o

movingPeaksSVM: Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o
	$(CXX) $(CXXFLAGS) -DSVM -o $@ Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o

movingPeaksCombined: Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o
	$(CXX) $(CXXFLAGS) -DCOMBINED -o $@ Lem.o movingPeaks.cpp SVM.o DecisionStump.o AdaBoost.o movpeaks.o

lemAdaBoost: Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DADABOOST -o $@ Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 

lemSVM: Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DSVM -o $@ Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 

lemCombined: Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 
	$(CXX) $(CXXFLAGS) -DCOMBINED -o $@ Lem.o main.cpp SVM.o DecisionStump.o AdaBoost.o 


lem: lemCombined
	rm -rf $@
	ln -s lemCombined $@

pack: clean
	zip lem.zip *.cpp *.h Makefile

clean:
	rm -rf $(OUTPUT)AdaBoost $(OUTPUT)SVM $(OUTPUT)Combined  movingSinAdaBoost movingSinSVM movingSinCombined movingPeaksAdaBoost movingPeaksSVM movingPeaksCombined $(OUTPUT) *.o 
