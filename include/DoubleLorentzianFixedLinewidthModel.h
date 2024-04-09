// Derived class for fitting a double lorentzian profile with a fixed linewidth
// Created by Enrico Corsaro @ INAF-OACT - March 2019
// e-mail: emncorsaro@gmail.com
// Header file "DoubleLorentzianFixedLinewidthModel.h"
// Implementations contained in "DoubleLorentzianFixedLinewidthModel.cpp"


#ifndef DOUBLELORENTZIANFIXEDLINEWIDTHMODEL_H
#define DOUBLELORENTZIANFIXEDLINEWIDTHMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class DoubleLorentzianFixedLinewidthModel : public Model
{
    public:
    
        DoubleLorentzianFixedLinewidthModel(RefArrayXd const covariates, const double linewidth, BackgroundModel &backgroundModel, 
        	const string modeVisibilityFileName); 
        ~DoubleLorentzianFixedLinewidthModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};
        void readModeVisibilityFromFile(const string inputFileName);

    protected:
    

    private:

        double linewidth;                       // The fixed value of the linewidth for the profile used
        double quadrupoleToRadialHeightRatio;	// The mode visibility of the quadrupole mode
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
