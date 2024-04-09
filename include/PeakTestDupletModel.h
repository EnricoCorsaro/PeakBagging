// Derived class for fitting a duplet (two lorentzian profiles) with fixed background.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestDupletModel.h"
// Implementations contained in "PeakTestDupletModel.cpp"

#ifndef PEAKTESTDUPLETMODEL_H
#define PEAKTESTDUPLETMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestDupletModel : public Model
{
    public:
    
        PeakTestDupletModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestDupletModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
