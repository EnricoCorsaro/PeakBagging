// Derived class for fitting a single sinc^2 profile plus a background level,
// useful for performing the peak significance test.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestSincBackgroundModel.h"
// Implementations contained in "PeakTestSincBackgroundModel.cpp"


#ifndef PEAKTESTSINCBACKGROUNDMODEL_H
#define PEAKTESTSINCBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestSincBackgroundModel : public Model
{
    public:
    
        PeakTestSincBackgroundModel(RefArrayXd const covariates, const double frequencyResolution, BackgroundModel &backgroundModel); 
        ~PeakTestSincBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        double frequencyResolution;             // The frequency resolution of the input dataset
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
