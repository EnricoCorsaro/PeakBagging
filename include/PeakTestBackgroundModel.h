// Derived class for fitting a background level for the peak significance test.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestBackgroundModel.h"
// Implementations contained in "PeakTestBackgroundModel.cpp"


#ifndef PEAKTESTBACKGROUNDMODEL_H
#define PEAKTESTBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestBackgroundModel : public Model
{
    public:
    
        PeakTestBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
