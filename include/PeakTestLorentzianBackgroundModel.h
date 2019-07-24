// Derived class for fitting a single lorentzian profile plus a background level,
// useful for performing the peak significance test.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestLorentzianBackgroundModel.h"
// Implementations contained in "PeakTestLorentzianBackgroundModel.cpp"


#ifndef PEAKTESTLORENTZIANBACKGROUNDMODEL_H
#define PEAKTESTLORENTZIANBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestLorentzianBackgroundModel : public Model
{
    public:
    
        PeakTestLorentzianBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestLorentzianBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
