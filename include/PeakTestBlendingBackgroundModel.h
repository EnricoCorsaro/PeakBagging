// Derived class for fitting two lorentzian profiles in a blended form.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestBlendingBackgroundModel.h"
// Implementations contained in "PeakTestBlendingBackgroundModel.cpp"


#ifndef PEAKTESTBLENDINGBACKGROUNDMODEL_H
#define PEAKTESTBLENDINGBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestBlendingBackgroundModel : public Model
{
    public:
    
        PeakTestBlendingBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestBlendingBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
