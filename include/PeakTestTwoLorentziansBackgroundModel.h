// Derived class for fitting two lorentzian profiles.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestTwoLorentziansBackgroundModel.h"
// Implementations contained in "PeakTestTwoLorentziansBackgroundModel.cpp"


#ifndef PEAKTESTTWOLORENTZIANSBACKGROUNDMODEL_H
#define PEAKTESTTWOLORENTZIANSBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestTwoLorentziansBackgroundModel : public Model
{
    public:
    
        PeakTestTwoLorentziansBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestTwoLorentziansBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
