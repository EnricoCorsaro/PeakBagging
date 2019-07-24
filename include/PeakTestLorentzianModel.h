// Derived class for fitting a single lorentzian profile.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestLorentzianModel.h"
// Implementations contained in "PeakTestLorentzianModel.cpp"


#ifndef PEAKTESTLORENTZIANMODEL_H
#define PEAKTESTLORENTZIANMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestLorentzianModel : public Model
{
    public:
    
        PeakTestLorentzianModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestLorentzianModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
