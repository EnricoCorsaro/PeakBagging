// Derived class for fitting a single lorentzian profile with rotationally
// split components.
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestLorentzianRotationModel.h"
// Implementations contained in "PeakTestLorentzianRotationModel.cpp"


#ifndef PEAKTESTLORENTZIANROTATIONMODEL_H
#define PEAKTESTLORENTZIANROTATIONMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestLorentzianRotationModel : public Model
{
    public:
    
        PeakTestLorentzianRotationModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestLorentzianRotationModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
