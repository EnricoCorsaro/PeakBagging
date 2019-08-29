// Derived class for fitting a single Lorentzian profile with rotationally
// split components (intended as a dipole splitting).
// Created by Enrico Corsaro @ INAF-OACT - July 2019
// e-mail: emncorsaro@gmail.com
// Header file "PeakTestDipoleRotationModel.h"
// Implementations contained in "PeakTestDipoleRotationModel.cpp"


#ifndef PEAKTESTDIPOLEROTATIONMODEL_H
#define PEAKTESTDIPOLEROTATIONMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PeakTestDipoleRotationModel : public Model
{
    public:
    
        PeakTestDipoleRotationModel(RefArrayXd const covariates, BackgroundModel &backgroundModel); 
        ~PeakTestDipoleRotationModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
