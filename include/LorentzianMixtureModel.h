// Derived class for building a mixture of lorentzian profiles
// Created by Enrico Corsaro @ IvS - 12 November 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "LorentzianMixtureModel.h"
// Implementations contained in "LorentzianMixtureModel.cpp"


#ifndef LORENTZIANMIXTUREMODEL_H
#define LORENTZIANMIXTUREMODEL_H

#include <iostream>
#include "Model.h"
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LorentzianMixtureModel : public Model
{
    public:
    
        LorentzianMixtureModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, BackgroundModel &backgroundModel);
        ~LorentzianMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);


    protected:


    private:

        int NprofileParameters;                         // Number of parameters determining the shape 
                                                        // of the mode profile (central frequency, height, linewidth, inclination angle, etc.)
        int Nmodes;                                     // Total number of modes to be fitted
        ArrayXd backgroundPrediction;                   // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;                       // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundParameters;                   // An array containing the configuring parameters for the background model adopted

}; 


#endif
