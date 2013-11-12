// Derived class for building a mixture of lorentzian profiles
// Created by Enrico Corsaro @ IvS - 12 November 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "LorentzianMixtureModel.h"
// Implementations contained in "LorentzianMixtureModel.cpp"


#ifndef LORENTZIANMIXTUREMODEL_H
#define LORENTZIANMIXTUREMODEL_H

#include <iostream>
#include "Model.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LorentzianMixtureModel : public Model
{
    public:
    
        LorentzianMixtureModel(const RefArrayXd covariates, const vector<int> &NparametersPerType);
        ~LorentzianMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);


    protected:



    private:

        int NglobalParameters;                          // Total number of global parameters to be fitted 
                                                        // (e.g. nuMax, white noise, period spacing, etc.)
        int NprofileParameters;                         // Number of parameters determining the shape 
                                                        // of the mode profile (central frequency, height, linewidth, inclination angle, etc.)
        int Nmodes;                                     // Total number of modes to be fitted


}; 


#endif
