// Derived class for building a mixture of lorentzian profiles
// Created by Enrico Corsaro @ IvS - 6 June 2014
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MixedModesMixtureModel.h"
// Implementations contained in "MixedModesMixtureModel.cpp"


#ifndef MIXEDMODESMIXTUREMODEL_H
#define MIXEDMODESMIXTUREMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MixedModesMixtureModel : public Model
{
    public:
    
        MixedModesMixtureModel(const RefArrayXd covariates, const int Nresolved, const int Nunresolved,
                               const double frequencyResolution, BackgroundModel &backgroundModel);
        ~MixedModesMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);


    protected:


    private:

        int Nresolved;                                  // Total number of p modes and resolved (or partially resolved) mixed modes to be fitted
        int Nunresolved;                                // Total number of unresolved mixed modes to be fitted
        double frequencyResolution;                     // The frequency bin size used to compute the sinc-square profile
        ArrayXd backgroundPrediction;                   // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;                       // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundParameters;                   // An array containing the configuring parameters for the background model adopted

}; 


#endif
