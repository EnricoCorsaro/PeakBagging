// Derived class for building a mixture of lorentzian profiles
// Created by Enrico Corsaro @ IvS - 2 March 2014
// e-mail: emncorsaro@gmail.com
// Header file "LorentzianMixtureRotationModel.h"
// Implementations contained in "LorentzianMixtureRotationModel.cpp"


#ifndef LORENTZIANMIXTUREROTATIONMODEL_H
#define LORENTZIANMIXTUREROTATIONMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LorentzianMixtureRotationModel : public Model
{
    public:
    
        LorentzianMixtureRotationModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, BackgroundModel &backgroundModel);
        ~LorentzianMixtureRotationModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);


    protected:


    private:

        int NprofileParameters;                         // Number of parameters determining the shape 
                                                        // of the mode profile (central frequency, height/amplitude, linewidth)
        int Nmodes;                                     // Total number of modes to be fitted
        int startingModeNumber;                         // The reference mode number for which to include the split components (between 0 and 2).
                                                        // If set to zero, the rotation split components will refer to the first mode to the 
                                                        // left of the series of peaks in one radial order.
        ArrayXd backgroundPrediction;                   // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;                       // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundParameters;                   // An array containing the configuring parameters for the background model adopted
}; 


#endif
