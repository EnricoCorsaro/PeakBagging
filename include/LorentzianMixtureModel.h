// Derived class for building a mixture of lorentzian profiles
// Created by Enrico Corsaro @ IvS - 6 June 2014
// e-mail: emncorsaro@gmail.com
// Header file "LorentzianMixtureModel.h"
// Implementations contained in "LorentzianMixtureModel.cpp"


#ifndef LORENTZIANMIXTUREMODEL_H
#define LORENTZIANMIXTUREMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LorentzianMixtureModel : public Model
{
    public:
    
        LorentzianMixtureModel(RefArrayXd const covariates, const int Nresolved, BackgroundModel &backgroundModel); 
        ~LorentzianMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:


    private:

        int Nresolved;                          // Total number of p modes and resolved, partially resolved, or unresolved mixed modes to be fitted
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data

}; 


#endif
