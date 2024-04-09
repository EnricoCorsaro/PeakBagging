// Derived class for building a mixture of lorentzian and sinc profiles with rotational effect included
// Created by Enrico Corsaro @ CEA - November 2015
// e-mail: emncorsaro@gmail.com
// Header file "LorentzianRotationMixtureModel.h"
// Implementations contained in "LorentzianRotationMixtureModel.cpp"


#ifndef LORENTZIANROTATIONMIXTUREMODEL_H
#define LORENTZIANROTATIONMIXTUREMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LorentzianRotationMixtureModel : public Model
{
    public:
    
        LorentzianRotationMixtureModel(RefArrayXd const covariates, const int Nresolved, string angularDegreesFileName, BackgroundModel &backgroundModel); 
        ~LorentzianRotationMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};
        ArrayXd computeModeVisibility(const int angularDegree, const double cosi);


    protected:

        void readAngularDegreesFromFile(const string inputFileName);


    private:

        int Nresolved;                          // Total number modes considered in the fit.
        string angularDegreesFileName;          // a string containing the file name of the ASCII file with the list of angular degrees for all
                                                // the modes considered in the fit.
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
        ArrayXd angularDegreeArray;             // An array of dimension Nresolved+Nunresolved containing the angular degree of the mode
}; 


#endif
