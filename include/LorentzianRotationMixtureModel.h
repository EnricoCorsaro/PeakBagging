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
    
        LorentzianRotationMixtureModel(RefArrayXd const covariates, const int NmodesWithRotation, const int NmodesWithoutRotation, 
                                       const double cosi, string angularDegreesFileName, BackgroundModel &backgroundModel); 
        ~LorentzianRotationMixtureModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        void computeModeVisibility();
        void createModeAmplitudeVisibility();


    protected:

        void readAngularDegreesFromFile(const string inputFileName);


    private:

        int NmodesWithRotation;                 // Total number of p modes in which the rotational effect is fitted
        int NmodesWithoutRotation;              // Total number of p modes in which the rotational effect is not fitted
        double cosi;                            // Solar rotation axis inclination angle
        string angularDegreesFileName;          // a string containing the file name of the ASCII file with the list of angular degrees for the modes
                                                // where rotation has to be fitted
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundParameters;           // An array containing the configuring parameters for the background model adopted
        ArrayXd angularDegreeArray;             // An array of dimension Nresolved+Nunresolved containing the angular degree of the mode
                                                // The elements must be sorted in the same order as the modes are defined in the prior list
        ArrayXd modeAmplitudeVisibilityArray;   // The visibility rescaled for the amplitude for each angular degree and azimuthal component 
                                                // following the input order of the modes with rotation
        ArrayXd modeVisibility1;                // Mode visibility for l=1 modes
        ArrayXd modeVisibility2;                // Mode visibility for l=2 modes
        ArrayXd modeVisibility3;                // Mode visibility for l=3 modes
        ArrayXd modeVisibility4;                // Mode visibility for l=4 modes
}; 


#endif
