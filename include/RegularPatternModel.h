// Derived class for building mono lorentzian profile models
// Created by Enrico Corsaro @ IvS - 13 June 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "RegularPatternModel.h"
// Implementations contained in "RegularPatternModel.cpp"


#ifndef REGULARPATTERNMODEL_H
#define REGULARPATTERNMODEL_H

#include <iostream>
#include "Model.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class RegularPatternModel : public Model
{
    public:
    
        RegularPatternModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, double nuMax);
        ~RegularPatternModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        double predictDeltaNuFromNuMax();
        double predictDeltaNu01FromDeltaNu(double DeltaNu);
        double predictDeltaNu02FromDeltaNu(double DeltaNu);
        double predictDeltaNu03FromDeltaNu(double DeltaNu);


    protected:
        
        double computeSinglePmodeFrequency(int relativeRadialOrder, double frequencyOfCentralRadialMode, 
                                           double DeltaNu, double deltaNu02, double deltaNu01, int angularDegree);
        // void computeMixedDipoleModeFrequencies(double frequencyOfCentralRadialMode, double DeltaNu, double DeltaPi1,
        //                                        double couplingFactor){};



    private:

        int NglobalParameters;                          // Total number of global parameters to be fitted 
                                                        // (e.g. nuMax, white noise, period spacing, etc.)
        int NprofileParameters;                         // Number of parameters determining the shape 
                                                        // of the mode profile (height, linewidth, inclination angle, etc.)
                                                        // N.B the mode frequency may not be considered as a profile parameter
        int NradialOrders;                              // Number of radial orders expected in the PSD
        int NangularDegrees;                            // Maximum number of angular degrees to be fitted
        int NpressureModes;                             // Total number of pressure modes to be fitted 
                                                        // N.B this number can differ from the product NangularDegrees*NradialOrders
        int NmixedModes;                                // Total number of mixed modes to be fitted
        double nuMax;                                   // The frequency of maximum power excess
        vector<int> relativeRadialOrders;               // Array of integers containing the radial order numbers with respect to
                                                        // the central one


}; // END class RegularPatternModel


#endif
