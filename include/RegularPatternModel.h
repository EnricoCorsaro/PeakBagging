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
    
        RegularPatternModel(const RefArrayXd covariates, const int NradialOrders);
        ~RegularPatternModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        int getNparameters();
       

    protected:
        
        double computeSinglePmodeFrequency(int relativeRadialOrder, double frequencyOfCentralRadialMode, 
                                           double DeltaNu, double deltaNu02, double epsilon, int angularDegree);
        void computePmodeFrequencies(vector<int> relativeRadialOrders, double frequencyOfCentralRadialMode, 
                                            double DeltaNu, double deltaNu02);
        // void computeMixedDipoleModeFrequencies(double frequencyOfCentralRadialMode, double DeltaNu, double DeltaPi1,
        //                                        vector<int> relativeRadialOrders, double couplingFactor){};



    private:

        int NradialOrders;                              // The number of radial orders of the model
        int Nparameters;                                // The number of free parameters of the model
        int radialOrderOfCentralRadialMode;             // The radial order of the central radial mode
        vector<double> frequenciesOfPmodes;
        vector<double> angularDegreesOfPmodes;
        vector<double> radialOrdersOfPmodes;
   

}; // END class RegularPatternModel


#endif
