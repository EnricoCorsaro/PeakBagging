// Derived class for building a sliding pattern from the asymptotic relation for p modes
// Created by Enrico Corsaro @ OACT - August 2019
// e-mail: emncorsaro@gmail.com
// Header file "AsymptoticPatternModel.h"
// Implementations contained in "AsymptoticPatternModel.cpp"


#ifndef SLIDINGPATTERNMODEL_H
#define SLIDINGPATTERNMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class AsymptoticPatternModel : public Model
{
    public:
    
        AsymptoticPatternModel(RefArrayXd const covariates, const double linewidth, const string asymptoticParametersFileName, BackgroundModel &backgroundModel); 
        ~AsymptoticPatternModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        void readAsymptoticParametersFromFile(const string inputFileName);

    protected:


    private:

        double linewidth;                       // The linewidth of the radial modes
        int Norders;                            // Total number of radial orders to build the sliding pattern
        double DeltaNu;                         // The large frequency separation
        double deltaNu02;                       // The small frequency separation e02
        double deltaNu01;                       // The small frequency separation d01
        double dipoleToRadialHeightRatio;       // The ratio between the assumed height of a dipole region of modes and that of a radial mode
        double quadrupoleToRadialHeightRatio;   // The ratio between the assumed height of a quadrupole region of modes and that of a radial mode
        double dipoleToRadialLinewidthRatio;    // The ratio between the assumed linewidth of a dipole region of modes and that of a radial mode
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data

}; 


#endif
