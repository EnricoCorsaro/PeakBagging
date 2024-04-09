// Derived class for building a sliding pattern from the asymptotic relation for p modes
// Created by Enrico Corsaro @ OACT - August 2019
// e-mail: emncorsaro@gmail.com
// Header file "SlidingPatternModel.h"
// Implementations contained in "SlidingPatternModel.cpp"


#ifndef SLIDINGPATTERNMODEL_H
#define SLIDINGPATTERNMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class SlidingPatternModel : public Model
{
    public:
    
        SlidingPatternModel(RefArrayXd const covariates, const double linewidth, const ArrayXd slidingPatternParameters, 
            const string asymptoticParametersFileName, const string gaussianEnvelopeParametersFileName, 
            BackgroundModel &backgroundModel); 
        ~SlidingPatternModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};
        void readAsymptoticParametersFromFile(const string inputFileName);
        void readGaussianEnvelopeParametersFromFile(const string inputFileName);


    protected:


    private:

        bool includeDipoleModes;                // Boolean variables to either include or not the individual modes in the sliding pattern
        bool includeQuadrupoleModes;
        bool includeOctupoleModes;

        double linewidth;                       // The linewidth of the radial modes
        double modeVisibilityComponent;         // The visibility of the m-component in units of height in the PSD 
        ArrayXd slidingPatternParameters;       // The parameters of the sliding pattern

        int Norders;                            // Total number of radial orders to build the sliding pattern
        double dipoleToRadialHeightRatio;       // The ratio between the assumed height of a dipole region of modes and that of a radial mode
        double quadrupoleToRadialHeightRatio;   // The ratio between the assumed height of a quadrupole region of modes and that of a radial mode
        double octupoleToRadialHeightRatio;     // The ratio between the assumed height of an octupole region of modes and that of a radial mode
        double dipoleToRadialLinewidthRatio;    // The ratio between the assumed linewidth of a dipole region of modes and that of a radial mode
        double nuMax;                           // The frequency of maximum oscillation power
        double standardDeviation;               // The width of the Gaussian oscillation envelope given as standard deviation
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data

}; 


#endif
