// Derived class for computing the global background level in stellar power spectra.
// Created by Enrico Corsaro @ IvS - 21 November 2013
// e-mail: emncorsaro@gmail.com
// Header file "BackgroundModel.h"
// Implementations contained in "BackgroundModel.cpp"


#ifndef BACKGROUNDMODEL_H
#define BACKGROUNDMODEL_H

#include <iostream>
#include "Model.h"
#include "Functions.h"
#include "File.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class BackgroundModel : public Model
{
    public:
    
        BackgroundModel(const RefArrayXd covariates, const string backgroundModelName);
        ~BackgroundModel();

        void readNyquistFrequencyFromFile(const string inputNyquistFrequencyFileName);
        void readConfiguringParametersFromFile(const string inputFileName);
        void writeBackgroundPredictionToFile(const string outputFileName);

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};
        
        void predict(RefArrayXd predictions);
        void predictThreeHarveyColor(RefArrayXd predictions);
        void predictThreeHarvey(RefArrayXd predictions);
        void predictTwoHarveyColor(RefArrayXd predictions);
        void predictTwoHarvey(RefArrayXd predictions);
        void predictOneHarveyColor(RefArrayXd predictions);
        void predictOneHarvey(RefArrayXd predictions);
        void predictOriginal(RefArrayXd predictions);
        void predictFlat(RefArrayXd predictions);

        string getBackgroundModelName();
        double getNyquistFrequency();
        ArrayXd getConfiguringParameters();
        ArrayXd getResponseFunction();

    protected:

        string backgroundModelName;
        double NyquistFrequency;
        ArrayXd configuringParameters;
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundPrediction;           // An array containing the background prediction of the input dataset


    private:

}; 


#endif
