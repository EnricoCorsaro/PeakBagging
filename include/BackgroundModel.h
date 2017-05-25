// Derived abstract class for global background fit to RG stars.
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
    
        BackgroundModel(const RefArrayXd covariates);
        ~BackgroundModel();

        void readNyquistFrequencyFromFile(const string inputFileName);
        double getNyquistFrequency();
        ArrayXd getConfiguringParameters();
        ArrayXd getResponseFunction();

        virtual void readConfiguringParametersFromFile(const string inputFileName) = 0;
        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters) = 0;
        virtual void predict(RefArrayXd predictions) = 0;
    


    protected:

        double NyquistFrequency;
        ArrayXd configuringParameters;
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data


    private:

}; 


#endif
