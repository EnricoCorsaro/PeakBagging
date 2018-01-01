// Derived concrete class for standard background model (no colored-noise component) for oscillating stars with photometry.
// Created by Enrico Corsaro @ CEA - January 2016
// e-mail: emncorsaro@gmail.com
// Header file "StandardBackgroundModel.h"
// Implementations contained in "StandardBackgroundModel.cpp"


#ifndef STANDARDBACKGROUNDMODEL_H
#define STANDARDBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class StandardBackgroundModel : public BackgroundModel
{
    public:
    
        StandardBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFromFile);
        ~StandardBackgroundModel();

        virtual void readConfiguringParametersFromFile(const string inputFileName);
        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
