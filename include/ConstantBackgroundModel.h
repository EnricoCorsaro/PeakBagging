// Derived concrete class for constant background model.
// Created by Enrico Corsaro @ OACT - October 2017
// e-mail: emncorsaro@gmail.com
// Header file "ConstantBackgroundModel.h"
// Implementations contained in "ConstantBackgroundModel.cpp"


#ifndef CONSTANTBACKGROUNDMODEL_H
#define CONSTANTBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class ConstantBackgroundModel : public BackgroundModel
{
    public:
    
        ConstantBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFromFile);
        ~ConstantBackgroundModel();

        virtual void readConfiguringParametersFromFile(const string inputFileName);
        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
