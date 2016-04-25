// Derived class for building mono lorentzian profile models
// Created by Enrico Corsaro @ IvS - 29 March 2014
// e-mail: emncorsaro@gmail.com
// Header file "IslandsModel.h"
// Implementations contained in "IslandsModel.cpp"


#ifndef ISLANDSMODEL_H
#define ISLANDSMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class IslandsModel : public Model
{
    public:
    
        IslandsModel(const RefArrayXd covariates, BackgroundModel &backgroundModel);
        ~IslandsModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);


    protected:
        

    private:

        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
        ArrayXd backgroundParameters;           // An array containing the configuring parameters for the background model adopted

};


#endif
