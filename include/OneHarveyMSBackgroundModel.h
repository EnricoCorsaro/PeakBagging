// Derived concrete class for global background fit to MS stars.
// Created by Enrico Corsaro @ IvS - 21 January 2014
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "OneHarveyMSBackgroundModel.h"
// Implementations contained in "OneHarveyMSBackgroundModel.cpp"


#ifndef ONEHARVEYMSBACKGROUNDMODEL_H
#define ONEHARVEYMSBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class OneHarveyMSBackgroundModel : public BackgroundModel
{
    public:
    
        OneHarveyMSBackgroundModel(const RefArrayXd covariates);
        ~OneHarveyMSBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
