// Derived concrete class for global background fit to RG stars.
// Created by Enrico Corsaro @ IvS - 21 November 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "RedGiantBackgroundModel.h"
// Implementations contained in "RedGiantBackgroundModel.cpp"


#ifndef REDGIANTBACKGROUNDMODEL_H
#define REDGIANTBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class RedGiantBackgroundModel : public BackgroundModel
{
    public:
    
        RedGiantBackgroundModel(const RefArrayXd covariates);
        ~RedGiantBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
