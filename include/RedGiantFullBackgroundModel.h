// Derived concrete class for global background fit to RG stars.
// Created by Enrico Corsaro @ IvS - 13 February 2015
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "RedGiantFullBackgroundModel.h"
// Implementations contained in "RedGiantFullBackgroundModel.cpp"


#ifndef REDGIANTFULLBACKGROUNDMODEL_H
#define REDGIANTFULLBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class RedGiantFullBackgroundModel : public BackgroundModel
{
    public:
    
        RedGiantFullBackgroundModel(const RefArrayXd covariates);
        ~RedGiantFullBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
