// Derived concrete class for global background fit to RG stars.
// Created by Enrico Corsaro @ IvS - 21 November 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "CanBackgroundModel.h"
// Implementations contained in "CanBackgroundModel.cpp"


#ifndef CANBACKGROUNDMODEL_H
#define CANBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class CanBackgroundModel : public BackgroundModel
{
    public:
    
        CanBackgroundModel(const RefArrayXd covariates);
        ~CanBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
