// Derived concrete class for global background fit to MS stars.
// Created by Enrico Corsaro @ IvS - 21 January 2014
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "PuntoBackgroundModel.h"
// Implementations contained in "PuntoBackgroundModel.cpp"


#ifndef PUNTOBACKGROUNDMODEL_H
#define PUNTOBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class PuntoBackgroundModel : public BackgroundModel
{
    public:
    
        PuntoBackgroundModel(const RefArrayXd covariates);
        ~PuntoBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
