// Derived concrete class for global background fit to photometric observations.
// Created by Enrico Corsaro @ CEA - January 2016
// e-mail: emncorsaro@gmail.com
// Header file "FullBackgroundModel.h"
// Implementations contained in "FullBackgroundModel.cpp"


#ifndef FULLBACKGROUNDMODEL_H
#define FULLBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class FullBackgroundModel : public BackgroundModel
{
    public:
    
        FullBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFileName);
        ~FullBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif