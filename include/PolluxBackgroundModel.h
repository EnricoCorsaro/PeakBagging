// Derived concrete class for global background fit to RG stars.
// Created by Enrico Corsaro @ IvS - 21 November 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "PolluxBackgroundModel.h"
// Implementations contained in "PolluxBackgroundModel.cpp"


#ifndef POLLUXBACKGROUNDMODEL_H
#define POLLUXBACKGROUNDMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class PolluxBackgroundModel : public BackgroundModel
{
    public:
    
        PolluxBackgroundModel(const RefArrayXd covariates);
        ~PolluxBackgroundModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters){};
        virtual void predict(RefArrayXd predictions);


    protected:


    private:

}; 


#endif
