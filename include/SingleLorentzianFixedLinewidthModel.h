// Derived class for fitting a single lorentzian profile with a fixed linewidth
// Created by Enrico Corsaro @ INAF-OACT - 8 January 2019
// e-mail: emncorsaro@gmail.com
// Header file "SingleLorentzianFixedLinewidthModel.h"
// Implementations contained in "SingleLorentzianFixedLinewidthModel.cpp"


#ifndef SINGLELORENTZIANFIXEDLINEWIDTHMODEL_H
#define SINGLELORENTZIANFIXEDLINEWIDTHMODEL_H

#include <iostream>
#include "BackgroundModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class SingleLorentzianFixedLinewidthModel : public Model
{
    public:
    
        SingleLorentzianFixedLinewidthModel(RefArrayXd const covariates, const double linewidth, BackgroundModel &backgroundModel); 
        ~SingleLorentzianFixedLinewidthModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:
    

    private:

        double linewidth;                       // The fixed value of the linewidth for the profile used
        ArrayXd backgroundPrediction;           // An array containing the prediction for the background in all the range of the covariates
        ArrayXd responseFunction;               // An array containing the apodization response function for the signal of the input data
}; 


#endif
