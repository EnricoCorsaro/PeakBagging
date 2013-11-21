#include "CanBackgroundModel.h"


// CanBackgroundModel::CanBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

CanBackgroundModel::CanBackgroundModel(const RefArrayXd covariates)
: BackgroundModel(covariates)
{
}










// CanBackgroundModel::CanBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

CanBackgroundModel::~CanBackgroundModel()
{

}










// CanBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a Background model based on a configuration
//      file that contains all the free parameters of the given model.
//      This is an overloaded function, and its implementation is only given
//      by the concrete derived class that implements the desired model for
//      modeling the background.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//
// OUTPUT:
//      void
//

void CanBackgroundModel::predict(RefArrayXd predictions)
{
    Nparameters = configuringParameters.size();


    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitude1 = configuringParameters(1);
    double amplitude2 = configuringParameters(2);
    double amplitude3 = configuringParameters(3);
    double noiseFrequency1 = configuringParameters(4);
    double noiseFrequency2 = configuringParameters(5);
    double noiseFrequency3 = configuringParameters(6);


    // Compute power law components and add them to the predictions

    predictions = 2*Functions::PI*amplitude1*amplitude1/(noiseFrequency1*(1.0 + pow(covariates/noiseFrequency1,4)));
    predictions += 2*Functions::PI*amplitude2*amplitude2/(noiseFrequency2*(1.0 + pow(covariates/noiseFrequency2,4)));
    predictions += 2*Functions::PI*amplitude3*amplitude3/(noiseFrequency3*(1.0 + pow(covariates/noiseFrequency3,4)));


    // Add flat noise level component

    predictions += flatNoiseLevel;           
}











