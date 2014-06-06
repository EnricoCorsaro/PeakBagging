#include "RedGiantBackgroundModel.h"


// RedGiantBackgroundModel::RedGiantBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

RedGiantBackgroundModel::RedGiantBackgroundModel(const RefArrayXd covariates)
: BackgroundModel(covariates)
{
    // Create response function modulating the sampling rate of input Kepler LC data

    double NyquistFrequency = 283.2116656017908;    // muHz
    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square(); 
}










// RedGiantBackgroundModel::RedGiantBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

RedGiantBackgroundModel::~RedGiantBackgroundModel()
{

}










// RedGiantBackgroundModel::predict()
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

void RedGiantBackgroundModel::predict(RefArrayXd predictions)
{
    Nparameters = configuringParameters.size();


    // Initialize global parameters

    double flatNoiseLevel = modelParameters(0);
    double amplitudeHarvey1 = modelParameters(1);
    double frequencyHarvey1 = modelParameters(2);
    double amplitudeHarvey2 = modelParameters(3);
    double frequencyHarvey2 = modelParameters(4);
    double amplitudeHarvey3 = modelParameters(5);
    double frequencyHarvey3 = modelParameters(6);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));
    predictions += zeta*amplitudeHarvey2*amplitudeHarvey2/(frequencyHarvey2*(1.0 + (covariates/frequencyHarvey2).pow(4)));
    predictions += zeta*amplitudeHarvey3*amplitudeHarvey3/(frequencyHarvey3*(1.0 + (covariates/frequencyHarvey3).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise level component

    predictions += flatNoiseLevel;
    
}











