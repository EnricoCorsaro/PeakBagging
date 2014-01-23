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
    // Create response function modulating the sampling rate of input Kepler LC data

    double NyquistFrequency = 283.2116656017908;    // muHz
    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square(); 
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
    double amplitudeHarvey1 = configuringParameters(1);
    double amplitudeHarvey2 = configuringParameters(2);
    double amplitudeHarvey3 = configuringParameters(3);
    double frequencyHarvey1 = configuringParameters(4);
    double frequencyHarvey2 = configuringParameters(5);
    double frequencyHarvey3 = configuringParameters(6);
    double heightOscillation = configuringParameters(7);
    double nuMax = configuringParameters(8);
    double sigma = configuringParameters(9);


    // Compute Harvey components and add them to the predictions

    predictions = 2*Functions::PI*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));
    predictions += 2*Functions::PI*amplitudeHarvey2*amplitudeHarvey2/(frequencyHarvey2*(1.0 + (covariates/frequencyHarvey2).pow(4)));
    predictions += 2*Functions::PI*amplitudeHarvey3*amplitudeHarvey3/(frequencyHarvey3*(1.0 + (covariates/frequencyHarvey3).pow(4)));


    // Compute Gaussian envelope for RGB Oscillations and add it to the predictions

    predictions += heightOscillation * exp(-1.0*(nuMax - covariates)*(nuMax - covariates)/(2.0 * sigma * sigma));
    
    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;        


    // Add flat noise level component

    predictions += flatNoiseLevel;

    
}











