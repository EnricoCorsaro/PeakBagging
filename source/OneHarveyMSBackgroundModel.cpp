#include "OneHarveyMSBackgroundModel.h"


// OneHarveyMSBackgroundModel::OneHarveyMSBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

OneHarveyMSBackgroundModel::OneHarveyMSBackgroundModel(const RefArrayXd covariates)
: BackgroundModel(covariates)
{
    // Create response function modulating the sampling rate of input Kepler SC data

    double NyquistFrequency = 8496.355743094671;    // muHz
    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square(); 
}










// OneHarveyMSBackgroundModel::OneHarveyMSBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

OneHarveyMSBackgroundModel::~OneHarveyMSBackgroundModel()
{

}










// OneHarveyMSBackgroundModel::predict()
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

void OneHarveyMSBackgroundModel::predict(RefArrayXd predictions)
{
    Nparameters = configuringParameters.size();


    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double heightPowerLaw = configuringParameters(1);
    double exponentPowerLaw = configuringParameters(2);
    double amplitudeHarvey1 = configuringParameters(3);
    double timescaleHarvey1 = configuringParameters(4);
    double exponentHarvey1 = configuringParameters(5);
    double heightOscillation = configuringParameters(6);
    double nuMax = configuringParameters(7);
    double sigma = configuringParameters(8);


    // Compute decaying power law component
    
    predictions = heightPowerLaw * covariates.pow(exponentPowerLaw);


    // Compute Harvey-like components and add them to the predictions

    predictions += 4*amplitudeHarvey1*amplitudeHarvey1 * (timescaleHarvey1 / 1.e6) / 
                   (1.0 + (2*Functions::PI*covariates*timescaleHarvey1/1.e6).pow(exponentHarvey1));


    // Compute Gaussian envelope for MS Oscillations and add it to the predictions

    predictions += heightOscillation * exp(-1.0*(nuMax - covariates)*(nuMax - covariates)/(2.0 * sigma * sigma));


    // Modulate the model by the response function (apodization)

    predictions *= responseFunction;


    // Add flat noise level component

    predictions += flatNoiseLevel;
}











