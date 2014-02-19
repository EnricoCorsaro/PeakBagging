#include "PuntoBackgroundModel.h"


// PuntoBackgroundModel::PuntoBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

PuntoBackgroundModel::PuntoBackgroundModel(const RefArrayXd covariates)
: BackgroundModel(covariates)
{
    // Create response function modulating the sampling rate of input Kepler SC data

    double NyquistFrequency = 8496.355743094671;    // muHz
    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square(); 
}










// PuntoBackgroundModel::PuntoBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

PuntoBackgroundModel::~PuntoBackgroundModel()
{

}










// PuntoBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a Background model based on a configuration
//      file that contains all the free parameters of the given model.
//      The model for Punto (KIC 9139163) consists in an exponential term
//      and a granulation component, overlaid to a flat noise level.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//
// OUTPUT:
//      void
//

void PuntoBackgroundModel::predict(RefArrayXd predictions)
{
    Nparameters = configuringParameters.size();


    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double heightPowerLaw = exp(configuringParameters(1));
    double exponentPowerLaw = configuringParameters(2);
    double amplitudeHarvey1 = configuringParameters(3);
    double timescaleHarvey1 = configuringParameters(4);
    double exponentHarvey1 = configuringParameters(5);


    // Compute decaying power law component
    
    predictions = heightPowerLaw * covariates.pow(exponentPowerLaw);


    // Compute Harvey-like components and add them to the predictions

    predictions += 4*amplitudeHarvey1*amplitudeHarvey1 * (timescaleHarvey1 / 1.e6) / 
                   (1.0 + (2*Functions::PI*covariates*timescaleHarvey1/1.e6).pow(exponentHarvey1));


    // Modulate the model by the response function (apodization)

    predictions *= responseFunction;


    // Add flat noise level component

    predictions += flatNoiseLevel;
}











