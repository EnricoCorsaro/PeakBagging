#include "FullBackgroundModel.h"


// FullBackgroundModel::FullBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      inputNyquestFrequencyFileName:      the string containing the file name of the input ASCII file with the
//                                          value of the Nyquist frequency to be adopted in the response function
//

FullBackgroundModel::FullBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFileName)
: BackgroundModel(covariates)
{
    // Create response function modulating the sampling rate of input Kepler LC data

    // NyquistFrequency = 8496.355743094671     muHz     // Kepler SC
    // NyquistFrequency = 283.2116656017908     muHz     // Kepler LC

    readNyquistFrequencyFromFile(inputNyquistFrequencyFileName);

    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square();
}










// FullBackgroundModel::FullBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

FullBackgroundModel::~FullBackgroundModel()
{

}










// FullBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a Background model based on a configuration
//      file that contains all the free parameters of the given model.
//      This is an overloaded function whose implementation is that of a 
//      background model for a red giant star, according to the findings
//      by Kallinger et al. 2014. This includes a colored noise for low-numax stars.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//
// OUTPUT:
//      void
//

void FullBackgroundModel::predict(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeNoise = configuringParameters(1);
    double frequencyNoise = configuringParameters(2);
    double amplitudeHarvey1 = configuringParameters(3);
    double frequencyHarvey1 = configuringParameters(4);
    double amplitudeHarvey2 = configuringParameters(5);
    double frequencyHarvey2 = configuringParameters(6);
    double amplitudeHarvey3 = configuringParameters(7);
    double frequencyHarvey3 = configuringParameters(8);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));
    predictions += zeta*amplitudeHarvey2*amplitudeHarvey2/(frequencyHarvey2*(1.0 + (covariates/frequencyHarvey2).pow(4)));
    predictions += zeta*amplitudeHarvey3*amplitudeHarvey3/(frequencyHarvey3*(1.0 + (covariates/frequencyHarvey3).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise and colored noise components

    predictions += flatNoiseLevel;
    predictions += 2.0*Functions::PI*amplitudeNoise*amplitudeNoise/(frequencyNoise*(1.0 + (covariates/frequencyNoise).pow(2)));
    
}
