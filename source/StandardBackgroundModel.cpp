#include "StandardBackgroundModel.h"


// StandardBackgroundModel::StandardBackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent
//      inputNyquestFrequencyFileName:      the string containing the file name of the input ASCII file with the
//                                          value of the Nyquist frequency to be adopted in the response function variable.
//

StandardBackgroundModel::StandardBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFileName)
: BackgroundModel(covariates)
{
    // NyquistFrequency = 8496.355743094671     muHz     // Kepler SC
    // NyquistFrequency = 283.2116656017908     muHz     // Kepler LC

    readNyquistFrequencyFromFile(inputNyquistFrequencyFileName);

    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square();
}










// StandardBackgroundModel::StandardBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

StandardBackgroundModel::~StandardBackgroundModel()
{

}










// StandardBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a Background model based on a configuration
//      file that contains all the free parameters of the given model.
//      This is an overloaded function whose implementation is that of a 
//      background model for a red giant star, according to the findings
//      by Kallinger et al. 2014.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//
// OUTPUT:
//      void
//

void StandardBackgroundModel::predict(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeHarvey1 = configuringParameters(1);
    double frequencyHarvey1 = configuringParameters(2);
    double amplitudeHarvey2 = configuringParameters(3);
    double frequencyHarvey2 = configuringParameters(4);
    double amplitudeHarvey3 = configuringParameters(5);
    double frequencyHarvey3 = configuringParameters(6);


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
