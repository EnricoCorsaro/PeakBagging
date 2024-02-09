#include "BackgroundModel.h"


// BackgroundModel::BackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

BackgroundModel::BackgroundModel(const RefArrayXd covariates, const string backgroundModelName)
: Model(covariates),
  backgroundModelName(backgroundModelName)
{
}










// BackgroundModel::~BackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

BackgroundModel::~BackgroundModel()
{

}











// BackgroundModel::readNyquistFrequencyFromFile()
//
// PURPOSE:
//      Reads the Nyquist frequency of the dataset from an input ASCII file.
//
// INPUT:
//      inputNyquistFrequencyFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void BackgroundModel::readNyquistFrequencyFromFile(const string inputNyquistFrequencyFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputNyquistFrequencyFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXd inputData;
    inputData.resize(Nrows);

    inputData = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    NyquistFrequency = inputData(0);

    inputFile.close();
}











// FullBackgroundModel::readConfiguringParametersFromFile()
//
// PURPOSE:
//      Reads the configuring parameters of the model from an input file
//      specified by its full path as a string. The values are then
//      store into a vector which is a protected data member and that can
//      be therefore accessed from a derived class that implements the model.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void BackgroundModel::readConfiguringParametersFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Ncols == 0.0)
    {
        cerr << "Wrong number of input configuring parameters for the background model." << endl;
        exit(EXIT_FAILURE);
    }

    configuringParameters.resize(Nrows);
    configuringParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);

    inputFile.close();
}












// FullBackgroundModel::writeBackgroundPredictionToFile()
//
// PURPOSE:
//      Stores the background level as a function of PSD frequency into a two-column
//      ASCII file.
//
// INPUT:
//      outputFileName:      a string specifying the full path (filename included) of the out file to create.
//
// OUTPUT:
//      void
//

void BackgroundModel::writeBackgroundPredictionToFile(const string outputFileName)
{
    backgroundPrediction.resize(covariates.size());
    BackgroundModel::predict(backgroundPrediction);
    
    ofstream outputFile;
    File::openOutputFile(outputFile, outputFileName);

    outputFile << "# Background level as a function of frequency and corrected by apodization." << endl;
    outputFile << "# log(Likelihood)" << endl;
    outputFile << scientific << setprecision(9);
    
    ArrayXXd backgroundLevel(covariates.size(),2);
    backgroundLevel.col(0) = covariates;
    backgroundLevel.col(1) = backgroundPrediction;

    File::arrayXXdToFile(outputFile, backgroundLevel);
    outputFile.close();
}













// BackgroundModel::predict()
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

void BackgroundModel::predict(RefArrayXd predictions)
{
    // Create response function modulating the sampling rate of input data

    // NyquistFrequency = 8496.355743094671     muHz     // Kepler SC
    // NyquistFrequency = 283.2116656017908     muHz     // Kepler LC

    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square();
    bool backgroundNameMatched = false;
    
    // Long-trend, meso-granulation, and granulation component included, with colored noise
    if (backgroundModelName == "ThreeHarveyColor")
    {
        predictThreeHarveyColor(predictions); 
        backgroundNameMatched = true;
    }
    
    // Long-trend, meso-granulation, and granulation component included, but no colored noise
    if (backgroundModelName == "ThreeHarvey")
    {
        predictThreeHarvey(predictions); 
        backgroundNameMatched = true;
    }
    
    // Meso-granulation and granulation components included, with colored noise
    if (backgroundModelName == "TwoHarveyColor")
    {
        predictTwoHarveyColor(predictions); 
        backgroundNameMatched = true;
    }
    
    // Meso-granulation and granulation components included, but no colored noise
    if (backgroundModelName == "TwoHarvey")
    {
        predictTwoHarvey(predictions);  
        backgroundNameMatched = true;
    }

    // Only meso-granulation component included, with colored noise
    if (backgroundModelName == "OneHarveyColor")
    {    
        predictOneHarveyColor(predictions); 
        backgroundNameMatched = true;
    }
    
    // Only meso-granulation component included, but no colored noise
    if (backgroundModelName == "OneHarvey")
    {    
        predictOneHarvey(predictions); 
        backgroundNameMatched = true;
    }

    // Only meso-granulation component included, but no colored noise
    if (backgroundModelName == "Original")
    {    
        predictOriginal(predictions); 
        backgroundNameMatched = true;
    }

    // Only Gaussian envelope and white noise
    if (backgroundModelName == "Flat")
    {    
        predictFlat(predictions); 
        backgroundNameMatched = true;
    }
   
    if (backgroundNameMatched == false)
    {
        cerr << "Cannot match background model name with implemented ones. Quitting program." << endl;
        exit(EXIT_FAILURE);
    }
}










// BackgroundModel::predictThreeHarveyColor()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      three Harvey profiles and with colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictThreeHarveyColor(RefArrayXd predictions)
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












// BackgroundModel::predictThreeHarvey()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      three Harvey profiles and without colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictThreeHarvey(RefArrayXd predictions)
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


    // Add flat noise

    predictions += flatNoiseLevel;
}












// BackgroundModel::predictTwoHarveyColor()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      two Harvey profiles and with colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictTwoHarveyColor(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeNoise = configuringParameters(1);
    double frequencyNoise = configuringParameters(2);
    double amplitudeHarvey1 = configuringParameters(3);
    double frequencyHarvey1 = configuringParameters(4);
    double amplitudeHarvey2 = configuringParameters(5);
    double frequencyHarvey2 = configuringParameters(6);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));
    predictions += zeta*amplitudeHarvey2*amplitudeHarvey2/(frequencyHarvey2*(1.0 + (covariates/frequencyHarvey2).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise and colored noise components

    predictions += flatNoiseLevel;
    predictions += 2.0*Functions::PI*amplitudeNoise*amplitudeNoise/(frequencyNoise*(1.0 + (covariates/frequencyNoise).pow(2)));
}










// BackgroundModel::predictTwoHarvey()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      two Harvey profiles and without colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictTwoHarvey(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeHarvey1 = configuringParameters(1);
    double frequencyHarvey1 = configuringParameters(2);
    double amplitudeHarvey2 = configuringParameters(3);
    double frequencyHarvey2 = configuringParameters(4);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));
    predictions += zeta*amplitudeHarvey2*amplitudeHarvey2/(frequencyHarvey2*(1.0 + (covariates/frequencyHarvey2).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise

    predictions += flatNoiseLevel;
}












// BackgroundModel::predictOneHarveyColor()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      one Harvey profile and with colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictOneHarveyColor(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeNoise = configuringParameters(1);
    double frequencyNoise = configuringParameters(2);
    double amplitudeHarvey1 = configuringParameters(3);
    double frequencyHarvey1 = configuringParameters(4);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise and colored noise components

    predictions += flatNoiseLevel;
    predictions += 2.0*Functions::PI*amplitudeNoise*amplitudeNoise/(frequencyNoise*(1.0 + (covariates/frequencyNoise).pow(2)));
}











// BackgroundModel::predictOneHarvey()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      one Harvey profile and without colored noise. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictOneHarvey(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeHarvey1 = configuringParameters(1);
    double frequencyHarvey1 = configuringParameters(2);


    // Compute Harvey components and add them to the predictions

    double zeta = 2.0*sqrt(2.0)/Functions::PI;
    predictions = zeta*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (covariates/frequencyHarvey1).pow(4)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise

    predictions += flatNoiseLevel;
}











// BackgroundModel::predictOriginal()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      one Harvey profile and without colored noise. This uses the
//      original Harvey law, with an exponent set to 2.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictOriginal(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);
    double amplitudeHarvey1 = configuringParameters(1);
    double frequencyHarvey1 = configuringParameters(2);


    // Compute Harvey components and add them to the predictions

    predictions = 4.0*amplitudeHarvey1*amplitudeHarvey1/(frequencyHarvey1*(1.0 + (2*Functions::PI*covariates/frequencyHarvey1).pow(2)));

    
    // Modulate the model by the response function (apodization)
    
    predictions *= responseFunction;           


    // Add flat noise

    predictions += flatNoiseLevel;
}











// BackgroundModel::predictFlat()
//
// PURPOSE:
//      Builds the predictions from a Background model based on 
//      just the flat noise level. 
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//

void BackgroundModel::predictFlat(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);


    // Add flat noise

    predictions += flatNoiseLevel;
}












// BackgroundModel::getBackgroundModelName()
//
// PURPOSE:
//      Gets the protected data member backgroundModelName.
//
// OUTPUT:
//      A string containing the name of the background model used.
//

string BackgroundModel::getBackgroundModelName()
{
    return backgroundModelName;
}











// BackgroundModel::getNyquistFrequency()
//
// PURPOSE:
//      Gets the protected data member NyquistFrequency.
//
// OUTPUT:
//      A double containing the Nyquist frequency for the given cadence adopted.
//

double BackgroundModel::getNyquistFrequency()
{
    return NyquistFrequency;
}










// BackgroundModel::getConfiguringParameters()
//
// PURPOSE:
//      Gets the protected data member configuringParameters.
//
// OUTPUT:
//      An eigen array containing the configuring parameters of the
//      background model.
//

ArrayXd BackgroundModel::getConfiguringParameters()
{
    return configuringParameters;
}













// BackgroundModel::getResponseFunction()
//
// PURPOSE:
//      Gets the protected data member responseFunction.
//
// OUTPUT:
//      An eigen array containing the apodization response function for the signal of the input data
//

ArrayXd BackgroundModel::getResponseFunction()
{
    return responseFunction;
}
