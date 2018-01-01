#include "ConstantBackgroundModel.h"


// ConstantBackgroundModel::ConstantBackgroundModel()
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

ConstantBackgroundModel::ConstantBackgroundModel(const RefArrayXd covariates, const string inputNyquistFrequencyFileName)
: BackgroundModel(covariates)
{
    // NyquistFrequency = 8496.355743094671     muHz     // Kepler SC
    // NyquistFrequency = 283.2116656017908     muHz     // Kepler LC

    readNyquistFrequencyFromFile(inputNyquistFrequencyFileName);

    ArrayXd sincFunctionArgument = (Functions::PI / 2.0) * covariates / NyquistFrequency;
    responseFunction = (sincFunctionArgument.sin() / sincFunctionArgument).square();
}










// ConstantBackgroundModel::ConstantBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

ConstantBackgroundModel::~ConstantBackgroundModel()
{

}










// ConstantBackgroundModel::readConfiguringParametersFromFile()
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

void ConstantBackgroundModel::readConfiguringParametersFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Nrows != 1)
    {
        cerr << "Wrong number of parameters that define the standard background model (1 parameter expected)." << endl;
        exit(EXIT_FAILURE);
    }

    configuringParameters.resize(Nrows);
    configuringParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);

    inputFile.close();
}








// ConstantBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a Background model based on a configuration
//      file that contains all the free parameters of the given model.
//      This is an overloaded function whose implementation is that of a 
//      constant background model accounting only for a white noise level.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//
// OUTPUT:
//      void
//

void ConstantBackgroundModel::predict(RefArrayXd predictions)
{
    // Initialize global parameters

    double flatNoiseLevel = configuringParameters(0);


    // Add flat noise level component

    predictions += flatNoiseLevel;
    
}
