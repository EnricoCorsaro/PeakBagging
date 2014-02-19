#include "BackgroundModel.h"


// BackgroundModel::BackgroundModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      NparametersPerType:     configuring numbers for the Peak Bagging model
//      nuMax:                  the frequency of maximum power excess
//

BackgroundModel::BackgroundModel(const RefArrayXd covariates)
: Model(covariates)
{
}










// BackgroundModel::BackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

BackgroundModel::~BackgroundModel()
{

}











// BackgroundModel::readConfiguringParametersFromFile()
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
    
    configuringParameters.resize(Nrows);
    configuringParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);

    inputFile.close();
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
