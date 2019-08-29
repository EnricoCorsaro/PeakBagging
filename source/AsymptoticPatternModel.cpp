#include "AsymptoticPatternModel.h"


// AsymptoticPatternModel::AsymptoticPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      linewidth:              the linewidth value to be used for the Lorentzian profile of radial modes
//      asymptoticParameters:   a string file containing the name of the list of asymptotic parameters
//                              to set up the sliding pattern.
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

AsymptoticPatternModel::AsymptoticPatternModel(RefArrayXd const covariates, const double linewidth, const string asymptoticParametersFileName, BackgroundModel &backgroundModel)
: Model(covariates),
linewidth(linewidth)
{
    // Set up asymptotic parameters

    readAsymptoticParametersFromFile(asymptoticParametersFileName);


    // Set up background model

    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// AsymptoticPatternModel::AsymptoticPatternModel()
//
// PURPOSE: 
//      Destructor.
//

AsymptoticPatternModel::~AsymptoticPatternModel()
{

}










// AsymptoticPatternModel::predict()
//
// PURPOSE:
//      Builds the predictions from a AsymptoticPatternModel model based to catch up
//      the correct position of the radial modes.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//      modelParameters:    one-dimensional array where each element
//                          contains the value of a free parameter of the model
//
// OUTPUT:
//      void
//
// NOTE:
//      The free parameters are to be given in the order
//      (i) Central radial mode frequency
//      (ii) Central radial mode height (as obtained from a smoothed PSD)
//

void AsymptoticPatternModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());
    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);

    
    // Add a Lorentzian profile for each mode in the sliding pattern
       
    // l = 0
    Functions::modeProfile(singleModePrediction, covariates, centralFrequency, height, linewidth);
    predictions += singleModePrediction;

    // l = 1
    Functions::modeProfile(singleModePrediction, covariates, centralFrequency + DeltaNu/2.0 - deltaNu01, height*dipoleToRadialHeightRatio, linewidth*dipoleToRadialLinewidthRatio);
    predictions += singleModePrediction;

    // l = 2
    Functions::modeProfile(singleModePrediction, covariates, centralFrequency - deltaNu02, height*quadrupoleToRadialHeightRatio, linewidth);
    predictions += singleModePrediction;


    for (int mode = 1; mode < (Norders-1)/2; ++mode)
    {
        // Add pattern to the right side 
        // l = 0
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu, height, linewidth);
        predictions += singleModePrediction;

        // l = 1
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01, height*dipoleToRadialHeightRatio, 
        linewidth*dipoleToRadialLinewidthRatio);
        predictions += singleModePrediction;

        // l = 2
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02, height*quadrupoleToRadialHeightRatio, linewidth);
        predictions += singleModePrediction; 
        
        // Add pattern to the left side
        // l = 0
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency - mode*DeltaNu, height, linewidth);
        predictions += singleModePrediction;

        // l = 1
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency - mode*DeltaNu + DeltaNu/2.0 - deltaNu01, height*dipoleToRadialHeightRatio, 
        linewidth*dipoleToRadialLinewidthRatio);
        predictions += singleModePrediction;

        // l = 2
        Functions::modeProfile(singleModePrediction, covariates, centralFrequency - mode*DeltaNu - deltaNu02, height*quadrupoleToRadialHeightRatio, linewidth);
        predictions += singleModePrediction; 
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}











// AsymptoticPatternModel::readAsymptoticParametersFromFile()
//
// PURPOSE:
//      Reads the asymptotic parameters to build the sliding pattern model from an input file
//      specified by its full path as a string. The values are then
//      stored into a vector which is a protected data member and that can
//      be therefore accessed from a derived class that implements the model.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void AsymptoticPatternModel::readAsymptoticParametersFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Ncols == 0.0)
    {
        cerr << "Wrong number of input asymptotic parameters for the asymptotic pattern model." << endl;
        cerr << "The input parameters requires are (1) N orders, (2) DeltaNu, (3) deltaNu02, (4) deltaNu01, " << endl;
        cerr << "(5) l=1/l=0 Height, (6) l=2/l=0 Height, (7) l=1/l=0 FWHM." << endl; 
        exit(EXIT_FAILURE);
    }

    ArrayXd asymptoticParameters(Nrows);
    asymptoticParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    
    inputFile.close();

    
    // Set up asymptotic parameters

    Norders = static_cast<int>(asymptoticParameters(0));
    DeltaNu = asymptoticParameters(1);
    deltaNu02 = asymptoticParameters(2);
    deltaNu01 = asymptoticParameters(3);
    dipoleToRadialHeightRatio = asymptoticParameters(4);
    quadrupoleToRadialHeightRatio = asymptoticParameters(5);
    dipoleToRadialLinewidthRatio = asymptoticParameters(6);
}
