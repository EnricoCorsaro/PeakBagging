#include "DoubleLorentzianFixedLinewidthModel.h"


// DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      linewidth:              the linewidth value to be used for the Lorentzian profile
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel(RefArrayXd const covariates, const double linewidth,
BackgroundModel &backgroundModel, const string modeVisibilityFileName)
: Model(covariates),
  linewidth(linewidth)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();

    // Set up mode visibility from input file

    readModeVisibilityFromFile(modeVisibilityFileName);
}










// DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel()
//
// PURPOSE: 
//      Destructor.
//

DoubleLorentzianFixedLinewidthModel::~DoubleLorentzianFixedLinewidthModel()
{

}










// DoubleLorentzianFixedLinewidthModel::predict()
//
// PURPOSE:
//      Builds the predictions from a DoubleLorentzianFixedLinewidth model to
//      construct a multi-modal posterior distribution (islands model).
//      
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
//      (i) Mode central frequency
//      (ii) Mode profile height
//      (iii) Small frequency separation d02
//

void DoubleLorentzianFixedLinewidthModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    
    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);
    double smallSeparation = modelParameters(2);


    // Add firt profile with reference frequency given by the centralFrequency

    Functions::modeProfile(singleModePrediction, covariates, centralFrequency, height, linewidth);
    predictions += singleModePrediction;
    

    // Add second profile with reference frequency given by the centralFrequency - smallSeparation
   
    Functions::modeProfile(singleModePrediction, covariates, centralFrequency - smallSeparation, height*quadrupoleToRadialHeightRatio, linewidth);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}









// DoubleLorentzianFixedLinewidthModel::readModeVisibilityFromFile()
//
// PURPOSE:
//      Reads the mode visibility from the same input file containing
//      the asymptotic parameters for the sliding pattern fit. The value is then
//      stored into private data members.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void DoubleLorentzianFixedLinewidthModel::readModeVisibilityFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Ncols == 0.0)
    {
        cerr << "Wrong number of input asymptotic parameters." << endl;
        cerr << "The input parameters required are (1) N orders, " << endl;
        cerr << "(2) l=1/l=0 Height, (3) l=2/l=0 Height, (4) l=3/l=0 Height, (5) l=1/l=0 FWHM, ." << endl; 
        exit(EXIT_FAILURE);
    }

    ArrayXd asymptoticParameters(Nrows);
    asymptoticParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    
    inputFile.close();

    
    // Set up mode visibility for l=2/l=0

    quadrupoleToRadialHeightRatio = asymptoticParameters(2);
}