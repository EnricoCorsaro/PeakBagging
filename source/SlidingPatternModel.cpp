#include "SlidingPatternModel.h"


// SlidingPatternModel::SlidingPatternModel()
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

SlidingPatternModel::SlidingPatternModel(RefArrayXd const covariates, const double linewidth, 
    const ArrayXd slidingPatternParameters, const string asymptoticParametersFileName, 
    const string gaussianEnvelopeParametersFileName, BackgroundModel &backgroundModel)
: Model(covariates),
linewidth(linewidth),
slidingPatternParameters(slidingPatternParameters)
{
    // Set up asymptotic parameters

    readAsymptoticParametersFromFile(asymptoticParametersFileName);
    readGaussianEnvelopeParametersFromFile(gaussianEnvelopeParametersFileName);
    includeDipoleModes = true;
    includeQuadrupoleModes = true;
    includeOctupoleModes = true;


    // Set up background model

    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// SlidingPatternModel::SlidingPatternModel()
//
// PURPOSE: 
//      Destructor.
//

SlidingPatternModel::~SlidingPatternModel()
{

}










// SlidingPatternModel::predict()
//
// PURPOSE:
//      Builds the predictions from a SlidingPatternModel model based to catch up
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
//      (iii) Small frequency separation deltaNu02
//      (iv) Small frequency separation deltaNu01
//      (v) Rotational splitting
//      (vi) cos i, with i the inclination angle
//

void SlidingPatternModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());
    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);

    double DeltaNu;
    double deltaNu02;
    double deltaNu01;
    double deltaNu13;
    double rotationalSplitting;
    double cosi;

    int NfreeParameters = modelParameters.size();

    if (slidingPatternParameters(0) != -99)
    {
        // DeltaNu is a fixed parameter

        DeltaNu = slidingPatternParameters(0);
    }
    else
    {
        // DeltaNu is a free parameter
        
        DeltaNu = modelParameters(2);
    }
 

    if (slidingPatternParameters(1) != -99.0)
    {
        // deltaNu02 is a fixed parameter
        
        if (slidingPatternParameters(1) == 99.0)
        {
            // Quadrupole mode excluded from the pattern
            
            includeQuadrupoleModes = false;
            deltaNu02 = 0;
        }
        else
        {
            // Quadrupole mode included in the pattern

            includeQuadrupoleModes = true;
            deltaNu02 = slidingPatternParameters(1);
        }
    }
    else
    {
        // deltaNu02 is a free parameter
        
        includeQuadrupoleModes = true;

        if (slidingPatternParameters(0) != -99)
        {
            // DeltaNu is a fixed parameter 

            deltaNu02 = modelParameters(2);
        }
        else
        {
            // DeltaNu is a fixed parameter 
            
            deltaNu02 = modelParameters(3);
        }
    }


    if (slidingPatternParameters(2) != -99.0)
    {
        // deltaNu01 is a fixed parameter

        if (slidingPatternParameters(2) == 99.0)
        {
            // Dipole mode excluded from the pattern

            includeDipoleModes = false;
            deltaNu01 = 0;
        }
        else
        {
            // Dipole mode included in the pattern

            includeDipoleModes = true;
            deltaNu01 = slidingPatternParameters(2);
        }
    }
    else
    {
        // deltaNu01 is a free parameter
        
        includeDipoleModes = true;

        if (slidingPatternParameters(0) != -99)
        {
            // DeltaNu is a fixed parameter

            if (slidingPatternParameters(1) != -99)
            {
                // deltaNu02 is a fixed parameter

                deltaNu01 = modelParameters(2);
            }
            else
            {
                // deltaNu02 is a free parameter

                deltaNu01 = modelParameters(3);
            }
        }
        else
        {
            // DeltaNu is a free parameter

            if (slidingPatternParameters(1) != -99)
            {
                // deltaNu02 is a fixed parameter
                
                deltaNu01 = modelParameters(3);
            }
            else
            {
                // deltaNu02 is a free parameter
                
                deltaNu01 = modelParameters(4);
            }
        }
    }


    if (slidingPatternParameters(3) != -99.0)
    {
        // deltaNu03 is a fixed parameter
        
        if (slidingPatternParameters(3) == 99.0)
        {
            // Octupole mode excluded from the pattern
            
            includeOctupoleModes = false;
            deltaNu13 = 0;
        }
        else
        {
            // Octupole mode included in the pattern

            includeOctupoleModes = true;
            deltaNu13 = slidingPatternParameters(3);
        }
    }
    else
    {
        // deltaNu03 is a free parameter
        
        includeOctupoleModes = true;

        if (slidingPatternParameters(4) == -99)
        {
            // Rotation is included

            deltaNu13 = modelParameters(NfreeParameters-3);
        }
        else
        {
            // Rotation is not included

            deltaNu13 = modelParameters(NfreeParameters-1);
        }
    }


    if (slidingPatternParameters(4) == -99.0)
    {
        // Rotational splitting and cosi are free parameters

        rotationalSplitting = modelParameters(NfreeParameters-2);
        cosi = modelParameters(NfreeParameters-1);
    }


    // Add a Lorentzian profile for each mode in the sliding pattern

    for (int mode = -(Norders-1)/2; mode <= (Norders-1)/2; ++mode)
    {
        // l = 0

        Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu, height, linewidth);
        predictions += singleModePrediction;

        if (slidingPatternParameters(3) == -99.0)
        {
            if (includeDipoleModes)
            {
                // l = 1, m = 0

                modeVisibilityComponent = cosi*cosi;
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01, 
                    height*dipoleToRadialHeightRatio*modeVisibilityComponent, linewidth*dipoleToRadialLinewidthRatio);
                predictions += singleModePrediction;
               
                // l = 1, m = + 1 / m = -1
                
                modeVisibilityComponent = 0.5*(1.0 - cosi*cosi);
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01 + rotationalSplitting, 
                    height*dipoleToRadialHeightRatio*modeVisibilityComponent, linewidth*dipoleToRadialLinewidthRatio);
                predictions += singleModePrediction;
                
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01 - rotationalSplitting, 
                    height*dipoleToRadialHeightRatio*modeVisibilityComponent, linewidth*dipoleToRadialLinewidthRatio);
                predictions += singleModePrediction;
            }


            if (includeQuadrupoleModes)
            {
                // l = 2, m = 0

                modeVisibilityComponent = 0.25*pow((3*cosi*cosi - 1.0),2);
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02, 
                    height*quadrupoleToRadialHeightRatio*modeVisibilityComponent, linewidth);
                predictions += singleModePrediction; 
           

                // l = 2    m = + 1 / - 1
                
                modeVisibilityComponent = (3.0/2.0)*cosi*cosi*(1.0 - cosi*cosi); 
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02 + rotationalSplitting, 
                    height*quadrupoleToRadialHeightRatio*modeVisibilityComponent, linewidth);
                predictions += singleModePrediction; 
      
                  Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02 - rotationalSplitting, 
                    height*quadrupoleToRadialHeightRatio*modeVisibilityComponent, linewidth);
                predictions += singleModePrediction; 

       
                // l = 2    m = + 2 / - 2
                
                modeVisibilityComponent = (3.0/8.0)*pow((1.0 - cosi*cosi),2); 
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02 + 2*rotationalSplitting, 
                    height*quadrupoleToRadialHeightRatio*modeVisibilityComponent, linewidth);
                predictions += singleModePrediction; 

                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02 - 2*rotationalSplitting, 
                    height*quadrupoleToRadialHeightRatio*modeVisibilityComponent, linewidth);
                predictions += singleModePrediction; 
            }
        }
        else
        {
            if (includeDipoleModes)
            {
                // l = 1

                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01, 
                    height*dipoleToRadialHeightRatio, linewidth*dipoleToRadialLinewidthRatio);
                predictions += singleModePrediction;
            }
            
            if (includeQuadrupoleModes)
            {
                // l = 2
                
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu - deltaNu02, height*quadrupoleToRadialHeightRatio, linewidth);
                predictions += singleModePrediction; 
            }
        }

        if (includeOctupoleModes)
        {
            // l = 3
            
            Functions::modeProfile(singleModePrediction, covariates, centralFrequency + mode*DeltaNu + DeltaNu/2.0 - deltaNu01 - deltaNu13, height*octupoleToRadialHeightRatio, linewidth);
            predictions += singleModePrediction; 
        }
    }


    // Adjust by Gaussian envelope modulation and by apodization of the signal
    
    predictions *= exp(-1.0*(nuMax - covariates)*(nuMax - covariates)/(2.0 * standardDeviation * standardDeviation));
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}











// SlidingPatternModel::readAsymptoticParametersFromFile()
//
// PURPOSE:
//      Reads the asymptotic parameters from an input file
//      specified by its full path as a string. The values are then
//      stored into private data members.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void SlidingPatternModel::readAsymptoticParametersFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Ncols == 0.0)
    {
        cerr << "Wrong number of input asymptotic parameters for the asymptotic pattern model." << endl;
        cerr << "The input parameters required are (1) N orders, " << endl;
        cerr << "(2) l=1/l=0 Height, (3) l=2/l=0 Height, (4) l=3/l=0 Height, (5) l=1/l=0 FWHM, ." << endl; 
        exit(EXIT_FAILURE);
    }

    ArrayXd asymptoticParameters(Nrows);
    asymptoticParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    
    inputFile.close();

    
    // Set up asymptotic parameters

    Norders = static_cast<int>(asymptoticParameters(0));
    dipoleToRadialHeightRatio = asymptoticParameters(1);
    quadrupoleToRadialHeightRatio = asymptoticParameters(2);
    octupoleToRadialHeightRatio = asymptoticParameters(3);
    dipoleToRadialLinewidthRatio = asymptoticParameters(4);
}














// SlidingPatternModel::readGaussianEnvelopeParametersFromFile()
//
// PURPOSE:
//      Reads the Gaussian envelope parameters from an input file
//      specified by its full path as a string. The values are then
//      stored into private data members.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void SlidingPatternModel::readGaussianEnvelopeParametersFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    if (Ncols == 0.0)
    {
        cerr << "Wrong number of input Gaussian envelope parameters for the asymptotic pattern model." << endl;
        cerr << "The input parameters required are (1) Height, (2) nuMax, (3) Gaussian standard deviation." << endl;
        exit(EXIT_FAILURE);
    }

    ArrayXd gaussianEnvelopeParameters(Nrows);
    gaussianEnvelopeParameters = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

    
    // Set up asymptotic parameters

    nuMax = gaussianEnvelopeParameters(1);
    standardDeviation = gaussianEnvelopeParameters(2);
}
