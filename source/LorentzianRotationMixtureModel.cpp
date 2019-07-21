#include "LorentzianRotationMixtureModel.h"


// LorentzianRotationMixtureModel::LorentzianRotationMixtureModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      angularDegreesFileName: a string containing the filename of the input list of angular degrees, one
//                              for each mode considered in the fit.
//      Nresolved:              total number of modes.
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

LorentzianRotationMixtureModel::LorentzianRotationMixtureModel(RefArrayXd const covariates, const int Nresolved, 
                                                               string angularDegreesFileName, BackgroundModel &backgroundModel)
: Model(covariates),
  Nresolved(Nresolved),
  angularDegreesFileName(angularDegreesFileName)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();


    // Read the angular degrees of the modes with rotational signature to be fitted

    readAngularDegreesFromFile(angularDegreesFileName);
}
    









// LorentzianRotationMixtureModel::LorentzianRotationMixtureModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianRotationMixtureModel::~LorentzianRotationMixtureModel()
{

}











// LorentzianRotationMixtureModel::predict()
//
// PURPOSE:
//      Builds the predictions from a LorentzianRotationMixture model based on a simple
//      inputting of the central frequencies for all the desired modes to be fitted.
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
//      (i) Mode central frequency (times the number of resolved peaks)
//      (ii) Mode profile amplitude (times the number of resolved peaks)
//      (iii) Mode profile linewidth (times the number of resolved peaks)
//      (iv) rotational splitting
//      (v) cos i, with i the inclination angle of the rotation axis
//

void LorentzianRotationMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());
    double rotationalSplitting = modelParameters(modelParameters.size() - 2);
    double cosi = modelParameters(modelParameters.size() - 1);

    // Add a Lorentzian profile with split components for each mode for which the rotational effect is fitted (only l=1 for this model)

    for (int mode = 0; mode < Nresolved; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(3*mode);
        double amplitude = modelParameters(3*mode + 1);
        double linewidth = modelParameters(3*mode + 2);

        
        // Obtain the angular degree of this mode
        
        double angularDegree = angularDegreeArray(mode);
        ArrayXd modeAmplitudeVisibility = computeModeVisibility(angularDegree, cosi);
  
        if (angularDegree != 0)
        {
            // Include first the m=0 component

            Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, 
                                                amplitude*modeAmplitudeVisibility(0), linewidth);
            predictions += singleModePrediction;

            
            // Then loop over the remaining 2*ell components

            for (int azimuthalNumber = 1; azimuthalNumber <= angularDegree; azimuthalNumber++)
            {
                // + m
                Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency + (rotationalSplitting*azimuthalNumber), 
                                                    amplitude*modeAmplitudeVisibility(azimuthalNumber), linewidth);
                predictions += singleModePrediction;
                
                // - m
                Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency - (rotationalSplitting*azimuthalNumber), 
                                                    amplitude*modeAmplitudeVisibility(azimuthalNumber), linewidth);
                predictions += singleModePrediction;
            }
        }
        else
        {
            // Add a Lorentzian profile for each mode for which no rotational effect is fitted
            
            Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
            predictions += singleModePrediction;
        }
        
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}











// LorentzianRotationMixturedModel::readAngularDegreesFromFile()
//
// PURPOSE:
//      Reads the angular degrees of the non radial modes to be fit in the peak bagging model from an input file
//      specified by its full path as a string. The values are then stored into an Eigen array.
//
// INPUT:
//      inputFileName:      a string specifying the full path (filename included) of the input file to read.
//
// OUTPUT:
//      void
//

void LorentzianRotationMixtureModel::readAngularDegreesFromFile(const string inputFileName)
{
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);

    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    
    angularDegreeArray.resize(Nrows);
    angularDegreeArray = File::arrayXXdFromFile(inputFile, Nrows, Ncols);

    inputFile.close();

    assert(angularDegreeArray.size() == Nresolved);
}












// LorentzianRotationMixtureModel::computeModeVisibility()
//
// PURPOSE: 
//      Computes the visibilities of each split component using the cosine of the 
//      inclination angle as input and the corresponding angular degree of interest.
//
// OUTPUT:
//      An array of mode visibilities for amplitudes (hence original mode visibilities under squareroot).
//

ArrayXd LorentzianRotationMixtureModel::computeModeVisibility(const int angularDegree, const double cosi)
{
    double modeVisibilityComponent;
    ArrayXd modeVisibility; 

    if (angularDegree == 1)
    {
        modeVisibility.resize(2); 
        
        // l = 1    m = 0
        modeVisibilityComponent = cosi*cosi;
        modeVisibility(0) = modeVisibilityComponent;
                
        // l = 1    m = -1 / +1
        modeVisibilityComponent = 0.5*(1.0 - cosi*cosi);
        modeVisibility(1) = modeVisibilityComponent;

    }       

    if (angularDegree == 2)
    {
        modeVisibility.resize(3); 
        
        // l = 2    m = 0
        modeVisibilityComponent = 0.25*pow((3*cosi*cosi - 1.0),2);
        modeVisibility(0) = modeVisibilityComponent;
                
        // l = 2    m = -1 / + 1
        modeVisibilityComponent = (3.0/2.0)*cosi*cosi*(1.0 - cosi*cosi); 
        modeVisibility(1) = modeVisibilityComponent;
                
        // l = 2    m = -2 / + 2
        modeVisibilityComponent = (3.0/8.0)*pow((1.0 - cosi*cosi),2); 
        modeVisibility(2) = modeVisibilityComponent;
        
    }       
    
    if (angularDegree == 3)
    {
        modeVisibility.resize(4); 
        
        // l = 3    m = 0
        modeVisibilityComponent = 0.25*cosi*cosi*pow((5.0*cosi*cosi - 3.0),2);
        modeVisibility(0) = modeVisibilityComponent;
                
        // l = 3    m = -1 / +1
        modeVisibilityComponent = (1.0/12.0)*(9.0/4.0)*pow((1.0 - 5.0*cosi*cosi),2)*(1.0 - cosi*cosi);
        modeVisibility(1) = modeVisibilityComponent;
                
        // l = 3    m = -2 / +2
        modeVisibilityComponent = 1.87500*cosi*cosi*pow((1.0 - cosi*cosi),2);
        modeVisibility(2) = modeVisibilityComponent;
                
        // l = 3    m = -3 / +3
        modeVisibilityComponent = 0.3125*pow((1.0 - cosi*cosi),3);
        modeVisibility(3) = modeVisibilityComponent;
        
    }
     
    if (angularDegree == 4)
    {
        modeVisibility.resize(5);

        // l = 4    m = 0
        modeVisibilityComponent = 0.0156250*pow((35*pow(cosi,4) - 30*cosi*cosi + 3),2);
        modeVisibility(0) = modeVisibilityComponent;
                    
        // l = 4    m = -1 / +1
        modeVisibilityComponent = 0.312500*cosi*cosi*pow((3.0 - 7.0*cosi*cosi),2)*(1.0 - cosi*cosi);
        modeVisibility(1) = modeVisibilityComponent;
                
        // l = 4    m = -2 / +2
        modeVisibilityComponent = 0.156250*pow((7*cosi*cosi - 1.0)*(1.0 - cosi*cosi),2);
        modeVisibility(2) = modeVisibilityComponent;
                
        // l = 4    m = -3 / +3
        modeVisibilityComponent = 2.1875000*cosi*cosi*pow((1.0 - cosi*cosi),3);
        modeVisibility(3) = modeVisibilityComponent;

        // l = 4    m = -4 / +4
        modeVisibilityComponent = 0.27343750*pow((1.0 - cosi*cosi),4);
        modeVisibility(4) = modeVisibilityComponent;
        
    }
        
    return modeVisibility.sqrt();
}

