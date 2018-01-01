#include "LorentzianRotationMixtureModel.h"


// LorentzianRotationMixtureModel::LorentzianRotationMixtureModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      NnonRadial:             the number of non radial p modes (to be fitted by a Lorentzian profile)
//      Nradial:                the number of radial modes (to be fitted by a Lorentzian profile)
//      cosi:                   the cosine of the inclination angle for the Sun (i = 90.0 degrees by default)
//      frequencyResolution:    the frequency bin size given by the dataset used
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

LorentzianRotationMixtureModel::LorentzianRotationMixtureModel(RefArrayXd const covariates, const int NmodesWithRotation, 
                                                               const int NmodesWithoutRotation, const double cosi, 
                                                               string angularDegreesFileName, BackgroundModel &backgroundModel)
: Model(covariates),
  NmodesWithRotation(NmodesWithRotation),
  NmodesWithoutRotation(NmodesWithoutRotation),
  cosi(cosi),
  angularDegreesFileName(angularDegreesFileName)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
    backgroundParameters = backgroundModel.getConfiguringParameters();


    // Read the angular degrees of the modes with rotational signature to be fitted

    readAngularDegreesFromFile(angularDegreesFileName);

    
    // Compute the mode visibilities for each angular degree and azimuthal number using the cosi given as input

    computeModeVisibility(); 
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
//      (iv) Mode central frequency (times the number of unresolved peaks)
//      (v) Mode profile height (times the number of unresolved peaks)
//

void LorentzianRotationMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());


    // Add a Lorentzian profile with split components for each mode for which the rotational effect is fitted

    for (int mode = 0; mode < NmodesWithRotation; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(4*mode);
        double amplitude = modelParameters(4*mode + 1);
        double linewidth = modelParameters(4*mode + 2);
        double rotationalSplitting = modelParameters(4*mode + 3);

        
        // Obtain the angular degree of this mode and find the corresponding starting element index in the visibility array
        
        double angularDegree = angularDegreeArray(mode);
        int startingIndex = static_cast<int>(angularDegree*(angularDegree - 1)/2 + angularDegree - 1);
   
        
        // Include first the m=0 component

        double amplitudeVisibilityComponent = modeAmplitudeVisibilityArray(startingIndex);
        Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, 
                                            amplitude*amplitudeVisibilityComponent, linewidth);
        predictions += singleModePrediction;

        
        // Then loop over the remaining 2*l components

        for (int azimuthalNumber = 1; azimuthalNumber <= angularDegree; azimuthalNumber++)
        {
            amplitudeVisibilityComponent = modeAmplitudeVisibilityArray(startingIndex + azimuthalNumber);       // m != 0
            
            // + m
            Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency + (rotationalSplitting*azimuthalNumber), 
                                                amplitude*amplitudeVisibilityComponent, linewidth);
            predictions += singleModePrediction;
            
            // - m
            Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency - (rotationalSplitting*azimuthalNumber), 
                                                amplitude*amplitudeVisibilityComponent, linewidth);
            predictions += singleModePrediction;
        }
        
    }


    // Add a Lorentzian profile for each mode for which no rotational effect is fitted

    for (int mode = 0; mode < NmodesWithoutRotation; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(4*NmodesWithRotation + 3*mode);
        double amplitude = modelParameters(4*NmodesWithRotation + 3*mode + 1);
        double linewidth = modelParameters(4*NmodesWithRotation + 3*mode + 2);

        Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
        predictions += singleModePrediction;
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

    if (angularDegreeArray.size() != NmodesWithRotation)
    {
        cerr << "The number of input angular degrees is not matching" << endl;
        cerr << "the number of non radial modes in the fit. Quitting program." << endl;
    }

}












// LorentzianRotationMixtureModel::computeModeVisibility()
//
// PURPOSE: 
//      Computes the visibilities of each split component using the cosine of the 
//      inclination angle defined by the constructor.
//
// OUTPUT:
//      void
//

void LorentzianRotationMixtureModel::computeModeVisibility()
{
    double modeVisibilityComponent;
    modeVisibility1.resize(2);
    modeVisibility2.resize(3);
    modeVisibility3.resize(4);
    modeVisibility4.resize(5);
    

    // l = 1    m = 0
    modeVisibilityComponent = cosi*cosi;
    modeVisibility1(0) = modeVisibilityComponent;
            
    // l = 1    m = -1 / +1
    modeVisibilityComponent = 0.5*(1.0 - cosi*cosi);
    modeVisibility1(1) = modeVisibilityComponent;
            

    // l = 2    m = 0
    modeVisibilityComponent = 0.25*pow((3*cosi*cosi - 1.0),2);
    modeVisibility2(0) = modeVisibilityComponent;
            
    // l = 2    m = -1 / + 1
    modeVisibilityComponent = (3.0/2.0)*cosi*cosi*(1.0 - cosi*cosi); 
    modeVisibility2(1) = modeVisibilityComponent;
            
    // l = 2    m = -2 / + 2
    modeVisibilityComponent = (3.0/8.0)*pow((1.0 - cosi*cosi),2); 
    modeVisibility2(2) = modeVisibilityComponent;
            
            
    // l = 3    m = 0
    modeVisibilityComponent = 0.25*cosi*cosi*pow((5.0*cosi*cosi - 3.0),2);
    modeVisibility3(0) = modeVisibilityComponent;
            
    // l = 3    m = -1 / +1
    modeVisibilityComponent = (1.0/12.0)*(9.0/4.0)*pow((1.0 - 5.0*cosi*cosi),2)*(1.0 - cosi*cosi);
    modeVisibility3(1) = modeVisibilityComponent;
            
    // l = 3    m = -2 / +2
    modeVisibilityComponent = 1.87500*cosi*cosi*pow((1.0 - cosi*cosi),2);
    modeVisibility3(2) = modeVisibilityComponent;
            
    // l = 3    m = -3 / +3
    modeVisibilityComponent = 0.3125*pow((1.0 - cosi*cosi),3);
    modeVisibility3(3) = modeVisibilityComponent;

            
    // l = 4    m = 0
    modeVisibilityComponent = 0.0156250*pow((35*pow(cosi,4) - 30*cosi*cosi + 3),2);
    modeVisibility4(0) = modeVisibilityComponent;
                
    // l = 4    m = -1 / +1
    modeVisibilityComponent = 0.312500*cosi*cosi*pow((3.0 - 7.0*cosi*cosi),2)*(1 - cosi*cosi);
    modeVisibility4(1) = modeVisibilityComponent;
            
    // l = 4    m = -2 / +2
    modeVisibilityComponent = 0.156250*pow((7*cosi*cosi - 1.0)*(1.0 - cosi*cosi),2);
    modeVisibility4(2) = modeVisibilityComponent;
            
    // l = 4    m = -3 / +3
    modeVisibilityComponent = 2.1875000*cosi*cosi*pow((1.0 - cosi*cosi),3);
    modeVisibility4(3) = modeVisibilityComponent;

    // l = 4    m = -4 / +4
    modeVisibilityComponent = 0.27343750*pow((1.0 - cosi*cosi),4);
    modeVisibility4(4) = modeVisibilityComponent;

    
    // Convert these visibilities into amplitude visibilities

    createModeAmplitudeVisibility();
}











// LorentzianRotationMixtureModel::createModeAmplitudeVisibility()
//
// PURPOSE: 
//      Stores the visibilities in amplitude for each angular degree and 
//      azimuthal number as obtained from the angularDegreeArray.
//
// OUTPUT:
//      void
//

void LorentzianRotationMixtureModel::createModeAmplitudeVisibility()
{
    int currentAngularDegree;
    int currentSize = 0;

    for (int i = 0; i < NmodesWithRotation; i++)
    {
        currentAngularDegree = static_cast<int>(angularDegreeArray(i));

        switch(currentAngularDegree)
        {
            case(1):          // l = 1 modes
                modeAmplitudeVisibilityArray.conservativeResize(currentSize + 2);
                modeAmplitudeVisibilityArray.segment(currentSize, 2) = modeVisibility1.sqrt();
                break;

            case(2):          // l = 2 modes
                modeAmplitudeVisibilityArray.conservativeResize(currentSize + 3);
                modeAmplitudeVisibilityArray.segment(currentSize, 3) = modeVisibility2.sqrt();
                break;

            case(3):          // l = 3 modes
                modeAmplitudeVisibilityArray.conservativeResize(currentSize + 4);
                modeAmplitudeVisibilityArray.segment(currentSize, 4) = modeVisibility3.sqrt();
                break;

            case(4):          // l = 4 modes
                modeAmplitudeVisibilityArray.conservativeResize(currentSize + 5);
                modeAmplitudeVisibilityArray.segment(currentSize, 5) = modeVisibility4.sqrt();
                break;
        }
       

        // Update the number of elements of the modeAmplitudeVisibilityArray for the next iteration

        currentSize = modeAmplitudeVisibilityArray.size();
    }
    

}
