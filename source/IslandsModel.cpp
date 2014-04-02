#include "IslandsModel.h"


// IslandsModel::IslandsModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

IslandsModel::IslandsModel(const RefArrayXd covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
    backgroundParameters = backgroundModel.getConfiguringParameters();
}










// IslandsModel::IslandsModel()
//
// PURPOSE: 
//      Destructor.
//

IslandsModel::~IslandsModel()
{

}










// IslandsModel::predict()
//
// PURPOSE:
//      Builds the predictions from an Islands model.
//      This implies the frequency position of the peaks is fitted with 
//      three parameters, which should then show a mode for each observed peak
//      in the parameter space.
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
//      (1) Position of the l = 1 mode (or modes)
//      (2) Spacing from l = 1 to adjacent l = 0 mode
//      (3) Spacing from l = 0 to adjacent l = 2 mode
//      (4) Amplitude for each mode
//      (6) Linewidth for each mode

void IslandsModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    Nparameters = modelParameters.size();
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());
    
    double dipoleModeFrequency = modelParameters(0);
    double spacing10 = modelParameters(1);
    double spacing20 = modelParameters(2);
    double amplitudeDipole = modelParameters(3);
    double linewidthDipole = modelParameters(4);
    double amplitudeQuadrupole = modelParameters(5);
    double linewidthQuadrupole = modelParameters(6);
    double amplitudeRadial = modelParameters(7);
    double linewidthRadial = modelParameters(8);

    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, dipoleModeFrequency, 
                                        amplitudeDipole, linewidthDipole);                
    predictions += singleModePrediction;

    double radialModeFrequency = dipoleModeFrequency + spacing10;
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, radialModeFrequency, 
                                        amplitudeRadial, linewidthRadial);                
    predictions += singleModePrediction;

    double quadrupoleModeFrequency = dipoleModeFrequency + spacing10 - spacing20;
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, quadrupoleModeFrequency, 
                                        amplitudeQuadrupole, linewidthQuadrupole);                
    predictions += singleModePrediction;


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
