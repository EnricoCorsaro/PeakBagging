#include "PeakTestLorentzianRotationModel.h"


// PeakTestLorentzianRotationModel::PeakTestLorentzianRotationModel()
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

PeakTestLorentzianRotationModel::PeakTestLorentzianRotationModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestLorentzianRotationModel::PeakTestLorentzianRotationModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestLorentzianRotationModel::~PeakTestLorentzianRotationModel()
{

}










// PeakTestLorentzianRotationModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestLorentzianRotation model to
//      construct a Lorentzian profile with rotationally split components
//      on top of a fixed background level. Note that the Lorentzian
//      profile here considered is assumed to be only a ell = 1 mode, hence
//      having at most three rotationally split components.
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
//      (i) Mode central frequency and data frequency range
//      (ii) Mode profile amplitude
//      (iii) Mode profile linewidth
//      (iv) Rotational splitting
//      (v) Inclination angle of stellar rotation axis
//

void PeakTestLorentzianRotationModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double amplitude = modelParameters(1);
    double linewidth = modelParameters(2);
    double rotationalSplitting = modelParameters(3);
    double cosi = modelParameters(4);


    // Include first the m=0 component

    double modeVisibilityComponent = cosi*cosi;
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, 
                                        amplitude*sqrt(modeVisibilityComponent), linewidth);
    predictions += singleModePrediction;

    
    // Then add the two remaining components
    // m = +1

    modeVisibilityComponent = 0.5*(1.0 - cosi*cosi);
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency + rotationalSplitting, 
                                        amplitude*sqrt(modeVisibilityComponent), linewidth);
    predictions += singleModePrediction;
    
    
    // m = -1
    
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency - rotationalSplitting, 
                                        amplitude*sqrt(modeVisibilityComponent), linewidth);
    predictions += singleModePrediction;

        
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
