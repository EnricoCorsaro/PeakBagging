#include "PeakTestLorentzianBackgroundModel.h"


// PeakTestLorentzianBackgroundModel::PeakTestLorentzianBackgroundModel()
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

PeakTestLorentzianBackgroundModel::PeakTestLorentzianBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestLorentzianBackgroundModel::PeakTestLorentzianBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestLorentzianBackgroundModel::~PeakTestLorentzianBackgroundModel()
{

}










// PeakTestLorentzianBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestLorentzianBackground model to
//      construct a Lorentzian profile on top of a varying background level.
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
//      (iv) Noise factor to change the level of background
//

void PeakTestLorentzianBackgroundModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double amplitude = modelParameters(1);
    double linewidth = modelParameters(2);
    double noiseFactor = modelParameters(3);

    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction*noiseFactor;
}
