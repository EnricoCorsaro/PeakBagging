#include "PeakTestLorentzianModel.h"


// PeakTestLorentzianModel::PeakTestLorentzianModel()
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

PeakTestLorentzianModel::PeakTestLorentzianModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestLorentzianModel::PeakTestLorentzianModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestLorentzianModel::~PeakTestLorentzianModel()
{

}










// PeakTestLorentzianModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestLorentzian model to
//      construct Lorentzian profile on top of a fixed background level.
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
//

void PeakTestLorentzianModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double amplitude = modelParameters(1);
    double linewidth = modelParameters(2);

    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
