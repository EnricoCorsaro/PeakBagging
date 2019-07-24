#include "PeakTestBackgroundModel.h"


// PeakTestBackgroundModel::PeakTestBackgroundModel()
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

PeakTestBackgroundModel::PeakTestBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestBackgroundModel::PeakTestBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestBackgroundModel::~PeakTestBackgroundModel()
{

}










// PeakTestBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestBackground model to
//      construct a varying background level.
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
//      (ii) Noise factor to change the level of background
//

void PeakTestBackgroundModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double noiseFactor = modelParameters(0);


    // Add background component

    predictions += backgroundPrediction*noiseFactor;
}
