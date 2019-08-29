#include "PeakTestSincBackgroundModel.h"


// PeakTestSincBackgroundModel::PeakTestSincBackgroundModel()
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

PeakTestSincBackgroundModel::PeakTestSincBackgroundModel(RefArrayXd const covariates, const double frequencyResolution, BackgroundModel &backgroundModel)
: Model(covariates),
frequencyResolution(frequencyResolution)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestSincBackgroundModel::PeakTestSincBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestSincBackgroundModel::~PeakTestSincBackgroundModel()
{

}










// PeakTestSincBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestSincBackground model to
//      construct a Sinc^2 profile on top of a varying background level.
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
//      (ii) Mode profile height
//      (iv) Noise factor to change the level of background
//

void PeakTestSincBackgroundModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);
    double noiseFactor = modelParameters(2);

    Functions::modeProfileSinc(singleModePrediction, covariates, centralFrequency, height, frequencyResolution);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction*noiseFactor;
}
