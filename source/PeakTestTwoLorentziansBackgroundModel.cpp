#include "PeakTestTwoLorentziansBackgroundModel.h"


// PeakTestTwoLorentziansBackgroundModel::PeakTestTwoLorentziansBackgroundModel()
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

PeakTestTwoLorentziansBackgroundModel::PeakTestTwoLorentziansBackgroundModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestTwoLorentziansBackgroundModel::PeakTestTwoLorentziansBackgroundModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestTwoLorentziansBackgroundModel::~PeakTestTwoLorentziansBackgroundModel()
{

}










// PeakTestTwoLorentziansBackgroundModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestTwoLorentziansBackground model to
//      construct two Lorentzian profiles on top of a varying background level.
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
//      (i) First mode central frequency
//      (ii) First mode profile amplitude
//      (iii) First mode profile linewidth
//      (iv) Second mode central frequency
//      (v) Second mode profile amplitude
//      (vi) Second mode profile linewidth
//      (vii) Noise factor to change the level of background
//

void PeakTestTwoLorentziansBackgroundModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency1 = modelParameters(0);
    double amplitude1 = modelParameters(1);
    double linewidth1 = modelParameters(2);
    double centralFrequency2 = modelParameters(3);
    double amplitude2 = modelParameters(4);
    double linewidth2 = modelParameters(5);
    double noiseFactor = modelParameters(6);

    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency1, amplitude1, linewidth1);
    predictions += singleModePrediction;
    
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency2, amplitude2, linewidth2);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction*noiseFactor;
}
