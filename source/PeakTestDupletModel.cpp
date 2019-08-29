#include "PeakTestDupletModel.h"


// PeakTestDupletModel::PeakTestDupletModel()
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

PeakTestDupletModel::PeakTestDupletModel(RefArrayXd const covariates, BackgroundModel &backgroundModel)
: Model(covariates)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// PeakTestDupletModel::PeakTestDupletModel()
//
// PURPOSE: 
//      Destructor.
//

PeakTestDupletModel::~PeakTestDupletModel()
{

}










// PeakTestDupletModel::predict()
//
// PURPOSE:
//      Builds the predictions from a PeakTestDuplet model to
//      construct two Lorentzian profiles in the form of a duplet on top of a constant background level.
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
//      (i) Left mode central frequency
//      (ii) Left mode profile amplitude
//      (iii) Left mode profile linewidth
//      (iv) Duplet frequency splitting
//      (v) Right mode profile amplitude
//      (vi) Right mode profile linewidth
//

void PeakTestDupletModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    // Initialize parameters of the oscillation mode

    double centralFrequency1 = modelParameters(0);
    double amplitude1 = modelParameters(1);
    double linewidth1 = modelParameters(2);
    double dupletSplitting = modelParameters(3);
    double amplitude2 = modelParameters(4);
    double linewidth2 = modelParameters(5);

    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency1, amplitude1, linewidth1);
    predictions += singleModePrediction;
    
    Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency1 + dupletSplitting, amplitude2, linewidth2);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
