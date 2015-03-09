#include "LorentzianMixtureModel.h"


// LorentzianMixtureModel::LorentzianMixtureModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      NparametersPerType:     configuring numbers for the Peak Bagging model.
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

LorentzianMixtureModel::LorentzianMixtureModel(const RefArrayXd covariates, const int Npeaks, BackgroundModel &backgroundModel)
: Model(covariates),
  Npeaks(Npeaks)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
    backgroundParameters = backgroundModel.getConfiguringParameters();
}










// LorentzianMixtureModel::LorentzianMixtureModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianMixtureModel::~LorentzianMixtureModel()
{

}










// LorentzianMixtureModel::predict()
//
// PURPOSE:
//      Builds the predictions from a LorentzianMixture model based on a simple
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
//      (i) Mode central frequency (times the number of peaks)
//      (ii) Mode profile amplitude (times the number of peaks)
//      (iii) Mode profile linewidth (times the number of peaks)

void LorentzianMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    Nparameters = modelParameters.size();
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    for (int peak = 0; peak < Npeaks; ++peak)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(3*peak);
        double amplitude = modelParameters(3*peak + 1);
        double linewidth = modelParameters(3*peak + 2);

        Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
        predictions += singleModePrediction;
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}











