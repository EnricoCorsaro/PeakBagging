#include "LorentzianMixtureModel.h"


// LorentzianMixtureModel::LorentzianMixtureModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      Nresolved:              the number of resolved modes (to be fitted by a Lorentzian profile)
//      frequencyResolution:    the frequency bin size given by the dataset used
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

LorentzianMixtureModel::LorentzianMixtureModel(RefArrayXd const covariates, const int Nresolved, BackgroundModel &backgroundModel)
: Model(covariates),
  Nresolved(Nresolved)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
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
//      (i) Mode central frequency (times the number of resolved peaks)
//      (ii) Mode profile amplitude (times the number of resolved peaks)
//      (iii) Mode profile linewidth (times the number of resolved peaks)
//

void LorentzianMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    
    // Add a Lorentzian profile for each resolved mode (both p modes and resolved mixed modes)

    for (int mode = 0; mode < Nresolved; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(3*mode);
        double amplitude = modelParameters(3*mode + 1);
        double linewidth = modelParameters(3*mode + 2);

        Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
        predictions += singleModePrediction;
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}

