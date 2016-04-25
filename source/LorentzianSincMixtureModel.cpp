#include "LorentzianSincMixtureModel.h"


// LorentzianSincMixtureModel::LorentzianSincMixtureModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      Nresolved:              the number of resolved modes (to be fitted by a Lorentzian profile)
//      Nunresolved:            the number of unresolved modes (to be fitted by a sinc-square profile)
//      frequencyResolution:    the frequency bin size given by the dataset used
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

LorentzianSincMixtureModel::LorentzianSincMixtureModel(RefArrayXd const covariates, const int Nresolved, const int Nunresolved, 
                                                       const double frequencyResolution, BackgroundModel &backgroundModel)
: Model(covariates),
  Nresolved(Nresolved),
  Nunresolved(Nunresolved),
  frequencyResolution(frequencyResolution)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
    backgroundParameters = backgroundModel.getConfiguringParameters();
}










// LorentzianSincMixtureModel::LorentzianSincMixtureModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianSincMixtureModel::~LorentzianSincMixtureModel()
{

}










// LorentzianSincMixtureModel::predict()
//
// PURPOSE:
//      Builds the predictions from a LorentzianSincMixture model based on a simple
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
//      (iv) Mode central frequency (times the number of unresolved peaks)
//      (v) Mode profile height (times the number of unresolved peaks)
//

void LorentzianSincMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
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

    
    // Add sinc-square profile for each unresolved mixed mode
    
    for (int mode = 0; mode < Nunresolved; ++mode)
    {
        double centralFrequency = modelParameters(3*Nresolved + 2*mode);
        double height = modelParameters(3*Nresolved + 2*mode + 1);            // mixedModesHeights(mode);
        
        Functions::modeProfileSinc(singleModePrediction, covariates, centralFrequency, height, frequencyResolution);
        predictions += singleModePrediction;
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}

