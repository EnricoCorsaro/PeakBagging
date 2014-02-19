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

LorentzianMixtureModel::LorentzianMixtureModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, BackgroundModel &backgroundModel)
: Model(covariates),
  NprofileParameters(NparametersPerType[0]),
  Nmodes(NparametersPerType[1])
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
//      (i) Mode central frequency (times the number of modes)
//      (ii) Mode profile amplitude (times the number of modes)
//      (iii) Mode profile linewidth (times the number of modes)

void LorentzianMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    Nparameters = modelParameters.size();
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    for (int mode = 0; mode < Nmodes; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        double centralFrequency = modelParameters(NprofileParameters*mode);
        double amplitude = modelParameters(NprofileParameters*mode + 1);
        double linewidth = modelParameters(NprofileParameters*mode + 2);

        Functions::modeProfileWithAmplitude(singleModePrediction, covariates, centralFrequency, amplitude, linewidth);
        predictions += singleModePrediction;
    }


    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}











