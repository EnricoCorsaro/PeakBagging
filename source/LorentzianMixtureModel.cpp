#include "LorentzianMixtureModel.h"


// LorentzianMixtureModel::RegularPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      NparametersPerType:     configuring numbers for the Peak Bagging model
//      nuMax:                  the frequency of maximum power excess
//

LorentzianMixtureModel::RegularPatternModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, double nuMax)
: Model(covariates),
  NglobalParameters(NparametersPerType[0]),
  NprofileParameters(NparametersPerType[1]),
  NradialOrders(NparametersPerType[2]),
  NangularDegrees(NparametersPerType[3]),
  NpressureModes(NparametersPerType[2]*NparametersPerType[3]),
  NmixedModes(NparametersPerType[4]),
  nuMax(nuMax)
{
    relativeRadialOrders.resize(NradialOrders);
    int startingRadialOrder = 0;

    if (static_cast<int>(fmod(NradialOrders,2)) == 0)
    {
        startingRadialOrder = -((NradialOrders/2) - 1);
    }
    else 
        if ((static_cast<int>(fmod(NradialOrders,2)) != 0) && (NradialOrders != 1))
        {
            startingRadialOrder = -((NradialOrders - 1)/2);
        }

    for (int radialOrder = 0; radialOrder < NradialOrders; ++radialOrder)
    {
        relativeRadialOrders[radialOrder] = startingRadialOrder + radialOrder;
    }

}










// LorentzianMixtureModel::RegularPatternModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianMixtureModel::~RegularPatternModel()
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
//      (2) White noise background (flat noise level)
//      (3) Mode central frequency
//      (5) Mode profile ln(height)
//      (6) Mode profile linewidth

void LorentzianMixtureModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    Nparameters = modelParameters.size();
    assert(Nparameters == (NglobalParameters + NprofileParameters*Nmodes));
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());

    for (int mode = 0; mode < Nmodes; ++mode)
    {
        // Initialize parameters of current mode with proper access to elements of total array of free parameters

        centralFrequency = modelParameters(NglobalParameters + mode);
        naturalLogarithmOfHeight = modelParameters(NglobalParameters + Nmodes + mode);
        linewidth = modelParameters(NglobalParameters + 2*Nmodes + mode);


        // Compute the prediction for the mode, provided the mode frequency is falling in the observed frequency range

        if ((centralFrequency > covariates.minCoeff()) && (centralFrequency < covariates.maxCoeff()))
        {
                Functions::modeProfile(singleModePrediction, covariates, centralFrequency, 
                                       exp(naturalLogarithmOfHeight), linewidth);
                predictions += singleModePrediction;
        }
        else
        {
            continue;
        }
    }


    // Add flat noise level component

    predictions += modelParameters(0);           
}











