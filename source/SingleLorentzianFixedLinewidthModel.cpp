#include "SingleLorentzianFixedLinewidthModel.h"


// SingleLorentzianFixedLinewidthModel::SingleLorentzianFixedLinewidthModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      linewidth:              the linewidth value to be used for the Lorentzian profile
//      backgroundModel:        an object of class BackgroundModel containing a
//                              model for the background of the star.
//

SingleLorentzianFixedLinewidthModel::SingleLorentzianFixedLinewidthModel(RefArrayXd const covariates, const double linewidth,
BackgroundModel &backgroundModel)
: Model(covariates),
  linewidth(linewidth)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// SingleLorentzianFixedLinewidthModel::SingleLorentzianFixedLinewidthModel()
//
// PURPOSE: 
//      Destructor.
//

SingleLorentzianFixedLinewidthModel::~SingleLorentzianFixedLinewidthModel()
{

}










// SingleLorentzianFixedLinewidthModel::predict()
//
// PURPOSE:
//      Builds the predictions from a SingleLorentzianFixedLinewidth model to
//      construct a multi-modal posterior distribution (islands model).
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
//      (i) Mode central frequency
//      (ii) Mode profile height 
//

void SingleLorentzianFixedLinewidthModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());


    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);

    Functions::modeProfile(singleModePrediction, covariates, centralFrequency, height, linewidth);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
