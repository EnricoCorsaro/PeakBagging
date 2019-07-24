#include "DoubleLorentzianFixedLinewidthModel.h"


// DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel()
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

DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel(RefArrayXd const covariates, const double linewidth,
BackgroundModel &backgroundModel)
: Model(covariates),
  linewidth(linewidth)
{
    backgroundPrediction.resize(covariates.size());
    backgroundModel.predict(backgroundPrediction);
    responseFunction = backgroundModel.getResponseFunction();
}










// DoubleLorentzianFixedLinewidthModel::DoubleLorentzianFixedLinewidthModel()
//
// PURPOSE: 
//      Destructor.
//

DoubleLorentzianFixedLinewidthModel::~DoubleLorentzianFixedLinewidthModel()
{

}










// DoubleLorentzianFixedLinewidthModel::predict()
//
// PURPOSE:
//      Builds the predictions from a DoubleLorentzianFixedLinewidth model to
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
//      (iii) Small frequency separation d02
//

void DoubleLorentzianFixedLinewidthModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    ArrayXd singleModePrediction = ArrayXd::Zero(covariates.size());


    // Initialize parameters of the oscillation mode

    double centralFrequency = modelParameters(0);
    double height = modelParameters(1);
    double smallSeparation = modelParameters(2);


    // Add firt profile with reference frequency given by the centralFrequency

    Functions::modeProfile(singleModePrediction, covariates, centralFrequency, height, linewidth);
    predictions += singleModePrediction;
    

    // Add second profile with reference frequency given by the centralFrequency - smallSeparation
    
    Functions::modeProfile(singleModePrediction, covariates, centralFrequency - smallSeparation, height, linewidth);
    predictions += singleModePrediction;

    
    // Adjust by apodization of the signal
    
    predictions *= responseFunction;


    // Add background component

    predictions += backgroundPrediction;
}
