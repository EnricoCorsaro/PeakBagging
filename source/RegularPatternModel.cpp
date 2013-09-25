#include "RegularPatternModel.h"


// RegularPatternModel::RegularPatternModel()
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

RegularPatternModel::RegularPatternModel(const RefArrayXd covariates, const vector<int> &NparametersPerType, double nuMax)
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










// RegularPatternModel::RegularPatternModel()
//
// PURPOSE: 
//      Destructor.
//

RegularPatternModel::~RegularPatternModel()
{

}










// RegularPatternModel::predict()
//
// PURPOSE:
//      Builds the predictions from a RegularPattern model based on the
//      Tassoul's asymptotic relation for p modes and on the regular patterns
//      and formalism proposed by Mosser et al. 2012, A&A, 540, A143 for mixed modes 
//      and Mosser et al. 2012, A&A, 548, A10 for rotational splittings.
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
//      (1) Reference Frequency of central radial order
//      (2) White noise background (flat noise level)
//      (3) Large frequency separation (DeltaNu)
//      (4) Small frequency separation 02 (deltaNu02)
//      (5) Small frequency separation 01 (deltaNu01)
//      (5) Mode profile heights
//      (6) Mode profile linewidths

void RegularPatternModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    Nparameters = modelParameters.size();
    ArrayXd pressureModePrediction = ArrayXd::Zero(covariates.size());
    //double linewidth = 0.15;
    cout << modelParameters << endl;

    for (int radialOrder = 0; radialOrder < NradialOrders; ++radialOrder)
    {
        for (int angularDegree = 0; angularDegree < NangularDegrees; ++angularDegree)
        {
            double frequencyOfPmode = computeSinglePmodeFrequency(relativeRadialOrders[radialOrder], modelParameters(0),
                                                                  modelParameters(NglobalParameters + radialOrder), 
                                                                  modelParameters(NglobalParameters + NradialOrders + radialOrder),
                                                                  modelParameters(NglobalParameters + 2*NradialOrders + radialOrder),
                                                                  angularDegree);
            
            if ((frequencyOfPmode > covariates.minCoeff()) && (frequencyOfPmode < covariates.maxCoeff()))
            {
                Functions::modeProfile(pressureModePrediction, covariates, frequencyOfPmode, 
                                       modelParameters(NglobalParameters + 3*NradialOrders + radialOrder + angularDegree), 
                                       modelParameters(NglobalParameters + 3*NradialOrders + NpressureModes + 
                                       radialOrder + angularDegree));
                                       //linewidth);
                predictions += pressureModePrediction;
            }
            else
            {
                continue;
            }
        }
    }
    
    predictions += modelParameters(1);           // Add flat noise level component

}










// RegularPatternModel::computeSinglePmodeFrequency()
//
// PURPOSE:
//      Computes the frequency of a single pressure mode 
//      based on the Tassoul's asymptotic relation. The formalism
//      adopted is that used by Mosser et al. 2012, A&A, 540, A143 and
//      Mosser et al. 2012, A&A, 548, A10
//
// INPUT:
//      relativeRadialOrder:            an integer specifying the radial
//                                      order to be considered, expressed as
//                                      the difference with respect to the
//                                      radial order of the central radial mode
//      referenceFrequency:             the frequency corresponding to the central radial mode
//      DeltaNu:                        the large frequency separation of p modes
//      deltaNu02:                      the small frequency separation 02 of p modes
//      deltaNu01:                      the small frequency separation 01 of p modes
//      epsilon:                        the phase shift constant
//      angularDegree:                  the angular degree of the mode
//
// OUTPUT:
//      A double containing the frequency of the desidered mode.
//

double RegularPatternModel::computeSinglePmodeFrequency(int relativeRadialOrder, double referenceFrequency, double DeltaNu, 
                                                        double deltaNu02, double deltaNu01, int angularDegree)
{
    // Compute the alpha factor used for the curvature correction
    // double alphaCorrection = 0.015 * pow(DeltaNu,-0.32);        // Mosser et al. 2012, A&A, 548, A10
   

    // Compute the first order term of the Tassoul's asymptotic relation in units of DeltaNu
    
    double linearTerm = 0.0;
    double startingFrequency = nuMax + referenceFrequency;
    
    switch (angularDegree)
    {
        case 0:        
        linearTerm = DeltaNu*relativeRadialOrder + startingFrequency;
        break;
    
        case 1:
        linearTerm = DeltaNu*(relativeRadialOrder + 0.5) + startingFrequency - deltaNu01;
        break;

        case 2:
        linearTerm = DeltaNu*(relativeRadialOrder) + startingFrequency - deltaNu02;
        break;

        case 3:
        linearTerm = DeltaNu*(relativeRadialOrder + 0.5) + startingFrequency - 0.280*DeltaNu;
        break;
    }


    // Compute the second order term of the Tassoul's asymptotic relation in units of DeltaNu
    // double quadraticTerm = alphaCorrection/2.0 * relativeRadialOrder * relativeRadialOrder;
    // return DeltaNu * (linearTerm + quadraticTerm);

    return linearTerm;
}










// RegularPatternModel::predictDeltaNuFromNuMax()
//
// PURPOSE:
//      Predict the average DeltaNu based on nuMax by means of
//      the scaling relation revised by Huber et al. (2011), ApJ 743, 143
//      for the subset of Kepler Red Giants.
//      
// INPUT:
//
// OUTPUT:
//      A double containing the value of DeltaNu 
//

double RegularPatternModel::predictDeltaNuFromNuMax()
{
 double alpha = 0.267;
 double beta = 0.760;
 
 return alpha*pow(nuMax,beta);
}











// RegularPatternModel::predictDeltaNu01FromDeltaNu()
//
// PURPOSE:
//      Predict the average deltaNu01 based on DeltaNu by means of
//      the scaling relation revised by Mosser et al. (2011), A&A 525, 9
//      
// INPUT:
//      DeltaNu:    a double specifying the DeltaNu to be used
//
// OUTPUT:
//      A double containing the value of deltaNu01 
//

double RegularPatternModel::predictDeltaNu01FromDeltaNu(double DeltaNu)
{
 double A = -0.056;
 double B = -0.002;
 
 return DeltaNu*(A + B*log10(DeltaNu));
}










// RegularPatternModel::predictDeltaNu02FromDeltaNu()
//
// PURPOSE:
//      Predict the average deltaNu02 based on DeltaNu by means of
//      the scaling relation revised by Mosser et al. (2011), A&A 525, 9
//      
// INPUT:
//      DeltaNu:    a double specifying the DeltaNu to be used
//
// OUTPUT:
//      A double containing the value of deltaNu02 
//

double RegularPatternModel::predictDeltaNu02FromDeltaNu(double DeltaNu)
{
 double A = 0.131;
 double B = -0.033;
 
 return DeltaNu*(A + B*log10(DeltaNu));
}










// RegularPatternModel::predictDeltaNu03FromDeltaNu()
//
// PURPOSE:
//      Predict the average deltaNu03 based on DeltaNu by means of
//      the scaling relation revised by Mosser et al. (2011), A&A 525, 9
//      
// INPUT:
//      DeltaNu:    a double specifying the DeltaNu to be used
//
// OUTPUT:
//      A double containing the value of deltaNu03 
//

double RegularPatternModel::predictDeltaNu03FromDeltaNu(double DeltaNu)
{
 double A = 0.280;
 
 return DeltaNu*A;
}
