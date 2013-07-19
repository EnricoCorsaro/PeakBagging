#include "RegularPatternModel.h"


// RegularPatternModel::RegularPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:     one-dimensional array containing the values
//                      of the independent variable.
//      NradialOrders:        the number of radial orders to be fitted
//

RegularPatternModel::RegularPatternModel(const RefArrayXd covariates, const int NradialOrders)
: Model(covariates),
  NradialOrders(NradialOrders)
{

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
//                          contains the value of one free parameter of the model
//
// OUTPUT:
//      void
//
// NOTE:
//      The free parameters are to be given in the order
//      (1) Frequency of central radial order
//      (2) Large frequency separation (DeltaNu)
//      (3) Small frequency separation 02 (deltaNu02)
//      (4) ...

void RegularPatternModel::predict(RefArrayXd predictions, const RefArrayXd modelParameters)
{
    Nparameters = modelParameters.size();

    double frequencyOfCentralRadialMode = modelParameters(0);
    double DeltaNu = modelParameters(1);
    double deltaNu02 = modelParameters(2);
    
    vector<int> relativeRadialOrders(NradialOrders);
    ArrayXd modes = ArrayXd::Zero(predictions.size());
    
    double height = 16.0;
    double linewidth = 0.4;

    // radialOrderOfCentralRadialMode = floor(frequencyOfCentralRadialMode / DeltaNu);
    double epsilon = frequencyOfCentralRadialMode / DeltaNu - fmod(frequencyOfCentralRadialMode, DeltaNu);
    
    for (int i=0; i < relativeRadialOrders.size(); ++i)
    {
        relativeRadialOrders[i] = -static_cast<int>(floor(NradialOrders/2.0)) + i;

        for (int j=0; j < 4; ++j)
        {
            double frequencyOfPmode = frequencyOfCentralRadialMode + relativeRadialOrders[i]*DeltaNu;
            // computeSinglePmodeFrequency(relativeRadialOrders[i], frequencyOfCentralRadialMode,
                                       //  DeltaNu, deltaNu02, epsilon, j);
            if (frequencyOfPmode >= covariates.minCoeff() && frequencyOfPmode <= covariates.maxCoeff())
            {
                Functions::modeProfile(modes, covariates, frequencyOfPmode, height, linewidth);
                predictions += modes;
            }
            else
                continue;
        }
    }

    predictions += 0.01;
    
    //computePmodeFrequencies(relativeRadialOrders, frequencyOfCentralRadialMode, DeltaNu, deltaNu02);


    /*
    for (int i = 0; i < frequenciesOfPmodes.size(); i++)
    {
        Functions::modeProfile(modes, covariates, frequenciesOfPmodes[i], height, linewidth);
        predictions += modes;
    }
    */

}







// RegularPatternModel::getNparameters()
//
// PURPOSE: 
//      Get the private data member Nparameters;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int RegularPatternModel::getNparameters()
{
    return Nparameters;
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
//      frequencyOfCentralRadialMode:   the frequency corresponding to the central radial mode
//      DeltaNu:                        the large frequency separation of p modes
//      deltaNu02:                      the small frequency separation 02 of p modes
//      epsilon:                        the phase shift constant
//      angularDegree:                  the angular degree of the mode
//
// OUTPUT:
//      A double containing the frequency of the desidered mode.
//

double RegularPatternModel::computeSinglePmodeFrequency(int relativeRadialOrder, double frequencyOfCentralRadialMode, 
                                                        double DeltaNu, double deltaNu02, double epsilon, int angularDegree)
{
    
    // Compute the alpha factor used for the curvature correction
    
    double alphaCorrection = 0.015 * pow(DeltaNu,-0.32);        // Mosser et al. 2012, A&A, 548, A10
   

    // Compute the constant d0l to be used in the asymptotic relation.
    // The values are derived from the univeral pattern proposed by Mosser et al. 2011, A&A, 525, L9

    double constantDell = 0.0;

    switch (angularDegree)
    {
        case 0:
        constantDell = 0.0;
        break;
        
        case 1:
        constantDell = -0.056 + (-0.002)*log(DeltaNu);
        break;

        case 2:
        constantDell = deltaNu02/DeltaNu;
        break;

        case 3:
        constantDell = 0.280;
        break;
    }

    // Compute the first order term of the Tassoul's asymptotic relation in units of DeltaNu

    double linearTerm = relativeRadialOrder + radialOrderOfCentralRadialMode + epsilon + angularDegree/2.0 - constantDell;
    
    
    // Compute the second order term of the Tassoul's asymptotic relation in units of DeltaNu
    
    double quadraticTerm = alphaCorrection/2.0 * relativeRadialOrder * relativeRadialOrder;

    return DeltaNu * (linearTerm + quadraticTerm);
}












// RegularPatternModel::computePmodeFrequencies()
//
// PURPOSE:
//      Computes the frequency of an entire set of p modes 
//      based on the Tassoul's asymptotic relation. The formalism
//      adopted is that used by Mosser et al. 2012, A&A, 540, A143 and
//      Mosser et al. 2012, A&A, 548, A10
//
// INPUT:
//      relativeRadialOrders:           a vector of integers specifying the radial
//                                      orders to be considered, expressed as
//                                      the difference with respect to the
//                                      radial order of the central radial mode.
//      frequencyOfCentralRadialMode:   the frequency corresponding to the central radial mode
//      DeltaNu:                        the large frequency separation of p modes
//      deltaNu02:                      the small frequency separation 02 of p modes
//
// OUTPUT:
//      The output frequencies, angular degress and radial orders are stored into vectors.
//

void RegularPatternModel::computePmodeFrequencies(vector<int> relativeRadialOrders, double frequencyOfCentralRadialMode, 
                                                         double DeltaNu, double deltaNu02)
{

    double epsilon = frequencyOfCentralRadialMode / DeltaNu - fmod(frequencyOfCentralRadialMode, DeltaNu);
    radialOrderOfCentralRadialMode = floor(frequencyOfCentralRadialMode / DeltaNu);

    for (int i = 0; i < relativeRadialOrders.size(); ++i)
    {
        for (int angularDegree = 0; angularDegree < 4; ++angularDegree)
        {
            double frequencyOfPmode = frequencyOfCentralRadialMode + DeltaNu/2.0 - deltaNu02; 
            // computeSinglePmodeFrequency(relativeRadialOrders[i], frequencyOfCentralRadialMode,
            //                                                      DeltaNu, deltaNu02, epsilon, angularDegree);
            frequenciesOfPmodes.push_back(frequencyOfPmode);
            angularDegreesOfPmodes.push_back(angularDegree);
            radialOrdersOfPmodes.push_back(relativeRadialOrders[i]+radialOrderOfCentralRadialMode);
        }
    }
}

