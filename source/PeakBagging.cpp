// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "PeakBagging.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "Functions.h"
#include "File.h"
#include "MultiEllipsoidSampler.h"
#include "KmeansClusterer.h"
#include "EuclideanMetric.h"
#include "Prior.h"
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "ExponentialLikelihood.h"
#include "LorentzianMixtureModel.h"
#include "CanBackgroundModel.h"
#include "FerozReducer.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
  

    // Check number of arguments for main function
    
    if (argc != 3)
    {
        cerr << "Usage: peakbagging <input file> <output directory>" << endl;
        exit(EXIT_FAILURE);
    }


    // Read data from input file specified

    ifstream inputFile(argv[1]);
    if (!inputFile.good())
    {
        cerr << "Error opening input file" << endl;
        exit(EXIT_FAILURE);
    }

    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

   
    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
   
   
    // --------------------------------------------------
    // ----- Zeroth step. Set up problem dimensions -----
    // --------------------------------------------------

    vector<int> NparametersPerType(2);                              // Vector containing the configuring parameters for the LorentzianMixture model
    int NprofileParameters = 3;
    NparametersPerType[0] = NprofileParameters;                     // Number of parameters determining the shape 
                                                                    // of the mode profile (central frequency, height, linewidth, inclination angle, etc.)
    int Nmodes = 2; 
    NparametersPerType[1] = Nmodes;                                 // Total number of modes to fit
    int Ndimensions = NprofileParameters*Nmodes;
    cout << "--------------------------------------" << endl;
    cout << "Inference problem has " << Ndimensions << " dimensions." << endl;
    cout << "--------------------------------------" << endl;
    
    
    // -------------------------------------------------------------------
    // ----- First step. Set up the models for the inference problem ----- 
    // -------------------------------------------------------------------
    
    // Chose a model for background and configure it, then do the same for the peakbagging model.
    
    CanBackgroundModel backgroundModel(covariates);
    string outputDirName(argv[2]);
    string backgroundConfiguringParameters = outputDirName + "CanModel_configuringParameters.txt";
    backgroundModel.readConfiguringParametersFromFile(backgroundConfiguringParameters);
    LorentzianMixtureModel model(covariates, NparametersPerType, backgroundModel);


    // -------------------------------------------------------
    // ----- Second step. Set up all prior distributions -----
    // -------------------------------------------------------
    
    // Total number of prior types

    vector<Prior*> ptrPriors(1);


    // Uniform Prior

    ArrayXd parametersMinima(Ndimensions);                      // Minima
    ArrayXd parametersMaxima(Ndimensions);                      // Maxima
    
    ArrayXd centralFrequencyMinima(Nmodes);                     // Central frequency of oscillation mode
    ArrayXd naturalLogarithmOfHeightsMinima(Nmodes);            // Natural logarithm of mode height
    ArrayXd linewidthsMinima(Nmodes);                           // Mode Linewidth (FWHM)
    ArrayXd centralFrequencyMaxima(Nmodes);                     // Central frequency of oscillation mode
    ArrayXd naturalLogarithmOfHeightsMaxima(Nmodes);            // Natural logarithm of mode height
    ArrayXd linewidthsMaxima(Nmodes);                           // Mode Linewidth (FWHM)

    
    centralFrequencyMinima(0) = 119.5;
    centralFrequencyMaxima(0) = 120.0;
    naturalLogarithmOfHeightsMinima(0) = 5.5;
    naturalLogarithmOfHeightsMaxima(0) = 6.2;
    linewidthsMinima(0) = 0.2;
    linewidthsMaxima(0) = 0.4;
    

    centralFrequencyMinima(1) = 120.4;
    centralFrequencyMaxima(1) = 122.5;
    naturalLogarithmOfHeightsMinima(1) = 5.8;
    naturalLogarithmOfHeightsMaxima(1) = 6.9;
    linewidthsMinima(1) = 0.28;
    linewidthsMaxima(1) = 0.46;
/*    
    centralFrequencyMinima(2) = 126.8;
    centralFrequencyMaxima(2) = 127.3;
    naturalLogarithmOfHeightsMinima(2) = 3.6;
    naturalLogarithmOfHeightsMaxima(2) = 6.0;
    linewidthsMinima(2) = 0.02;
    linewidthsMaxima(2) = 0.3;
    
    centralFrequencyMinima(3) = 130.3;
    centralFrequencyMaxima(3) = 130.9;
    naturalLogarithmOfHeightsMinima(3) = 3.4;
    naturalLogarithmOfHeightsMaxima(3) = 5.2;
    linewidthsMinima(3) = 0.05;
    linewidthsMaxima(3) = 0.20;

    centralFrequencyMinima(4) = 131.7;
    centralFrequencyMaxima(4) = 132.8;
    naturalLogarithmOfHeightsMinima(4) = 5.6;
    naturalLogarithmOfHeightsMaxima(4) = 6.9;
    linewidthsMinima(4) = 0.2;
    linewidthsMaxima(4) = 0.45;

    centralFrequencyMinima(5) = 137.0;
    centralFrequencyMaxima(5) = 138.5;
    naturalLogarithmOfHeightsMinima(5) = 3.6;
    naturalLogarithmOfHeightsMaxima(5) = 6.0;
    linewidthsMinima(5) = 0.05;
    linewidthsMaxima(5) = 0.3;

    centralFrequencyMinima(0) = 142.5;
    centralFrequencyMaxima(0) = 144.0;
    naturalLogarithmOfHeightsMinima.fill(4.6);
    naturalLogarithmOfHeightsMaxima.fill(8.01);
    linewidthsMinima.fill(0.01);
    linewidthsMaxima.fill(0.5);

*/

    parametersMinima.segment(0, Nmodes) = centralFrequencyMinima;
    parametersMinima.segment(Nmodes, Nmodes) = naturalLogarithmOfHeightsMinima;
    parametersMinima.segment(2*Nmodes, Nmodes) = linewidthsMinima;
    parametersMaxima.segment(0, Nmodes) = centralFrequencyMaxima;
    parametersMaxima.segment(Nmodes, Nmodes) = naturalLogarithmOfHeightsMaxima;
    parametersMaxima.segment(2*Nmodes, Nmodes) = linewidthsMaxima;

    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    

    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    ExponentialLikelihood likelihood(observations, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the K-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 6;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 1000;                     // Maximum (and initial) number of live points evolving within the nested sampling process. 
    int minNobjects = 1000;                          // Minimum number of live points allowed in the computation
    int maxNdrawAttempts = 10000;                   // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = static_cast<int>(initialNobjects*0.5);    // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = static_cast<int>(initialNobjects*0.05);        // Clustering is only happening every N iterations.
    double initialEnlargementFraction = 1.20;       // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.8;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.05;                // Termination factor for nested sampling process.

    
    // Save configuring parameters into an ASCII file

    ofstream outputFile; 
    string fullPath = outputDirName + "peakBagging_configuringParameters.txt";
    File::openOutputFile(outputFile, fullPath);
    File::configuringParametersToFile(outputFile, initialNobjects, minNobjects, minNclusters, maxNclusters, NinitialIterationsWithoutClustering,
                                     NiterationsWithSameClustering, maxNdrawAttempts, initialEnlargementFraction, shrinkingRate, terminationFactor);
    outputFile.close();

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);


    // Choose which reducer of the live points to adopt and include it in the run.

    double toleranceOnEvidence = 0.01;
    FerozReducer livePointsReducer(nestedSampler, toleranceOnEvidence);
    nestedSampler.run(livePointsReducer, terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


    // -------------------------------------------------------
    // ----- Last step. Save the results in output files -----
    // -------------------------------------------------------
    
    Results results(nestedSampler);
    results.writeParametersToFile(outputDirName + "parameter");
    results.writeLogLikelihoodToFile(outputDirName + "likelihoodDistribution.txt");
    results.writeEvidenceInformationToFile(outputDirName + "evidenceInformation.txt");
    results.writePosteriorProbabilityToFile(outputDirName + "posteriorDistribution.txt");
    results.writeParametersSummaryToFile(outputDirName + "parameterSummary.txt");

    
    // END

    return EXIT_SUCCESS;
}
