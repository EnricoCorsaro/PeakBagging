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
#include "FerozReducer.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
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
   
   
    // Zero step - Set up problem dimensions

    vector<int> NparametersPerType(5);                              // Vector containing the configuring parameters for the LorentzianMixture model
    int NglobalParameters = 1;
    NparametersPerType[0] = NglobalParameters;                      // Total number of global parameters to be fitted 
                                                                    // (e.g. referenceFrequency for central radial order, flat noise level, 
                                                                    // period spacing, etc.)
    int NprofileParameters = 3;
    NparametersPerType[1] = NprofileParameters;                     // Number of parameters determining the shape 
                                                                    // of the mode profile (central frequency, height, linewidth, inclination angle, etc.)
    int Nmodes = 1; 
    NparametersPerType[2] = Nmodes;                                 // Total number of modes to fit
    
    
    // First step - Set up a model for the inference problem
    
    LorentzianMixtureModel model(covariates, NparametersPerType);


    // Second step - Setting Prior distribution
    
    vector<Prior*> ptrPriors(1);
    

    // Total number of free parameters (dimensions) of the problem

    int Ndimensions = NglobalParameters + NprofileParameters*Nmodes;


    // Uniform Prior

    ArrayXd parametersMinima(Ndimensions);                      // Minima
    
    double noiseLevelMinimum = 5;                               // Flat noise level
    parametersMinima(0) = noiseLevelMinimum;

    ArrayXd centralFrequencyMinima(Nmodes);                     // Central frequency of oscillation mode
    centralFrequencyMinima(0) = 200.0;
    
    ArrayXd naturalLogarithmOfHeightsMinima(Nmodes);            // Natural logarithm of mode height
    naturalLogarithmOfHeightsMinima.fill(4.6);
    
    ArrayXd linewidthsMinima(Nmodes);                           // Mode Linewidth (FWHM)
    linewidthsMinima.fill(0.05);
    
    parametersMinima.segment(NglobalParameters, Nmodes) = centralFrequencyMinima;
    parametersMinima.segment(NglobalParameters + Nmodes, Nmodes) = naturalLogarithmOfHeightsMinima;
    parametersMinima.segment(NglobalParameters + 2*Nmodes, Nmodes) = linewidthsMinima;
    

    ArrayXd parametersMaxima(Ndimensions);                      // Maxima
    
    double noiseLevelMaximum = 20;                              // Flat noise level
    parametersMaxima(0) = noiseLevelMaximum;

    ArrayXd centralFrequencyMaxima(Nmodes);                     // Central frequency of oscillation mode
    centralFrequencyMaxima(0) = 201.0;
    
    ArrayXd naturalLogarithmOfHeightsMaxima(Nmodes);            // Natural logarithm of mode height
    naturalLogarithmOfHeightsMaxima.fill(8.01);
    
    ArrayXd linewidthsMaxima(Nmodes);                           // Mode Linewidth (FWHM)
    linewidthsMaxima.fill(0.5);
    
    parametersMaxima.segment(NglobalParameters, Nmodes) = centralFrequencyMaxima;
    parametersMaxima.segment(NglobalParameters + Nmodes, Nmodes) = naturalLogarithmOfHeightsMaxima;
    parametersMaxima.segment(NglobalParameters + 2*Nmodes, Nmodes) = linewidthsMaxima;

    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    //*/


    // Normal Prior
    /*

    ArrayXd parametersMean(Ndimensions);                    // Mean
    
    double noiseLevelMean = 5;                              // Flat noise level
    parametersMean(0) = noiseLevelMean;

    ArrayXd centralFrequencyMean(Nmodes);                   // Central frequency of oscillation mode
    centralFrequencyMean(0) = 200.0;
    
    ArrayXd naturalLogarithmOfHeightsMean(Nmodes);          // Natural logarithm of mode height
    naturalLogarithmOfHeightsMean.fill(4.6);
    
    ArrayXd linewidthsMean(Nmodes);                         // Mode Linewidth (FWHM)
    linewidthsMean.fill(0.05);
    
    parametersMean.segment(NglobalParameters, Nmodes) = centralFrequencyMean;
    parametersMean.segment(NglobalParameters + Nmodes, Nmodes) = naturalLogarithmOfHeightsMean;
    parametersMean.segment(NglobalParameters + 2*Nmodes, Nmodes) = linewidthsMean;
    

    ArrayXd parametersSDV(Ndimensions);                     // SDV
    
    double noiseLevelSDV = 20;                              // Flat noise level
    parametersSDV(0) = noiseLevelSDV;

    ArrayXd centralFrequencySDV(Nmodes);                    // Central frequency of oscillation mode
    centralFrequencySDV(0) = 201.0;
    
    ArrayXd naturalLogarithmOfHeightsSDV(Nmodes);           // Natural logarithm of mode height
    naturalLogarithmOfHeightsSDV.fill(8.01);
    
    ArrayXd linewidthsSDV(Nmodes);                          // Mode Linewidth (FWHM)
    linewidthsSDV.fill(0.5);
    
    parametersSDV.segment(NglobalParameters, Nmodes) = centralFrequencySDV;
    parametersSDV.segment(NglobalParameters + Nmodes, Nmodes) = naturalLogarithmOfHeightsSDV;
    parametersSDV.segment(NglobalParameters + 2*Nmodes, Nmodes) = linewidthsSDV;

    NormalPrior normalPrior(parametersMean, parametersSDV);
    ptrPriors[0] = &normalPrior;

    
    // Super-Gaussian Prior
    ArrayXd parametersWOP(Ndimensions);                     // WOP
    
    double noiseLevelWOP = 20;                              // Flat noise level
    parametersWOP(0) = noiseLevelWOP;

    ArrayXd centralFrequencyWOP(Nmodes);                    // Central frequency of oscillation mode
    centralFrequencyWOP(0) = 201.0;
    
    ArrayXd naturalLogarithmOfHeightsWOP(Nmodes);           // Natural logarithm of mode height
    naturalLogarithmOfHeightsWOP.fill(8.01);
    
    ArrayXd linewidthsWOP(Nmodes);                          // Mode Linewidth (FWHM)
    linewidthsWOP.fill(0.5);
    
    parametersWOP.segment(NglobalParameters, Nmodes) = centralFrequencyWOP;
    parametersWOP.segment(NglobalParameters + Nmodes, Nmodes) = naturalLogarithmOfHeightsWOP;
    parametersWOP.segment(NglobalParameters + 2*Nmodes, Nmodes) = linewidthsWOP;
    
    SuperGaussianPrior superGaussianPrior(parametersMean, parametersSDV, parametersWOP);
    ptrPriors[0] = &superGaussianPrior;
*/
    

    // Third step - Set up the likelihood function to be used
    
    ExponentialLikelihood likelihood(observations, model);
    

    // Fourth step - Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 6;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Start nested sampling process
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 10000;                    // Maximum (and initial) number of live points evolving within the nested sampling process. 
    int minNobjects = 1000;                         // Minimum number of live points allowed in the computation
    int maxNdrawAttempts = 10000;                   // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = static_cast<int>(initialNobjects*10*0.05);    // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = static_cast<int>(initialNobjects*10*0.005);        // Clustering is only happening every X iterations.
    double initialEnlargementFraction = 1.50;       // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.8;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.05;                // Termination factor for nested sampling process.

    
    // Save configuring parameters into an ASCII file

    ofstream outputFile;
    string fullPath = "peakBagging_configuringParameters.txt";
    File::openOutputFile(outputFile, fullPath);
    File::configuringParametersToFile(outputFile, initialNobjects, minNobjects, minNclusters, maxNclusters, NinitialIterationsWithoutClustering,
                                     NiterationsWithSameClustering, maxNdrawAttempts, initialEnlargementFraction, shrinkingRate, terminationFactor);
    outputFile.close();


    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);

    double toleranceOnEvidence = 0.01;
    FerozReducer livePointsReducer(nestedSampler, toleranceOnEvidence);

    nestedSampler.run(livePointsReducer, terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


    // Save the results in output files
    
    Results results(nestedSampler);
    string outputDirName(argv[2]);
    results.writeParametersToFile(outputDirName + "parameter");
    results.writeLogLikelihoodToFile(outputDirName + "likelihoodDistribution.txt");
    results.writeEvidenceInformationToFile(outputDirName + "evidenceInformation.txt");
    results.writePosteriorProbabilityToFile(outputDirName + "posteriorDistribution.txt");
    results.writeParametersSummaryToFile(outputDirName + "parameterSummary.txt");

    return EXIT_SUCCESS;
}
