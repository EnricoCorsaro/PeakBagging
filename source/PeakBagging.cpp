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
#include "RegularPatternModel.h"
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

    vector<int> NparametersPerType(5);                              // Vector containing the configuring parameters for the RegularPattern model
    int NglobalParameters = 2;
    NparametersPerType[0] = NglobalParameters;                      // Total number of global parameters to be fitted 
                                                                    // (e.g. referenceFrequency for central radial order, flat noise level, 
                                                                    // period spacing, etc.)
    int NprofileParameters = 2;
    NparametersPerType[1] = NprofileParameters;                     // Number of parameters determining the shape 
                                                                    // of the mode profile (height, linewidth, inclination angle, etc.)
                                                                    // N.B the mode frequency may not be considered as a profile parameter
    int NradialOrders = 1; 
    NparametersPerType[2] = NradialOrders;                          // Number of radial orders expected in the PSD
    
    int NangularDegrees = 3;
    NparametersPerType[3] = NangularDegrees;                        // Maximum number of angular degrees to be fitted
                                                                    // N.B this number can differ from the product NangularDegrees*NradialOrders
    int NmixedModes = 0;
    NparametersPerType[4] = 0;                                      // Total number of mixed modes to be fitted
  
    double nuMax = 133.0 ;      // microHz
    
    
    // First step - Set up a model for the inference problem
    
    RegularPatternModel model(covariates, NparametersPerType, nuMax);
    double referenceDeltaNu = model.predictDeltaNuFromNuMax();

    // Second step - Setting Prior distribution
    
    vector<Prior*> ptrPriors(1);
    

    // Total number of free parameters (dimensions) of the problem

    int Ndimensions = NglobalParameters + 3*NradialOrders + NprofileParameters*NradialOrders*NangularDegrees;


    // Uniform Prior

    ArrayXd parametersMinima(Ndimensions);                              // Minima
    
    double referenceFrequencyMinimum = -referenceDeltaNu/2. - referenceDeltaNu*0.1;            // Frequency of central radial mode
    parametersMinima(0) = referenceFrequencyMinimum;
    
    double noiseLevelMinimum = 10;                                      // Flat noise level
    parametersMinima(1) = noiseLevelMinimum;

    ArrayXd DeltaNuMinima(NradialOrders);                               // Large frequency spacing DeltaNu
    DeltaNuMinima.fill(10.0);
    
    ArrayXd deltaNu02Minima(NradialOrders);                             // Small frequency spacing deltaNu02
    deltaNu02Minima.fill(0.6);
    
    ArrayXd deltaNu01Minima(NradialOrders);                             // Small frequency spacing deltaNu01
    deltaNu01Minima.fill(-0.6);
    
    ArrayXd naturalLogarithmOfHeightsMinima(NradialOrders*NangularDegrees);              // Natural logarithm of Heights
    naturalLogarithmOfHeightsMinima.fill(4.6);
    
    ArrayXd linewidthsMinima(NradialOrders*NangularDegrees);            // Linewidths
    linewidthsMinima.fill(0.05);
    
    parametersMinima.segment(NglobalParameters, NradialOrders) = DeltaNuMinima;
    parametersMinima.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02Minima;
    parametersMinima.segment(NglobalParameters + 2*NradialOrders, NradialOrders) = deltaNu01Minima;
    parametersMinima.segment(NglobalParameters + 3*NradialOrders, NradialOrders*NangularDegrees) = naturalLogarithmOfHeightsMinima;
    parametersMinima.segment(NglobalParameters + 3*NradialOrders + NradialOrders*NangularDegrees, NradialOrders*NangularDegrees) = linewidthsMinima;
    
    ArrayXd parametersMaxima(Ndimensions);                              // Maxima
    
    double referenceFrequencyMaximum = referenceDeltaNu/2.;
    parametersMaxima(0) = referenceFrequencyMaximum;
    
    double noiseLevelMaximum = 30;
    parametersMaxima(1) = noiseLevelMaximum; 
    
    ArrayXd DeltaNuMaxima(NradialOrders);
    DeltaNuMaxima.fill(12.0);
    
    ArrayXd deltaNu02Maxima(NradialOrders);
    deltaNu02Maxima.fill(1.8);
    
    ArrayXd deltaNu01Maxima(NradialOrders);
    deltaNu01Maxima.fill(0.1);
    
    ArrayXd naturalLogarithmOfHeightsMaxima(NradialOrders*NangularDegrees);
    naturalLogarithmOfHeightsMaxima.fill(8.01);
    
    ArrayXd linewidthsMaxima(NradialOrders*NangularDegrees);
    linewidthsMaxima.fill(0.5);
    
    parametersMaxima.segment(NglobalParameters, NradialOrders) = DeltaNuMaxima;
    parametersMaxima.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02Maxima;
    parametersMaxima.segment(NglobalParameters + 2*NradialOrders, NradialOrders) = deltaNu01Maxima;
    parametersMaxima.segment(NglobalParameters + 3*NradialOrders, NradialOrders*NangularDegrees) = naturalLogarithmOfHeightsMaxima;
    parametersMaxima.segment(NglobalParameters + 3*NradialOrders + NradialOrders*NangularDegrees, NradialOrders*NangularDegrees) = linewidthsMaxima;
    
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    //*/


    // Normal Prior
    /*
    ArrayXd parametersMean(Ndimensions);
    ArrayXd parametersSDV(Ndimensions);
    ArrayXd DeltaNuMean = ArrayXd::Zero(NradialOrders);
    ArrayXd DeltaNuSDV = ArrayXd::Zero(NradialOrders);
    ArrayXd deltaNu02Mean = ArrayXd::Zero(NradialOrders);
    ArrayXd deltaNu02SDV = ArrayXd::Zero(NradialOrders);
    ArrayXd heightsMean = ArrayXd::Zero(NradialOrders*NangularDegrees);
    ArrayXd heightsSDV = ArrayXd::Zero(NradialOrders*NangularDegrees);
    double nuMaxMean = 26.0;
    double nuMaxSDV = 30.0;
    double noiseLevelMean = 0.1;
    double noiseLevelSDV = 20;
    
    parametersMean(0) = nuMaxMean;
    parametersMean(1) = noiseLevelMean; 
    DeltaNuMean += 2.5;
    deltaNu02Mean += 0.1;
    heightsMean += 100.0;        
    parametersMean.segment(NglobalParameters, NradialOrders) = DeltaNuMean;
    parametersMean.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02Mean;
    parametersMean.segment(NglobalParameters + 2*NradialOrders, NradialOrders*NangularDegrees) = heightsMean;
    
    parametersSDV(0) = nuMaxSDV;
    parametersSDV(1) = noiseLevelSDV; 
    DeltaNuSDV += 2.5;
    deltaNu02SDV += 0.1;
    heightsSDV += 100.0;        
    parametersSDV.segment(NglobalParameters, NradialOrders) = DeltaNuSDV;
    parametersSDV.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02SDV;
    parametersSDV.segment(NglobalParameters + 2*NradialOrders, NradialOrders*NangularDegrees) = heightsSDV;
    
    NormalPrior normalPrior(parametersMean, parametersSDV);
    ptrPriors[0] = &normalPrior;
    */ 

    // Super-Gaussian Prior
    /*
    ArrayXd parametersMean(Ndimensions);
    ArrayXd parametersSDV(Ndimensions);
    ArrayXd parametersWOP(Ndimensions);
    ArrayXd DeltaNuMean = ArrayXd::Zero(NradialOrders);
    ArrayXd DeltaNuSDV = ArrayXd::Zero(NradialOrders);
    ArrayXd DeltaNuWOP = ArrayXd::Zero(NradialOrders);
    ArrayXd deltaNu02Mean = ArrayXd::Zero(NradialOrders);
    ArrayXd deltaNu02SDV = ArrayXd::Zero(NradialOrders);
    ArrayXd deltaNu02WOP = ArrayXd::Zero(NradialOrders);
    ArrayXd heightsMean = ArrayXd::Zero(NradialOrders*NangularDegrees);
    ArrayXd heightsSDV = ArrayXd::Zero(NradialOrders*NangularDegrees);
    ArrayXd heightsWOP = ArrayXd::Zero(NradialOrders*NangularDegrees);
    double nuMaxMean = 26.0;
    double nuMaxSDV = 30.0;
    double nuMaxWOP = 30.0;
    double noiseLevelMean = 0.1;
    double noiseLevelSDV = 20;
    double noiseLevelWOP = 20;
    
    parametersMean(0) = nuMaxMean;
    parametersMean(1) = noiseLevelMean; 
    DeltaNuMean += 2.5;
    deltaNu02Mean += 0.1;
    heightsMean += 100.0;        
    parametersMean.segment(NglobalParameters, NradialOrders) = DeltaNuMean;
    parametersMean.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02Mean;
    parametersMean.segment(NglobalParameters + 2*NradialOrders, NradialOrders*NangularDegrees) = heightsMean;
    
    parametersSDV(0) = nuMaxSDV;
    parametersSDV(1) = noiseLevelSDV; 
    DeltaNuSDV += 2.5;
    deltaNu02SDV += 0.1;
    heightsSDV += 100.0;        
    parametersSDV.segment(NglobalParameters, NradialOrders) = DeltaNuSDV;
    parametersSDV.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02SDV;
    parametersSDV.segment(NglobalParameters + 2*NradialOrders, NradialOrders*NangularDegrees) = heightsSDV;
    
    parametersWOP(0) = nuMaxWOP;
    parametersWOP(1) = noiseLevelWOP; 
    DeltaNuWOP += 2.5;
    deltaNu02WOP += 0.1;
    heightsWOP += 100.0;        
    parametersWOP.segment(NglobalParameters, NradialOrders) = DeltaNuWOP;
    parametersWOP.segment(NglobalParameters + NradialOrders, NradialOrders) = deltaNu02WOP;
    parametersWOP.segment(NglobalParameters + 2*NradialOrders, NradialOrders*NangularDegrees) = heightsWOP;
    
    SuperGaussianPrior superGaussianPrior(parametersMean, parametersSDV, parametersWOP);
    ptrPriors[0] = &superGaussianPrior;
*/


    

    // Third step - Set up the likelihood function to be used
    
    ExponentialLikelihood likelihood(observations, model);
    

    // Fourth step - Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 20;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Start nested sampling process
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int Nobjects = 200;                            // Number of active points evolving within the nested sampling process. 
    int maxNdrawAttempts = 10000;                   // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = 1000;  // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = 50;         // Clustering is only happening every X iterations.
    double initialEnlargementFraction = 0.50;       // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.8;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.05;                // Termination factor for nested sampling process.

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        Nobjects, initialEnlargementFraction, shrinkingRate);
    nestedSampler.run(terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


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
