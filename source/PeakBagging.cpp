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
#include "PuntoBackgroundModel.h"
#include "PowerlawReducer.h"
#include "FerozReducer.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
    
    if (argc != 3)
    {
        cerr << "Usage: peakbagging <input file> <output directory>" << endl;
        exit(EXIT_FAILURE);
    }

    
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    string inputFileName(argv[1]);
    string outputDirName(argv[2]);
    string outputPathPrefix = outputDirName;

    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

   
    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);

    
    // Trim input dataset if desired

    double lowerFrequency = 1635.0;        // muHz
    double upperFrequency = 1725.0;        // muHz
    vector<int> trimIndices = Functions::findArrayIndicesWithinBoundaries(covariates, lowerFrequency, upperFrequency);
    int Nbins = trimIndices.size();
    ArrayXd trimmedArray(Nbins);
    trimmedArray = covariates.segment(trimIndices[0],Nbins);
    covariates.resize(Nbins);
    covariates = trimmedArray;
    trimmedArray = observations.segment(trimIndices[0],Nbins);
    observations.resize(Nbins);
    observations = trimmedArray;

    cerr << " Frequency range: [" << covariates.minCoeff() << ", " << covariates.maxCoeff() << "] muHz" << endl;
    cerr << endl;


    // --------------------------------------------------
    // ----- Zeroth step. Set up problem dimensions -----
    // --------------------------------------------------

    vector<int> NparametersPerType(2);                              // Vector containing the configuring parameters for the LorentzianMixture model
    int NprofileParameters = 3;
    NparametersPerType[0] = NprofileParameters;                     // Number of parameters determining the shape 
                                                                    // of the mode profile (central frequency, height, linewidth, inclination angle, etc.)
    int Nmodes = 3; 
    NparametersPerType[1] = Nmodes;                                 // Total number of modes to fit
    int Ndimensions = NprofileParameters*Nmodes;

    
    
    // -------------------------------------------------------------------
    // ----- First step. Set up the models for the inference problem ----- 
    // -------------------------------------------------------------------
    
    // Chose a model for background and configure it, then do the same for the peakbagging model.
    
    PuntoBackgroundModel backgroundModel(covariates);             // CAN model by Kallinger et al. 2010

    string backgroundConfiguringParameters = outputDirName + "Punto_configuringParameters.txt";
    backgroundModel.readConfiguringParametersFromFile(backgroundConfiguringParameters);
    LorentzianMixtureModel model(covariates, NparametersPerType, backgroundModel);


    // -------------------------------------------------------
    // ----- Second step. Set up all prior distributions -----
    // -------------------------------------------------------
    
    // Total number of prior types

    int NpriorTypes = 1;
    vector<Prior*> ptrPriors(NpriorTypes);


    // Uniform Prior

    ArrayXd parametersMinima(Ndimensions);                      // Minima
    ArrayXd parametersMaxima(Ndimensions);                      // Maxima
    
    ArrayXd centralFrequencyMinima(Nmodes);                     // Central frequency of oscillation mode
    ArrayXd centralFrequencyMaxima(Nmodes);                     
    ArrayXd amplitudeMinima(Nmodes);                         // Natural logarithm of mode height
    ArrayXd amplitudeMaxima(Nmodes);            
    ArrayXd linewidthsMinima(Nmodes);                           // Mode Linewidth (FWHM)
    ArrayXd linewidthsMaxima(Nmodes);                           


    // l = 1
    centralFrequencyMinima(0) = 1650.0;
    centralFrequencyMaxima(0) = 1670.0;
    amplitudeMinima(0) = 5.5;
    amplitudeMaxima(0) = 8.0;
    linewidthsMinima(0) = 4.0;
    linewidthsMaxima(0) = 10.0;
    
    centralFrequencyMinima(1) = 1693.0;
    centralFrequencyMaxima(1) = 1700.0;
    amplitudeMinima(1) = 3.0;
    amplitudeMaxima(1) = 8.0;
    linewidthsMinima(1) = 4.0;
    linewidthsMaxima(1) = 8.0;
    
    centralFrequencyMinima(2) = 1700.1;
    centralFrequencyMaxima(2) = 1710.0;
    amplitudeMinima(2) = 5.5;
    amplitudeMaxima(2) = 8.0;
    linewidthsMinima(2) = 4.0;
    linewidthsMaxima(2) = 10.0;

    for (int i=0; i < Nmodes; ++i)
    {
        parametersMinima.segment(Nmodes*i, NprofileParameters) << centralFrequencyMinima(i), amplitudeMinima(i), linewidthsMinima(i); 
        parametersMaxima.segment(Nmodes*i, NprofileParameters) << centralFrequencyMaxima(i), amplitudeMaxima(i), linewidthsMaxima(i); 
    }
    
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    string fullPathHyperParameters = outputPathPrefix + "hyperParametersUniform.txt";
    uniformPrior.writeHyperParametersToFile(fullPathHyperParameters);


    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    ExponentialLikelihood likelihood(observations, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the K-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    int minNclusters = 1;
    int maxNclusters = 6;
    int Ntrials = 10;
    double relTolerance = 0.01;

    EuclideanMetric myMetric;
    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 1000;                      // Maximum (and initial) number of live points evolving within the nested sampling process. 
    int minNobjects = 1000;                          // Minimum number of live points allowed in the computation
    int maxNdrawAttempts = 10000;                   // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = 1000;    // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = 50;        // Clustering is only happening every N iterations.
    double initialEnlargementFraction = 2.00;       // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.02;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.01;                // Termination factor for nested sampling process.

    
    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    
    double tolerance = 1.e2;
    double exponent = 0.4;
    PowerlawReducer livePointsReducer(nestedSampler, tolerance, exponent, terminationFactor);

    nestedSampler.run(livePointsReducer, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, 
                      maxNdrawAttempts, terminationFactor, outputPathPrefix);

    nestedSampler.outputFile << "# List of configuring parameters used for the ellipsoidal sampler and X-means" << endl;
    nestedSampler.outputFile << "# Row #1: Minimum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #2: Maximum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #3: Initial Enlargement Fraction" << endl;
    nestedSampler.outputFile << "# Row #4: Shrinking Rate" << endl;
    nestedSampler.outputFile << minNclusters << endl;
    nestedSampler.outputFile << maxNclusters << endl;
    nestedSampler.outputFile << initialEnlargementFraction << endl;
    nestedSampler.outputFile << shrinkingRate << endl;
    nestedSampler.outputFile.close();


    // -------------------------------------------------------
    // ----- Last step. Save the results in output files -----
    // -------------------------------------------------------
   
    Results results(nestedSampler);
    results.writeParametersToFile("parameter");
    results.writeLogLikelihoodToFile("logLikelihood.txt");
    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    results.writePosteriorProbabilityToFile("posteriorDistribution.txt");
    results.writeLogWeightsToFile("logWeights.txt");

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

    
    // END

    return EXIT_SUCCESS;
}
