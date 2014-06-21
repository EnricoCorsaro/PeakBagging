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
#include "LorentzianSincMixtureModel.h"
#include "RedGiantBackgroundModel.h"
#include "PowerlawReducer.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
    
    if (argc != 3)
    {
        cerr << "Usage: peakbagging <KIC ID> <output sub-directory>" << endl;
        exit(EXIT_FAILURE);
    }

    
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    string KICID(argv[1]);
    string baseInputDirName = "/Users/eco/asismology/github/PeakBagging/data/";
    string inputFileName = baseInputDirName + "KIC" + KICID + ".txt";
    string outputSubDirName(argv[2]);
    string baseOutputDirName = "/Users/eco/asismology/github/PeakBagging/results/KIC" + KICID + "/";
    string outputDirName = baseOutputDirName + outputSubDirName + "/";
    string outputPathPrefix = outputDirName + "peakbagging_";

    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

    // Creating frequency and PSD arrays
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    
    // Read input frequency range of the PSD
    inputFileName = outputDirName + "frequencyRange.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXd frequencyRange;
    frequencyRange = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

    double lowerFrequency = frequencyRange(0);        // muHz
    double upperFrequency = frequencyRange(1);        // muHz
    
    // Trim input dataset in the given frequency range
    vector<int> trimIndices = Functions::findArrayIndicesWithinBoundaries(covariates, lowerFrequency, upperFrequency);
    int Nbins = trimIndices.size();
    ArrayXd trimmedArray(Nbins);
    trimmedArray = covariates.segment(trimIndices[0],Nbins);
    covariates.resize(Nbins);
    covariates = trimmedArray;
    trimmedArray = observations.segment(trimIndices[0],Nbins);
    observations.resize(Nbins);
    observations = trimmedArray;

    cerr << " Frequency range: [" 
         << setprecision(4) 
         << covariates.minCoeff() 
         << ", " 
         << covariates.maxCoeff() 
         << "] muHz" << endl;
    cerr << endl;


    // -------------------------------------------------------
    // ----- First step. Set up all prior distributions -----
    // -------------------------------------------------------
    
    // Uniform Prior
    unsigned long Nparameters;              // Number of parameters for which prior distributions are defined
    
    // ---- Read prior hyper parameters for resolved modes -----
    inputFileName = outputDirName + "resolvedPeaks_hyperParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    
    ArrayXXd hyperParameters;
    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    int Nresolved = Nparameters/3;
    ArrayXd hyperParametersMinima1 = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima1 = hyperParameters.col(1);
    // ---------------------------------------------------------


    // --- Read prior hyper parameters for unresolved modes ----
    inputFileName = outputDirName + "unresolvedPeaks_hyperParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    
    hyperParameters.setZero();
    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    int Nunresolved = Nparameters/2;
    ArrayXd hyperParametersMinima2 = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima2 = hyperParameters.col(1);
    // ---------------------------------------------------------


    int Ndimensions = 3*Nresolved + 2*Nunresolved;              // Total number of dimensions of the peak bagging model
    double frequencyResolution = covariates(1)-covariates(0);   // The frequency bin width of the input data in muHz

    int NpriorTypes = 1;                                        // Total number of prior types included in the computation
    vector<Prior*> ptrPriors(NpriorTypes);
    
    ArrayXd parametersMinima(Ndimensions);                      // Minima
    ArrayXd parametersMaxima(Ndimensions);                      // Maxima

    parametersMinima << hyperParametersMinima1, hyperParametersMinima2;
    parametersMaxima << hyperParametersMaxima1, hyperParametersMaxima2; 
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    
    // -------------------------------------------------------------------
    // ----- Second step. Set up the models for the inference problem ----- 
    // -------------------------------------------------------------------
    
    // Chose a model for background and configure it, then do the same for the peakbagging model.
    
    RedGiantBackgroundModel backgroundModel(covariates);             // Model by Kallinger et al. 2014

    string backgroundConfiguringParameters = baseOutputDirName + "backgroundParameters.txt";
    backgroundModel.readConfiguringParametersFromFile(backgroundConfiguringParameters);
    LorentzianSincMixtureModel model(covariates, Nresolved, Nunresolved, frequencyResolution, backgroundModel);


    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    ExponentialLikelihood likelihood(observations, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the K-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------


    inputFileName = outputDirName + "Xmeans_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    ArrayXd configuringParameters;
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    int minNclusters = configuringParameters(0);
    int maxNclusters = configuringParameters(1);
    int Ntrials = 10;
    double relTolerance = 0.01;

    EuclideanMetric myMetric;
    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------

    inputFileName = outputDirName + "NSMC_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    configuringParameters.setZero();
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    bool printOnTheScreen = true;                       // Print results on the screen
    int initialNobjects = configuringParameters(0);     // Initial number of live points 
    int minNobjects = configuringParameters(1);         // Minimum number of live points 
    int maxNdrawAttempts = configuringParameters(2);    // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = configuringParameters(3); // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = configuringParameters(4);       // Clustering is only happening every N iterations.
    double initialEnlargementFraction = configuringParameters(5);   // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = configuringParameters(6);        // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                            // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                            // of the ellipsoids.
    double terminationFactor = configuringParameters(7);    // Termination factor for nested sampling process.

    
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

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

    
    // END

    return EXIT_SUCCESS;
}
