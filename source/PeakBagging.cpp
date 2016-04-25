// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ CEA - January 2016
// e-mail: emncorsaro@gmail.com
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
#include "StandardBackgroundModel.h"
#include "FullBackgroundModel.h"
#include "PowerlawReducer.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
    
    if (argc != 4)
    {
        cerr << "Usage: peakbagging <KIC ID> <output sub-directory> <run number>" << endl;
        exit(EXIT_FAILURE);
    }

    
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    string KICID(argv[1]);
    string runNumber(argv[3]);


    // Read the local path for the working session from an input ASCII file
    ifstream inputFile;
    File::openInputFile(inputFile, "localPath.txt");
    File::sniffFile(inputFile, Nrows, Ncols);
    vector<string> myLocalPath;
    myLocalPath = File::vectorStringFromFile(inputFile, Nrows);
    inputFile.close();


    // Set up some string paths used in the computation
    string baseInputDirName = myLocalPath[0] + "data/";
    string inputFileName = baseInputDirName + "KIC" + KICID + ".txt";
    string outputSubDirName(argv[2]);
    string baseOutputDirName = myLocalPath[0] + "results/KIC" + KICID + "/";
    string outputDirName = baseOutputDirName + outputSubDirName + "/";
    string outputPathPrefix = outputDirName + runNumber + "/peakbagging_";

   
    // Read the input dataset
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

    cout << " Frequency range: [" 
         << setprecision(4) 
         << covariates.minCoeff() 
         << ", " 
         << covariates.maxCoeff() 
         << "] muHz" << endl;
    cout << endl;

    
    // -------------------------------------------------------
    // ----- First step. Set up all prior distributions -----
    // -------------------------------------------------------
    
    unsigned long Nparameters;              // Number of parameters for which prior distributions are defined
    
    // ---- Read prior hyper parameters for resolved modes -----
    inputFileName = outputDirName + "resolvedPeaks_hyperParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    ArrayXXd hyperParameters;
    
    if (fmod(Nparameters,3) != 0.0)
    {
        cerr << "Wrong number of input prior hyper-parameters for resolved peaks." << endl;
        exit(EXIT_FAILURE);
    }
   
    if (Ncols == 1)
    {
        Ncols = 2;
        hyperParameters.conservativeResize(Nparameters,Ncols);
    }

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
  
    if (fmod(Nparameters,2) != 0.0)
    {
        cerr << "Wrong number of input prior hyper-parameters for unresolved peaks." << endl;
        exit(EXIT_FAILURE);
    }

    if (Ncols == 1)
    {
        Ncols = 2;
        hyperParameters.conservativeResize(Nparameters,Ncols);
    }
   
    hyperParameters.setZero();
    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    int Nunresolved = Nparameters/2;
    ArrayXd hyperParametersMinima2 = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima2 = hyperParameters.col(1);
    // ---------------------------------------------------------
   

    int Ndimensions = 3*Nresolved + 2*Nunresolved;              // Total number of dimensions of the peak bagging model
    double frequencyResolution = covariates(1)-covariates(0);   // The frequency bin width of the input data in muHz

    // Uniform Prior
    int NpriorTypes = 1;                                        // Total number of prior types included in the computation
    vector<Prior*> ptrPriors(NpriorTypes);
    
    ArrayXd parametersMinima(Ndimensions);                      // Minima
    ArrayXd parametersMaxima(Ndimensions);                      // Maxima
    parametersMinima << hyperParametersMinima1, hyperParametersMinima2;
    parametersMaxima << hyperParametersMaxima1, hyperParametersMaxima2; 
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    
    // -------------------------------------------------------------------
    // ---- Second step. Set up the models for the inference problem ----- 
    // -------------------------------------------------------------------
    
    // Chose a model for background and configure it, then do the same for the peakbagging model.
    inputFileName = baseOutputDirName + "NyquistFrequency.txt";
    StandardBackgroundModel backgroundModel(covariates, inputFileName);             // Model by Kallinger et al. 2014 - No colored-noise component
    // FullBackgroundModel backgroundModel(covariates, inputFileName);             // Model by Kallinger et al. 2014 - Colored-noise component included
    string backgroundConfiguringParameters = baseOutputDirName + "backgroundParameters.txt";
    backgroundModel.readConfiguringParametersFromFile(backgroundConfiguringParameters);
    LorentzianSincMixtureModel model(covariates, Nresolved, Nunresolved, frequencyResolution, backgroundModel);


    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    ExponentialLikelihood likelihood(observations, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the X-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    inputFileName = outputDirName + "Xmeans_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);

    if (Nparameters != 2)
    {
        cerr << "Wrong number of input parameters for X-means algorithm." << endl;
        exit(EXIT_FAILURE);
    }

    ArrayXd configuringParameters;
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    int minNclusters = configuringParameters(0);
    int maxNclusters = configuringParameters(1);
    
    if ((minNclusters <= 0) || (maxNclusters <= 0) || (maxNclusters < minNclusters))
    {
        cerr << "Minimum or maximum number of clusters cannot be <= 0, and " << endl;
        cerr << "minimum number of clusters cannot be larger than maximum number of clusters." << endl;
        exit(EXIT_FAILURE);
    }
    
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

    if (Nparameters != 8)
    {
        cerr << "Wrong number of input parameters for NSMC algorithm." << endl;
        exit(EXIT_FAILURE);
    }
    

    bool printOnTheScreen = true;                       // Print results on the screen
    int initialNobjects = configuringParameters(0);     // Initial number of live points 
    int minNobjects = configuringParameters(1);         // Minimum number of live points 
    
    if (Ndimensions >= 20)
    {
        initialNobjects = 1000;
        minNobjects = 1000;
    }

    int maxNdrawAttempts = configuringParameters(2);    // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = configuringParameters(3); // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = configuringParameters(4);       // Clustering is only happening every N iterations.
    double initialEnlargementFraction = configuringParameters(5);   //0.267*pow(Ndimensions,0.643);   // configuringParameters(5);   // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = configuringParameters(6);        // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                            // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                            // of the ellipsoids.
                                                            
    if ((shrinkingRate > 1) || (shrinkingRate) < 0)
    {
        cerr << "Shrinking Rate for ellipsoids must be in the range [0, 1]. " << endl;
        exit(EXIT_FAILURE);
    }

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
    nestedSampler.outputFile << "# Local working path used: " + myLocalPath[0] << endl;
    nestedSampler.outputFile << "# KIC ID: " + KICID << endl;
    nestedSampler.outputFile << "# Run Directory: " + outputSubDirName << endl;
    nestedSampler.outputFile << "# Run Number: " + runNumber << endl;
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
 
    cerr << "Process #" << runNumber << " has been completed." << endl;

    return EXIT_SUCCESS;
}
