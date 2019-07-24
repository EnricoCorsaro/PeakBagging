// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ INAF-OACT - January 2019
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
#include "ExponentialLikelihood.h"
#include "SingleLorentzianFixedLinewidthModel.h"
#include "DoubleLorentzianFixedLinewidthModel.h"
#include "LorentzianMixtureModel.h"
#include "LorentzianRotationMixtureModel.h"
#include "PeakTestBackgroundModel.h"
#include "PeakTestLorentzianBackgroundModel.h"
#include "PeakTestTwoLorentziansBackgroundModel.h"
#include "PeakTestLorentzianModel.h"
#include "PeakTestLorentzianRotationModel.h"
#include "BackgroundModel.h"
#include "PowerlawReducer.h"
#include "Results.h"
#include "Ellipsoid.h"
#include "PrincipalComponentProjector.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
   
    if (argc != 10)
    {
        cerr << "Usage: ./peakbagging <Catalog ID> <Star ID> <output sub-directory> <run number> <background model> <input prior base filename> " 
        << "<input linewidth> <Peak-Test flag> <Background level output flag>" << endl;
        exit(EXIT_FAILURE);
    }

    
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    string CatalogID(argv[1]);
    string StarID(argv[2]);
    string runNumber(argv[4]);
    string backgroundModelName(argv[5]);
    string inputPriorBaseName(argv[6]);
    string inputLinewidth(argv[7]);
    string inputPeakTestFlag(argv[8]);
    string inputBackgroundLevelOutputFlag(argv[9]);
    double peakLinewidth = stod(inputLinewidth);
    int peakTestFlag = stoi(inputPeakTestFlag);
    int backgroundLevelOutputFlag = stoi(inputBackgroundLevelOutputFlag);


    // Read the local path for the working session from an input ASCII file
    ifstream inputFile;
    File::openInputFile(inputFile, "localPath.txt");
    File::sniffFile(inputFile, Nrows, Ncols);
    vector<string> myLocalPath;
    myLocalPath = File::vectorStringFromFile(inputFile, Nrows);
    inputFile.close();


    // Set up some string paths used in the computation
    string baseInputDirName = myLocalPath[0] + "data/";
    string inputFileName = baseInputDirName + CatalogID + StarID + ".txt";
    string outputSubDirName(argv[3]);
    string baseOutputDirName = myLocalPath[0] + "results/" + CatalogID + StarID + "/";
    string outputDirName = baseOutputDirName + outputSubDirName + "/";
    string outputPathPrefix = outputDirName + runNumber + "/peakbagging_";
   
    cerr << "------------------------------------------------------ " << endl;
    if (peakTestFlag == 1)
    {
            cout << " Performing peak significance test for " + CatalogID + StarID << endl;
    } 
    else
    {
        if (peakLinewidth >= 0.0)
        {
            cout << " Multi-modal Peak Bagging analysis for " + CatalogID + StarID << endl;

            if (peakLinewidth > 0) cout << " Using input linewidth " + inputLinewidth + " microHz." << endl;
            if (peakLinewidth == 0) cout << " Using frequencyResolution as input linewidth." << endl;
        }
        else
        {
            if (peakLinewidth == -1)
            {
            cout << " Uni-modal Peak Bagging analysis without rotation for " + CatalogID + StarID << endl;
            }
            else
            {
            cout << " Uni-modal Peak Bagging analysis with rotation for " + CatalogID + StarID << endl;
            }
        }
    }
    cerr << "------------------------------------------------------ " << endl;
    cout << endl;

   
    // Read the input dataset
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();


    // Creating frequency and PSD arrays
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    double frequencyResolution = covariates(1) - covariates(0);
    double lowerFrequency;
    double upperFrequency;


    // If the fit has to be done in a uni-modal high-dimensional manner, then trim the dataset according to an input range.
    // If a peak testing has to be performed, then skip this part.
    if (peakTestFlag != 1)
    {
        if (peakLinewidth < 0.0)
        {
            // Read input frequency range of the PSD
            inputFileName = outputDirName + "frequencyRange_" + runNumber + ".txt";
            File::openInputFile(inputFile, inputFileName);
            File::sniffFile(inputFile, Nrows, Ncols);
            ArrayXd frequencyRange;
            frequencyRange = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
            inputFile.close();

            lowerFrequency = frequencyRange(0);        // muHz
            upperFrequency = frequencyRange(1);        // muHz

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
        }
    }


    // -------------------------------------------------------
    // ----- First step. Set up all prior distributions ------
    // -------------------------------------------------------
    
    unsigned long Nparameters;              // Number of parameters for which prior distributions are defined

    // ---- Read prior hyper parameters for resolved modes -----
    inputFileName = outputDirName + inputPriorBaseName + "_" + runNumber + ".txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    ArrayXXd hyperParameters;
  
    if (peakTestFlag == 1)
    {
        if (Nparameters != 2 && Nparameters != 3 && Nparameters != 4 && Nparameters != 5 && Nparameters != 7)
        {
            cerr << "When performing a peak significance test, the code allows for the following cases: " << endl;
            cerr << "(A) - Only background (Nparameters = 1). " << endl;
            cerr << "(B) - Background plus Lorentzian profile (Nparameters = 4). " << endl;
            cerr << "(C) - Background plus two Lorentzian profiles (Nparameters = 7). " << endl;
            cerr << "(D) - One Lorentzian profile with fixed background (Nparameters = 3). " << endl;
            cerr << "(E) - One Lorentzian profile with rotationally split components and fixed background (Nparameters = 5). " << endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        if (peakLinewidth >= 0.0)
        {
            if (Nparameters != 2 && Nparameters != 3)
            {
                cerr << "Wrong number of input prior hyper-parameters." << endl;
                cerr << "There are only two free parameters for the single-peak model and three for the double-peak model." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if (peakLinewidth == -1)
            {
                if (fmod(Nparameters,3) != 0.0)
                {
                    cerr << "Wrong number of input prior hyper-parameters for resolved peaks." << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                if (fmod((Nparameters-2),3) != 0.0)
                {
                    cerr << "Wrong number of input prior hyper-parameters for resolved peaks." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
   
    if (Ncols == 1)
    {
        Ncols = 2;
        hyperParameters.conservativeResize(Nparameters,Ncols);
    }

    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    ArrayXd hyperParametersMinima1 = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima1 = hyperParameters.col(1);
    // ---------------------------------------------------------

    int Ndimensions = Nparameters;              // Total number of dimensions of the peak bagging model

    // Uniform Prior
    
    int NpriorTypes = 1;                                        // Total number of prior types included in the computation
    vector<Prior*> ptrPriors(NpriorTypes);
    
    ArrayXd parametersMinima;                      // Minima for prior PDF
    ArrayXd parametersMaxima;                      // Maxima for prior PDF
    
    if (peakTestFlag == 1)
    {
        if (Ndimensions == 2)
        {
            // In this case the peak test contains only the background factor (A)

            parametersMinima.resize(Ndimensions-1);
            parametersMaxima.resize(Ndimensions-1);
            parametersMinima << hyperParametersMinima1.segment(1,Ndimensions-1);
            parametersMaxima << hyperParametersMaxima1.segment(1,Ndimensions-1); 
        }
        
        if (Ndimensions == 3 || Ndimensions == 4 || Ndimensions == 5 || Ndimensions == 6)
        {
            // In this case the peak test contains cases (B), (C), (D), (E)

            parametersMinima.resize(Ndimensions);
            parametersMaxima.resize(Ndimensions);
            parametersMinima << hyperParametersMinima1;
            parametersMaxima << hyperParametersMaxima1; 
        }        
    }
    else
    {
        parametersMinima.resize(Ndimensions);
        parametersMaxima.resize(Ndimensions);
        parametersMinima << hyperParametersMinima1;
        parametersMaxima << hyperParametersMaxima1; 
    }
   
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    string fullPathHyperParameters = outputPathPrefix + "hyperParametersUniform.txt";
    uniformPrior.writeHyperParametersToFile(fullPathHyperParameters);

    
    // If the fit has to be done in a multi-modal low-dimensional manner, then trim the dataset according to the input prior for the frequency.
    // Do this also in the case of a peak significance test. 
    
    if (peakLinewidth >= 0.0 || peakTestFlag == 1) 
    {
        // Trim input dataset in the given frequency range
        
        lowerFrequency = hyperParameters(0,0);        // muHz
        
        if (Ndimensions == 7)
        {
            // For the peak test model having two Lorentzian profiles
            // take as the lower frequency of the range the lower
            // prior bound of the first frequency centroid and as 
            // the upper frequency of the range, the upper prior
            // bound of the frequency centroid of the second
            // Lorentzian profile (the one at higher frequency).

            upperFrequency = hyperParameters(3,1);        // muHz
        }
        else
        {
            upperFrequency = hyperParameters(0,1);        // muHz
        }

        vector<int> trimIndices = Functions::findArrayIndicesWithinBoundaries(covariates, lowerFrequency, upperFrequency);
        int Nbins = trimIndices.size();
        ArrayXd trimmedArray(Nbins);
        trimmedArray = covariates.segment(trimIndices[0],Nbins);
        covariates.resize(Nbins);
        covariates = trimmedArray;
        trimmedArray = observations.segment(trimIndices[0],Nbins);
        observations.resize(Nbins);
        observations = trimmedArray;
    }

    cout << " Frequency range: [" << setprecision(4) << covariates.minCoeff() << ", " 
         << covariates.maxCoeff() << "] muHz" << endl;
    cout << " Frequency resolution: " << setprecision(4) << frequencyResolution << " muHz" << endl;
    cout << endl;
    

    // -------------------------------------------------------------------
    // ---- Second step. Set up the background model ---------------------
    // -------------------------------------------------------------------

    BackgroundModel backgroundModel(covariates, backgroundModelName);
    inputFileName = baseOutputDirName + "NyquistFrequency.txt";
    backgroundModel.readNyquistFrequencyFromFile(inputFileName);
    string backgroundConfiguringParameters = baseOutputDirName + "backgroundParameters.txt";
    backgroundModel.readConfiguringParametersFromFile(backgroundConfiguringParameters);
    
    if (backgroundLevelOutputFlag == 1)
    {
        string outputBackgroundLevelFileName = baseOutputDirName + "backgroundLevel.txt";
        backgroundModel.writeBackgroundPredictionToFile(outputBackgroundLevelFileName);
    }
   

    // -------------------------------------------------------------------
    // ---- Third step. Set up the peak bagging model --------------------
    // -------------------------------------------------------------------
   
    Model *model = nullptr;

    if (peakTestFlag == 1)
    {
        // Set the model to perform the peak significance test. This can be one of the
        // four models (A), (B), (C), (D), (E), depending on the number of
        // dimensions as set out from the input prior list.
       
        if (Nparameters == 2)
        {
            model = new PeakTestBackgroundModel(covariates, backgroundModel);
        }
        
        if (Nparameters == 3)
        {
            model = new PeakTestLorentzianModel(covariates, backgroundModel);
        }

        if (Nparameters == 4)
        {
            model = new PeakTestLorentzianBackgroundModel(covariates, backgroundModel);
        }

        if (Nparameters == 5)
        {
            model = new PeakTestLorentzianRotationModel(covariates, backgroundModel);
        }

        if (Nparameters == 7)
        {
            model = new PeakTestTwoLorentziansBackgroundModel(covariates, backgroundModel);
        }
    }
    else
    {
        if (peakLinewidth >= 0)
        {
            if (peakLinewidth > 0)
            {
                // Set the linewidth of the profile for the multi-modal fitting to an input value > 0.

                if (Nparameters == 2) 
                {
                    model = new SingleLorentzianFixedLinewidthModel(covariates, peakLinewidth, backgroundModel);
                }
                else
                {
                    model = new DoubleLorentzianFixedLinewidthModel(covariates, peakLinewidth, backgroundModel);
                }
            }
            else
            {
                // Use the frequency resolution as peak linewidth for the multi-modal fitting. Useful in case
                // narrow (unresolved) mixed modes are expected.

                if (Nparameters == 2) 
                {
                    model = new SingleLorentzianFixedLinewidthModel(covariates, frequencyResolution, backgroundModel);
                }
                else
                {
                    model = new DoubleLorentzianFixedLinewidthModel(covariates, frequencyResolution, backgroundModel);
                }
            }
        }
        else
        {
            int Nresolved;

            if (peakLinewidth == -1)
            {
                // No input linewidth has to be used (< 0). In this case perform the standard uni-modal fitting 
                // by assuming a mixture of Lorentzian profiles without rotational splitting.

                Nresolved = Ndimensions/3;
                model = new LorentzianMixtureModel(covariates, Nresolved, backgroundModel);
            }
            else
            {
                // No input linewidth has to be used (< 0). In this case perform the standard uni-modal fitting 
                // by assuming a mixture of Lorentzian profiles with rotational splitting.
                // This requires an input list of angular degrees to be fed into the peak bagging model.

                Nresolved = (Ndimensions-2)/3;
                string angularDegreesFileName = outputDirName + "angularDegrees_" + runNumber + ".txt";
                model = new LorentzianRotationMixtureModel(covariates, Nresolved, angularDegreesFileName, backgroundModel);
            }
        }
    }
    

    // -----------------------------------------------------------------
    // ---- Fourth step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    ExponentialLikelihood likelihood(observations, *model);
    

    // -------------------------------------------------------------------------------
    // ----- Fifth step. Set up the X-means clusterer using an Euclidean metric ------
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
    double relTolerance = 0.01;     // k-means
  
    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;
    EuclideanMetric myMetric;
    
    KmeansClusterer clusterer(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------

    inputFileName = outputDirName + "NSMC_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    configuringParameters.setZero();
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    if (Nparameters > 9 || Nparameters < 8)
    {
        cerr << "Wrong number of input parameters for NSMC algorithm." << endl;
        cerr << "There must be either 8 or 9 parameters (the latter specified only if " << 
        "a fixed number of nested iterations is required, i.e. for the multi-modal fit." << endl;
        exit(EXIT_FAILURE);
    }
    

    // Print results on the screen
    
    bool printOnTheScreen = true;                  


    // Initial number of live points
    
    int initialNobjects = configuringParameters(0);

    
    // Minimum number of live points 
    
    int minNobjects = configuringParameters(1);
    
    
    // Maximum number of attempts when trying to draw a new sampling point
    
    int maxNdrawAttempts = configuringParameters(2);

    
    // The first N iterations, we assume that there is only 1 cluster
    
    int NinitialIterationsWithoutClustering = configuringParameters(3);

    
    // Clustering is only happening every N iterations.
    
    int NiterationsWithSameClustering = configuringParameters(4);

    
    // Fraction by which each axis in an ellipsoid has to be enlarged
    // It can be a number >= 0, where 0 means no enlargement. configuringParameters(5)
    // Calibration from Corsaro et al. (2018)
    
    double initialEnlargementFraction = 0.369*pow(Ndimensions,0.574);    

    
    // Exponent for remaining prior mass in ellipsoid enlargement fraction.
    // It is a number between 0 and 1. The smaller the slower the shrinkage of the ellipsoids.
    
    double shrinkingRate = configuringParameters(6);        
                                                                                                                    
    
    // Termination factor for nested sampling process.                                 
    
    double terminationFactor = configuringParameters(7);

    
    // Total maximum number of nested iterations required to carry out the computation.
    // This is used only in the multi-modal approach.
    
    int maxNiterations = 0; 
    if (peakLinewidth >= 0 && peakTestFlag != 1)
    {
        maxNiterations = configuringParameters(8);
    }
    

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, clusterer, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    
    double tolerance = 1.e2;
    double exponent = 0.4;
    PowerlawReducer livePointsReducer(nestedSampler, tolerance, exponent, terminationFactor);

    nestedSampler.run(livePointsReducer, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, 
                      maxNdrawAttempts, terminationFactor, maxNiterations, outputPathPrefix);

    nestedSampler.outputFile << "# List of configuring parameters used for the ellipsoidal sampler and X-means" << endl;
    nestedSampler.outputFile << "# Row #1: Minimum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #2: Maximum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #3: Initial Enlargement Fraction" << endl;
    nestedSampler.outputFile << "# Row #4: Shrinking Rate" << endl;
    nestedSampler.outputFile << minNclusters << endl;
    nestedSampler.outputFile << maxNclusters << endl;
    nestedSampler.outputFile << initialEnlargementFraction << endl;
    nestedSampler.outputFile << shrinkingRate << endl;
    nestedSampler.outputFile << "# Other information on the run" << endl;
    nestedSampler.outputFile << "# Row #1: Local working path used" << endl;
    nestedSampler.outputFile << "# Row #2: Catalog and Star ID" << endl;
    nestedSampler.outputFile << "# Row #3: Run Directory" << endl;
    nestedSampler.outputFile << "# Row #4: Run Number" << endl;
    nestedSampler.outputFile << "# Row #5: Peak significance test (from A to E, 0 if not a test)" << endl;
    nestedSampler.outputFile << "# Row #6: Background model adopted" << endl;
    nestedSampler.outputFile << "# Row #7: Peak linewidth in muHz (0 if uni-modal peak bagging or peak test)" << endl;
    nestedSampler.outputFile << myLocalPath[0] << endl;
    nestedSampler.outputFile << CatalogID + StarID << endl;
    nestedSampler.outputFile << outputSubDirName << endl;
    nestedSampler.outputFile << runNumber << endl;

    if (peakTestFlag == 1)
    {
        if (Ndimensions == 2)
        {
            nestedSampler.outputFile << "A" << endl;
        }
        
        if (Ndimensions == 4)
        {
            nestedSampler.outputFile << "B" << endl;
        }

        if (Ndimensions == 7)
        {
            nestedSampler.outputFile << "C" << endl;
        }

        if (Ndimensions == 3)
        {
            nestedSampler.outputFile << "D" << endl;
        }

        if (Ndimensions == 5)
        {
            nestedSampler.outputFile << "E" << endl;
        }
    }
    else
    {
        nestedSampler.outputFile << 0 << endl;
    }

    nestedSampler.outputFile << backgroundModelName << endl;
    
    if (peakLinewidth < 0 || peakTestFlag == 1)
    {
        nestedSampler.outputFile << setprecision(2) << 0.0 << endl;
    }
    else
    {
        if (peakLinewidth > 0)
        {
            nestedSampler.outputFile << setprecision(3) << peakLinewidth << endl;
        }
        else
        {
            nestedSampler.outputFile << setprecision(3) << frequencyResolution << endl;
        }
    }

    nestedSampler.outputFile.close();


    // -------------------------------------------------------
    // ----- Last step. Save the results in output files -----
    // -------------------------------------------------------
   
    Results results(nestedSampler);

    if (peakTestFlag != 1)
    {
        results.writeParametersToFile("parameter");
        results.writeLogLikelihoodToFile("logLikelihood.txt");
        results.writePosteriorProbabilityToFile("posteriorDistribution.txt");
        results.writeLogEvidenceToFile("logEvidence.txt");
        results.writeLogMeanLiveEvidenceToFile("logMeanLiveEvidence.txt");
    }

    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    

    // Print out parameter estimates only in the case of a uni-modal high-dimensional fit.

    if ((peakLinewidth < 0.0 && peakTestFlag != 1) || (peakTestFlag == 1 && Ndimensions == 5))
    {
        double credibleLevel = 68.3;
        bool writeMarginalDistributionToFile = true;
        results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);
    }

    cerr << "Process #" << runNumber << " under subdir: " + outputSubDirName + " has been completed." << endl;

    return EXIT_SUCCESS;
}
