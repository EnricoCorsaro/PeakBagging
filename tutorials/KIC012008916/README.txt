Tutorial for the peak fitting and peak significance test of the red giant star KIC 12008916.
To run the tutorial follow the procedure:
1 - Move the files localPath.txt and multirun into PeakBagging/build/
2 - Edit the path inside the localPath.txt to match your local working path
3 - Move the file KIC012008916.txt into PeakBagging/data/
4 - Move the remaining files into PeakBagging/results/KIC012008916/
5 - Go to PeakBagging/build/ and edit the localPath.txt file by inserting with your local working path
6 - Execute the code for the first run (including all peaks) by using the command line ./peakbagging KIC 012008916 pb 0 ThreeHarvey prior_hyperParameters -1 0 0
7 - Execute the code for the second run (excluding the peak to be tested) by using the command line ./peakbagging KIC 012008916 pb 0A ThreeHarvey prior_hyperParameters -1 0 0
8 - Compare the output Bayesian evidences from the two models to check whether the peak is significant or not
