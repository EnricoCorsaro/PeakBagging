# Tutorial for the peak fitting and peak significance test of the red giant star KIC 12008916

In this tutorial you will find an example of format for data and configuring files for the same Kepler red giant used in the Background extension tutorial. The oscillation fitting model adopted is represented by a sum of Lorentzian profiles. The figure below provides an example of the resulting fit (in green) to the selected portion of power spectrum of this star using the PeakBagging extension to DIAMONDS. The peak significance test will let you know whether the third oscillation peak (starting from left) is significant or not according to the Bayesian model comparison.

![PeakBagging fit](https://raw.githubusercontent.com/EnricoCorsaro/PeakBagging/master/tutorials/KIC012008916_PeakBagging.png)

To run the tutorial follow the procedure:

1. Move the file `KIC012008916/localPath.txt` into `PeakBagging/build/`
2. Edit the path inside `localPath.txt` to match your local working path
3. Move the file `KIC012008916/KIC012008916.txt` into `PeakBagging/data/`
4. Create the empty directories labeled `0` and `0A` inside `KIC012008916/pb/`
5. Move the folder `KIC012008916` (and all its content) under `PeakBagging/results/`
6. Go to `PeakBagging/build/`
7. Execute the code for the first run (including all peaks) by using the command line 
`./peakbagging KIC 012008916 pb 0 ThreeHarvey prior_hyperParameters -1 0 0`
8. Execute the code for the second run (excluding the peak to be tested) by using the command line 
`./peakbagging KIC 012008916 pb 0A ThreeHarvey prior_hyperParameters -1 0 0`
9. Compare the output Bayesian evidences from the two models to check whether the peak is significant or not
