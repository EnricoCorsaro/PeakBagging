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
```bash 
./peakbagging KIC 012008916 pb 0 ThreeHarvey prior_hyperParameters -1 0 0 0
```
8. Execute the code for the second run (excluding the peak to be tested) by using the command line
```bash
./peakbagging KIC 012008916 pb 0A ThreeHarvey prior_hyperParameters -1 0 0 0
```
9. Compare the output Bayesian evidences from the two models to check whether the peak is significant or not
10. Once the computation is completed, you can plot the results with Python by using `plot_peakbagging.py` provided in the tutorials folder. The routine is set up in a way that it will plot the results from your folder `0`. Please make sure that all paths set inside the Python routines match correctly with your actual working paths for PeakBagging.


# Tutorial for the multi-modal fitting of the red giant star KIC 12008916

In this tutorial you will be able to perform a multi-modal fit using a single Lorentzian profile over the same frequency range of the previous tutorial. The figure below shows the resulting multi-modal sampling obtained by DIAMONDS, which can be plotted using the output peakbagging_parameter000.txt file of the computation, namely the one containing the sampling evolution of the frequency centroid parameter of the Lorentzian profile used for the fit.

**WARNING**: make sure to run this tutorial only once you have already executed the steps #1, #2, #3, #5 from the previous tutorial.

<img width="500" src="https://raw.githubusercontent.com/EnricoCorsaro/PeakBagging/master/tutorials/KIC012008916_Islands.png"/>
</p>

To run the tutorial follow the procedure:

1. Create an empty directory labeled `0` inside `KIC012008916/isla/`
2. Go to `PeakBagging/build/`
3. Execute the code for the multi-modal functionality by using the command line
```bash
./peakbagging KIC 012008916 isla 0 ThreeHarvey prior_hyperParameters 0.034 0 0 0
```

In this case 0.034 is the linewidth (in microHz) used for the Lorentzian profile, while the remaining flags are kept deactivated. The plot shows how a single Lorentzian profile is capable of reproducing the oscillation peak structures that are present in the data within the inspected frequency range (compare this sampling with the plot presented in the first tutorial).