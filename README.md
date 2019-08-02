# PeakBagging extension to DIAMONDS

<p align="center">
<a href="https://github.com/EnricoCorsaro/PeakBagging"><img src="https://img.shields.io/badge/GitHub-PeakBagging-yellow"/></a>
<a href="https://github.com/EnricoCorsaro/PeakBagging/blob/master/LICENSE.txt"><img src="https://img.shields.io/badge/license-CC%20BY--SA-blue"/></a>
<a href='https://diamonds.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/diamonds/badge/?version=latest' alt='Documentation Status' /></a>
<a href="https://github.com/EnricoCorsaro/PeakBagging/issues"><img src="https://img.shields.io/github/issues-closed/EnricoCorsaro/PeakBagging"/></a>
<img width="500" src="https://raw.githubusercontent.com/EnricoCorsaro/DIAMONDS/master/docs/figures/DIAMONDS_LOGO_WHITE.png"/>
</p>

### Authors
- [Enrico Corsaro](mailto:enrico.corsaro@inaf.it)
- [Joris De Ridder](mailto:joris.deridder@kuleuven.be)

### Short description
<div align="justify">
The PeakBagging extension to DIAMONDS is an extension of the Bayesian inference code DIAMONDS to perform the fitting of the oscillation peaks in the power spectrum of a solar-like  or red giant star. The tool can also be used for the peak significance test by means of the Bayesian model comparison. It currently implements several modalities to perform a fit, including the multi-modal approach, e.g. see [Corsaro 2019](https://www.frontiersin.org/articles/10.3389/fspas.2019.00021/full) for more details. In addition, background models can be changed at runtime by choosing the reference name for the background that is desired.
</div>

### Download & Installation
The procedure to retrieve the PeakBagging extension is identical to that of DIAMONDS (see [diamonds.readthedocs.io](http://diamonds.readthedocs.io/) for detailed information), so you can either clone the repository or simply download it as a ZIP file. In this second option, by unpacking the PeakBagging-master.zip file you will find a folder labeled Background-master and containing a structure similar to that of the folder Diamonds. First you need to rename the folder as **PeakBagging**, and place it in the same working directory of Diamonds (not inside the Diamonds folder!). This extension needs to be compiled separately from Diamonds, but only after you have compiled Diamonds first. Diamonds is used as a library for this extension. The compilation commands are the same as for Diamonds.  

**IMPORTANT**: Before proceeding with the compilation of the PeakBagging code make sure you put the PeakBagging folder at the same path level of that of Diamonds. This means that the PeakBagging folder has not to be placed inside the Diamonds folder, but inside the parent directory where you placed Diamonds.

### Documentation
Please make sure you read the documentation at [diamonds.readthedocs.io](http://diamonds.readthedocs.io/) before installing and using the code. This extension requires that the DIAMONDS code is first installed in your system. The installation of the PeakBagging extension is the same as that done for DIAMONDS.

### Tutorials
To run the tutorials provided in the package, please follow the guidelines presented in [tutorials/README.md](https://github.com/EnricoCorsaro/PeakBagging/blob/master/tutorials/README.md)
