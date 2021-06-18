# CSBP
 Matlab code of the Compressive Sensing BackProjection for earthquake source imaging.

**Compressive Sensing BackProjection (CSBP) is a high resolution backprojection method in frequency domain. It seeks for the sparse representation for the backprojection image using L1-norm inversion method. The details about the method can be found in the following papers, please cite them if you are using this code for your research.**

*Yao, Huajian, et al. "Compressive sensing of the Tohoku‐Oki Mw 9.0 earthquake: Frequency‐dependent rupture modes." Geophysical Research Letters 38.20 (2011).*

*Yin, Jiuxun, Marine A. Denolle, and Huajian Yao. "Spatial and temporal evolution of earthquake dynamics: case study of the Mw 8.3 Illapel Earthquake, Chile." Journal of Geophysical Research: Solid Earth 123.1 (2018): 344-367.*

Yao et al. (2011) is about the basic CSBP method. Yin et al. (2018) presents an improved CSBP method applying grid-refinement to improve the efficicency and resolution of the method. In Yin et al. (2018), a waveform realignment technique is also proposed, but that work best for unilateral rupture. So for general use, this package only include the one with grid-refinement.


==============================================================================

This CSBP code is based on Matlab, and we use a CVX package to solve for the sparsity inversion problem (CVX has both Matlab and Python version). In this package, I'm using the 2011 Tohoku earthquake as an example. Steps to run the Compressive Sensing Backprojection are shown below:

**Step 0: Install the CVX package (Step0_install_cvx_package.m)**


**Step 1: Prepare the input data. (Step1_Prepare_sacdata.m)**

Here the input is the event SAC data, which includes waveforms 60s before and 300s after the P-arrival of the 2011 Tohoku earthquake recorded by USArray (TA). The data has gone through some basic preprocessing including removing instrumental response, de-mean and de-trend. I often directly download the event data from IRIS Wilber3 for convenience. But starting from the continuous raw data is easy to implement.

The preparation includes combining the array data, filtering, resampling, and realigning the data based on the P wave arrival. Finally the processed data is output to a working folder, where the CSBP will start. In this code, the output is a binary format defined by Huajian, but basically the processed data is a matrix of array waveform data and some information about the array. For now we don't have to worry about the format.


**Step 2: Run the CSBP. (Step2_CSBP_AG.m)**

Before run the CSBP, we need to first prepare an inputfile file. Here I simply write a .txt file to include the necessary parameters. The parameter file is in the folder CS_input_files. Parameters include location of epicenter, grid information in the source region, information about time window and frequency, and some station selection thresholds. All of the parameters are commented with explanations. 

Both CSBP (srcGridSpec) and beamforming (srcGridBeam) results are organized as 4D arrays ([nlon, nlat, nfrequency, time-window]) and will be output to the folder ResultsFile in the working folder as srcCVXresult.mat. Other information such as arrays of latitude, longitude, frequencies and time windows are also included in the same .mat file.


**Step 3: Visualize the CSBP results. (Step3_Plot_results.m)**

A plotting script to visualize the CSBP results (and also Beamforming results) is provided here. Images of the CSBP results at different single frequencies (ftest = xxx), or summation of all frequencies within a frequency band (sumfreqrange =1) can be made.

Because CSBP method is based on sparse inversion, the direct results from CSBP are very sparse points in the source region that present the optimal locations where the waves are radiated. This is unlike the IMAGES from other BP methods. To look at the spatial distribution of BP images, we can choose to apply a spatial Gaussian smoothing (indexSpatialSmooth, SmoothRadius) to the CSBP results. The results in neighboring time window can also be averaged to present a more continuous results. But all these post-processings are totally optional.
