Description of example 0D (not included in the paper):

Set of 150 focal mechanisms describing normal faulting regime grouped in a single grid point [0 0]. Therefore, the problem is 0D.
**************************************************

Input data is saved in file 'INPUT_example0D.mat' as well as in 'INPUT_example0D.txt'  (same content)

The corresponding MATLAB path should contain at least msatsi.m and the executable files from corresponding platform
(i.e. satsi2d_tradeoff.exe,satsi2d.exe,bootmech2d.exe,bootuncert.exe)

First load the data from 'INPUT_EXAMPLE0D.mat' (or 'INPUT_example0D.txt') into workspace

***************************************************
MSATSI INVERSION:

To run the inversion for a single grid point, we can type:
load('INPUT_EXAMPLE0D.mat');
[OUT] = msatsi('example0',example0D,'PTPlots','on','BootstrapResamplings',2000)

Since it is only one grid point, a damped inversion cannot be performed (there are no neighboring grid points). Therefore, the 'Damping' is automatically settled to 'off' (no need to write it). As a consequence, there is no tradeoff curve output.

In this case, we have activated the P and T axis for this grid point. The plot is contained in the folder /example0.

****************************************************
MSATSI_PLOT VISUALIZATION:

By using the routine msatsi_plot, we can obtain different pictures.

1# By typing  [Hs, Hr] = msatsi_plot('example0','stereonet')
   We obtain two pictures (one for stress and one for R values)
   Examples of these are saved as 'example0_stereonet_S.png' and 'example0_stereonet_R.png'
   
2# By typing [Hs, Hr] = msatsi_plot('example0','stereonet','ConfidenceIntervals','bootstraps')
   We obtain two analogous pictures as the case before but this time the uncertainties are plotted as colored dots directly from the bootstrap resampling.
   Examples of these are saved as 'example0_stereonet_bootstraps_S.png' and 'example0_stereonet_bootstraps_R.png'
   
