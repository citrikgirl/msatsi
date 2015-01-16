Description of example 1:
Temporal stress field orientation changes at NAFZ (co-seismic period)
This data is part of the studies from Bohnhoff et al., 2006,GJI  and Ickrath et al.,2013, GJI. Here, it is used as example to ilustrate the method.

The focal mechanisms are divided into 27 grid points that vary over the Y coordinate. Therefore, the problem is 1D. Each grid has 70 focal mechanisms.
**************************************************

Input data is saved in file 'INPUT_example1D.mat' as well as in 'INPUT_example1D.txt'  (same content)

The corresponding MATLAB path should contain at least msatsi.m and the executable files from corresponding platform
(i.e. satsi2d_tradeoff.exe,satsi2d.exe,bootmech2d.exe,bootuncert.exe)

First load the data from 'INPUT_EXAMPLE1D.mat' (or 'INPUT_example1D.txt') into workspace

***************************************************
MSATSI INVERSION:

To obtain results from the coseismic period shown in the paper, type:
load('INPUT_EXAMPLE1D.mat');
[OUT] = msatsi('example1',example1D,'BootstrapResamplings',2000)

****************************************************
MSATSI_PLOT VISUALIZATION:

By using the routine msatsi_plot, we can obtain different pictures.

1# By typing  [Hs, Hr] = msatsi_plot('example1','stereonet','ConfidenceIntervals','bootstraps','SText','on')
   We obtain two pictures for each grid point (one for stress and one for R values)
   Examples of these are saved as 'example1_stereonet_bootstrap_S.png' and 'example1_stereonet_bootstrap_R.png'

2# By typing [Hs, Hr] = msatsi_plot('example1','stereonet','ConfidenceIntervals','off')
   We obtain only two pictures regardless of grid point number. The pictures contain the trajectory described by the best solution along the grid points.
   Examples of these are saved as: 'example1_stereonet_off_S.png' and 'example1_stereonet_off_R.png'

3# By typing only: [Hs, Hr] = msatsi_plot('example1','stereonet')
   We would obtain the same previous pictures but this time together with the uncertainties in each case.

4# By typing: [Hs, Hr] = msatsi_plot('example1','profile','Title',{'This is the example n°1'})
   We obtain three pictures: one for variations with trend for s1,s2 and s3 along the grid points, second analogous picture for variations in plunge, and a third one with variations in R. Uncertainties are represented with error bars.With the property 'Title', we have written the title we wished on top of each figure. 
   Examples of these are saved as: 'example1_profile_TR.png', 'example1_profile_PL.png', 'example1_profile_R.png'

5# By typing: [Hs, Hr] = msatsi_plot('example1','profile','ConfidenceIntervals','Bootstraps')
   We obtain analogous pictures as the case before, this time with the uncertainties plotted as points from the bootstrap resampling. This example was the one shown in the example from the paper.