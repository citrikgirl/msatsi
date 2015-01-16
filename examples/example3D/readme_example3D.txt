Description of example 3:
Synthetic test varying the parameters of dip_direction, dip_angle and rake (3D)

The focal mechanisms are divided into 75 grid points that vary over X, Y and Z coordinates (3x5x5). Therefore, the problem is 3D. Each grid has 100 focal mechanisms.
**************************************************

Input data is saved in file 'INPUT_example3D.mat' as well as in 'INPUT_example3D.txt'  (same content)

The corresponding MATLAB path should contain at least msatsi.m and the executable files from corresponding platform
(i.e. satsi4d_tradeoff.exe,satsi4d.exe,bootmech4d.exe,bootuncert.exe)

First load the data from 'INPUT_EXAMPLE3D.mat' (or 'INPUT_example3D.txt') into workspace

***************************************************
MSATSI INVERSION:

To obtain results shown in the paper, type:
load('INPUT_EXAMPLE3D.mat');
[OUT] = msatsi('example3',example3D,'BootstrapResamplings',2000,'Damping','off')

****************************************************
MSATSI_PLOT VISUALIZATION:

By using the routine msatsi_plot, we can obtain different pictures.

1# By typing  [Hs, Hr] = msatsi_plot('example3','stereovolume','Xlabel',{'dip direction'},'YLabel',{'dip angle'},'ZLabel',{'Rake'})
   We obtain two pictures(one for stress and one for R values)
   Examples of these are saved as 'example3_stereovolume_S.png' and 'example3_stereovolume_R.png'
   This command was used to prepare the picture shown in the paper.

2# By typing  [Hs, Hr] = msatsi_plot('example3','stereovolume','ConfidenceIntervals','Bootstraps')
   We obtain two pictures. The picture for stress shows this time the uncertainties as bootstrap resampling. R plot keeps the same as in the case before.
   Examples of the stress plot is saved as: 'example3_stereovolume_bootstraps_S.png'

3# One useful property for 3D is 'Slice'. In this example, we will represent only one 2D slice of our 3D results.
   For example, if we are interested in the inversions were dip direction = 45°, the corresponding X coordinate is 1 (because it is the second grid point). Therefore, we have to specify the condition 'X==1'.
   By typing:  [Hs, Hr] = msatsi_plot('example3','stereomap','Slice','X==1')
   We obtain the stereomap (2D) plot along the specified parameter. Two plots are created, one for stress and one for R values.
   Examples of this plots are saved as: 'example3_stereomap_slice_X=1_S.png' and example3_stereomap_slice_X=1_R.png