Description of example 2:
Surface stress field orientation at The Geysers geothermal reservoir

The focal mechanisms are divided into 20 grid points that vary over the X and Y coordinate. Therefore, the problem is 2D. Each grid has a different n° of focal mechanisms (from 35 to 448).
**************************************************

Input data is saved in file 'INPUT_example2D.mat' as well as in 'INPUT_example2D.txt'  (same content)

The corresponding MATLAB path should contain at least msatsi.m and the executable files from corresponding platform
(i.e. satsi2d_tradeoff.exe,satsi2d.exe,bootmech2d.exe,bootuncert.exe)

First load the data from 'INPUT_EXAMPLE2D.mat' (or 'INPUT_example2D.txt') into workspace

***************************************************
MSATSI INVERSION:

To obtain the results of the surface distribution shown in the paper, type:
load('INPUT_EXAMPLE2D.mat');
[OUT] = msatsi('example2',example2D,'BootstrapResamplings',2000,'MinEventsNode',30)

****************************************************
MSATSI_PLOT VISUALIZATION:

By using the routine msatsi_plot, we can obtain different pictures.

1# By typing just [Hs, Hr] = msatsi_plot('example2','stereomap')
   We obtain only two pictures regardless of grid point number. The pictures contain one stereonet for each grid point in a map view.
   Examples of these are saved as: 'example2_stereomap_intervals_S.png' and 'example2_stereomap_intervals_R.png'

2# To create the figure from the example shown in the paper we can first define the grid by the cells:
   xgridlab = num2cell([511:2:525]);
   ygridlab = num2cell([4289:2:4299]);
   Then, we type: [Hs, Hr] = msatsi_plot('example2','stereomap','ConfidenceIntervals','Bootstraps','XgridLabel',xgridlab,'YGridLabel',ygridlab)
   We obtain two pictures regardless of grid point number. (one for stress and one for R values)
   Examples of these are saved as 'example2_stereomap_bootstrap_S.png' and 'example2_stereomap_bootstrap_R.png'

3# To place the stress inversion in any arbitrary position (e.g. their coordinates), we first create matrix with X' Y' (the new values to the grids ordered as OUT.GRID)
   arb_grid = [(20:-1:1)',(5:1:24)'];
   Then, we type:  [Hs, Hr] = msatsi_plot('example2','stereomap','ArbitraryGrid',arb_grid,'ScaleFactor',0.5)
   By playing with the value of scale factor we will achieve increase/decrease size of circles and/or overlapping of them
   Examples of these pictures are saved as 'example2_stereomap_arbitrary_grid_S.png' 'and example2_stereomap_arbitrary_grid_R.png'
   
4# By typing: [Hs, Hr] = msatsi_plot('example2','wsm')
   We obtain the map with the direction of maximum horizontal stresses for each grid. This is the second picture shown in Fig 6b.
   Example of this is saved as: 'example2_wsm.png'
   
5# By typing: [Hs, Hr] = msatsi_plot('example2','wsm','IntermediateStresses','on')
   We obtain the same picture as the previous example but this time also the intermediate stress inversion results (i.e. those not included as NF,TF or SS).
   Example of this is saved as: 'example2_wsm.png'
   