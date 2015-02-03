
% Script to test msatsi
% 0D
file0 = 'examples\example0D\INPUT_example0D.mat';
load(file0);
[OUT] = msatsi('test0',example0D,'BootstrapResamplings', 20,'PTPlots', 'on');

% 1D
file1 = 'examples\example1D\INPUT_example1D.mat';
load(file1);
[OUT] = msatsi('test1',example1D,'BootstrapResamplings', 20,'PTPlots', 'off');

% 2D
file2 = 'examples\example2D\INPUT_example2D.mat';
load(file2);
[OUT] = msatsi('test2',example2D,'BootstrapResamplings', 20,'PTPlots', 'off');

% 3D
file3 = 'examples\example3D\INPUT_example3D.mat';
load(file3);
[OUT] = msatsi('test3',example3D,'BootstrapResamplings', 20,'PTPlots', 'off');
