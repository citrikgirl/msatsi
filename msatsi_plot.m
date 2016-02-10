function [Hs, Hr] = msatsi_plot(projectname, type, varargin)
%MSATSI_PLOT  Plot output of MSATSI program.
%   Auxillary routine allowing to display stress inversion results using a
%   variety of 0D-3D plots. For the full documentation of msatsi_plot.m routine,
%   see http://induced.pl/msatsi website or look into the documentation files
%   provided with MSATSI package.
%
%   msatsi_plot(projectname, type) displays results of 'projectname' stress
%     inversion using plot type 'type'. 'projectname' is either the name of a
%     directory where project files are stored (this directory is created by
%     calling msatsi.m program) or alternatively a structure array obtained
%     as an output of msatsi.m routing. The 'type' variable is a string
%     denoting the type of plot. Valid types of plot are 'stereonet',
%     'stereomap', 'profile', 'stereovolume', and 'wsm'.
%
%     'stereonet' allows displaying 0D-1D data using stereographic (Schmidt or
%       Wullf) projection. In case of 1D data, the trace of stress rotation is
%       plotted.
%
%     'stereomap' is a plot for 2D data only allowing to show the spatial
%       distribution of stress orientatiom.
%
%     'profile' produces plot(s) with changes in stress field orientation
%       with consequtive grid points.
%
%     'wsm' generates spatial distribution of 1D/2D output data in World
%       Stress Map fashion.
%
%
%   msatsi_plot(projectname, type, 'PropertyName', PropertyValue, ...)
%     Allows to specify/modify the detailes of plots created. See documentation
%     of MSATSI package for details.
%
%   See also MSATSI

%   Copyright 2013-2014 Patricia Martínez-Garzón <patricia@gfz-potsdam.de>
%                       Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%   $Revision: 1.0.8 $  $Date: 2015.02.04 $

%   Development info:
%     1.0.8 SilentMode and SaveImage options added. Small corrections to existing code.
%           Correction to the best solution of satsi2d.

global scale_factor;
global stereo_projection;
global grid_color;
global grid_step_azimuth;
global grid_step_plunge;
global Out;
global AxesColor;
global AxesEnabled;
global AXESI;
global confidence_intervals;
global silent_mode;
global save_image;
global project_name;

project_name = projectname;

% If Out is character array, read satsi.mat from the specified
% directory and put it into the Out structure array.
if ischar(projectname)
  load([projectname '/' projectname '_OUT.mat']);
  Out = OUT;  clear OUT;
end

% Input parser section
p = inputParser;
p.addParamValue('ConfidenceIntervals', 'intervals', @(x)any(strcmpi(x,{'intervals','off','bootstraps'})));
p.addParamValue('Slice', '', @(x) ischar(x));
p.addParamValue('Grid', 'on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('GridColor', [0.5 0.5 0.5]); % [0.8 0.8 0.8]
p.addParamValue('GridStepAzimuth', 15);
p.addParamValue('GridStepPlunge', 15);
p.addRequired('Out', @(x) isstruct(x));
p.addParamValue('Projection', 'wulff', @(x)any(strcmpi(x,{'schmidt','wulff'})));
p.addParamValue('ScaleFactor', 0.5, @(x) isscalar(x) && x > 0);
p.addParamValue('Stereonet', 'on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('S1Color', [1 0 0]);  %  [0.2 0.2 0.2]
p.addParamValue('S2Color', [0 1 0]);  %  [0.4 0.4 0.4]
p.addParamValue('S3Color', [0 0 1]); %   [0.6 0.6 0.6]
p.addParamValue('XGridLabel', {}, @(x)iscell(x));
p.addParamValue('YGridLabel', {}, @(x)iscell(x));
p.addParamValue('ZGridLabel', {}, @(x)iscell(x));
p.addParamValue('XLabel', {}, @(x)iscell(x));
p.addParamValue('YLabel', {}, @(x)iscell(x));
p.addParamValue('ZLabel', {}, @(x)iscell(x));
p.addParamValue('S1PatchColor',[1.0 0.7 0.7]); %  [0.2 0.2 0.2]
p.addParamValue('S2PatchColor',[0.7 1.0 0.7]); %  [0.4 0.4 0.4]
p.addParamValue('S3PatchColor',[0.7 0.7 1.0]); %  [0.6 0.6 0.6]
p.addParamValue('ArbitraryGrid',[],@(x)ismatrix(x));
p.addParamValue('S1','on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('S2','on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('S3','on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('IntermediateStresses','off', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('Symbol', '+', @(x)any(strcmpi(x,{'+','o','*','.','x','s','d','^','v','<','>','p','h','none'})));
p.addParamValue('SText', 'off', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('Title', 'off');
p.addParamValue('RPlot', 'on', @(x)any(strcmpi(x,{'off','on'})));
p.addParamValue('SilentMode', 'off', @(x)any(strcmpi(x,{'off','on'})));
p.addParamValue('SaveImage', 'off', @(x)any(strcmpi(x,{'off','on'})));
p.addRequired('type', @(x)any(strcmpi(x, {'stereonet','stereomap', 'profile', 'stereovolume','wsm'})));

% Parse input parameters.
p.parse(Out, type, varargin{:});

% Get input parameters.
% dimension = p.Results.Dimension;
scale_factor = p.Results.ScaleFactor;
view_stereonet = p.Results.Stereonet;
stereo_projection = p.Results.Projection;
view_grid = p.Results.Grid;
grid_color = p.Results.GridColor;
grid_step_azimuth = p.Results.GridStepAzimuth;
grid_step_plunge = p.Results.GridStepPlunge;
s1color = p.Results.S1Color; s2color = p.Results.S2Color; s3color = p.Results.S3Color;
s1 = p.Results.S1; s2 = p.Results.S2; s3 = p.Results.S3;
s1_patch_color = p.Results.S1PatchColor; s2_patch_color = p.Results.S2PatchColor; s3_patch_color = p.Results.S3PatchColor;
intermediate = p.Results.IntermediateStresses;
confidence_intervals = p.Results.ConfidenceIntervals;
symbol = p.Results.Symbol;
stext = p.Results.SText;
title_text = p.Results.Title;
R_plot = p.Results.RPlot;
cross_section = p.Results.Slice;
xgridlabel = p.Results.XGridLabel;
ygridlabel = p.Results.YGridLabel;
zgridlabel = p.Results.ZGridLabel;
x_label = p.Results.XLabel;
y_label = p.Results.YLabel;
z_label = p.Results.ZLabel;
silent_mode = p.Results.SilentMode; if strcmpi(silent_mode,'off'), silent_mode = false; else silent_mode = true; end
save_image = p.Results.SaveImage; if strcmpi(save_image,'off'), save_image = false; else save_image = true; end
GRID_ARB = p.Results.ArbitraryGrid;

%---- Determine colors and whether or not to plot stress axes.
AxesColor = {s1color, s2color, s3color};
AxesEnabled = {true, true, true};
AXESI = [];
naxes = 3;
if strcmpi(s1,'off'), AxesEnabled{1} = false; naxes = naxes - 1; else AXESI = [AXESI 1]; end;
if strcmpi(s2,'off'), AxesEnabled{2} = false; naxes = naxes - 1; else AXESI = [AXESI 2]; end;
if strcmpi(s3,'off'), AxesEnabled{3} = false; naxes = naxes - 1; else AXESI = [AXESI 3]; end;
if naxes == 0
  error('No stress axes is enabled, plotting aborted.');
end
PATCHCOLORS = {s1_patch_color, s2_patch_color, s3_patch_color};
Out.GRID = Out.GRID(:,1:end-1);

% Cut the dimension if necessary.
if ~isempty(cross_section)
  grd_str = strrep(cross_section,'X','*(:,1)');
  grd_str = strrep(grd_str,'Y','*(:,2)');
  grd_str = strrep(grd_str,'Z','*(:,3)');
  grd_str = strrep(grd_str,'T','*(:,4)');
  grd_str = strrep(grd_str,'*','Out.GRID');
  
  slb_str = strrep(cross_section,'X','*(:,1)');
  slb_str = strrep(slb_str,'Y','*(:,2)');
  slb_str = strrep(slb_str,'Z','*(:,3)');
  slb_str = strrep(slb_str,'T','*(:,4)');
  slb_str = strrep(slb_str,'*','Out.SLBOOT_TENSOR');
  
  sle_str = strrep(cross_section,'X','*(:,1)');
  sle_str = strrep(sle_str,'Y','*(:,2)');
  sle_str = strrep(sle_str,'Z','*(:,3)');
  sle_str = strrep(sle_str,'T','*(:,4)');
  sle_str = strrep(sle_str,'*','Out.BOOTST_EXT');
  
  % Determine which rows to preserve.
  I_GRD = eval(grd_str);
  I_SLB = eval(slb_str);
  I_SLE = eval(sle_str);
  
  % Cut matrices
  Out.SLBOOT_TENSOR = Out.SLBOOT_TENSOR(I_SLB,:);
  Out.SLBOOT_TRPL = Out.SLBOOT_TRPL(I_SLB,:);
  Out.SUMMARY_TABLE = Out.SUMMARY_TABLE(I_GRD,:);
  Out.BOOTST_EXT = Out.BOOTST_EXT(I_SLE,:);
  Out.GRID = Out.GRID(I_GRD,:);
end

% Total number of dimension (2 or 4)
GRID = Out.GRID;
ndim_total = size(GRID,2);

% Effective number of dimensions and number of unique grids in each dim.
N_GRID_UNIQUE = zeros(1,size(GRID,2));
ndim_eff = 0;   % effective number of dimensions (0-4)
COLGRID_EFF = [];
for i=1:size(GRID,2)
  N_GRID_UNIQUE(i) = length(unique(GRID(:,i)));
  if  N_GRID_UNIQUE(i) > 1
    ndim_eff = ndim_eff + 1;
    COLGRID_EFF = [COLGRID_EFF; i];  %#ok<AGROW>
  end
end
COLGRID_EFF_old = COLGRID_EFF; % COLGRID_EFF_old is used for EXTENDED SLBOOT DATA MATRIX ONLY.

if ndim_total == 4 && ndim_eff == 2
  % We have to cut from 3D to 2D
  ndim_total = 2;
  Out.SLBOOT_TENSOR = Out.SLBOOT_TENSOR(:,[COLGRID_EFF' 5:end]);
  Out.SLBOOT_TRPL = Out.SLBOOT_TRPL(:,[COLGRID_EFF' 5:end]);
  Out.GRID = Out.GRID(:,COLGRID_EFF');
  GRID = Out.GRID;
  COLGRID_EFF_old = COLGRID_EFF;
  COLGRID_EFF = [1 2];
end

GridLabel = {'x','y','z','t'};

%==== Preparation of common data useful in plotting various stress char. =
% TODO: Order of grid points is important for some plots!
% Prepare trend & plunge of stresses depending on the type of data to read

TR_BEST = round(Out.SUMMARY_TABLE(:,[4,10,16])); % We round quatities to optimize the bootstrap plots
PL_BEST = round(Out.SUMMARY_TABLE(:,[7,13,19]));
PHI_BEST = Out.SUMMARY_TABLE(:,1);
R_BEST = 1 - PHI_BEST;
TKO_BEST = (pi/2 - PL_BEST * pi/180);
AZM_BEST = TR_BEST * pi/180;

if strcmpi(confidence_intervals,'intervals') || strcmpi(confidence_intervals,'off')
  
  TR_MIN = round(Out.SUMMARY_TABLE(:,[4,10,16]+1));
  PL_MIN = round(Out.SUMMARY_TABLE(:,[7,13,19]+1));
  PHI_MIN = Out.SUMMARY_TABLE(:,2); R_MAX = 1 - PHI_MIN;
  TKO_MIN = (pi/2 - PL_MIN * pi/180);
  AZM_MIN = TR_MIN * pi/180;
  
  TR_MAX = round(Out.SUMMARY_TABLE(:,[4,10,16]+2));
  PL_MAX = round(Out.SUMMARY_TABLE(:,[7,13,19]+2));
  PHI_MAX = Out.SUMMARY_TABLE(:,3); R_MIN = 1 - PHI_MAX;
  TKO_MAX = (pi/2 - PL_MAX * pi/180);
  AZM_MAX = TR_MAX * pi/180;
end

if strcmpi(confidence_intervals,'bootstraps') || strcmpi(type, 'wsm')
  TR_BOOT  = cell(size(GRID,1),1);
  PL_BOOT  = cell(size(GRID,1),1);
  PHI_BOOT = cell(size(GRID,1),1);
  R_BOOT   = cell(size(GRID,1),1);
  SBOOT    = cell(size(GRID,1),1);
  
  % Base on the number of effective dimensions we have to select the
  % appropriate samples from BOOTST_EXT for each grid point.
  if ndim_eff <= 2
    for i=1:size(GRID,1)
      if ndim_eff == 0
        II = Out.BOOTST_EXT(:,1) > -1; % dummy case
      elseif ndim_eff == 1
        x = GRID(i,COLGRID_EFF(1));
        II = Out.BOOTST_EXT(:,COLGRID_EFF_old(1)) == x;
      elseif ndim_eff == 2
        x = GRID(i,COLGRID_EFF(1));
        y = GRID(i,COLGRID_EFF(2));
        II = Out.BOOTST_EXT(:,COLGRID_EFF_old(1)) == x ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(2)) == y;
      end
      % TR and PL are rounded to optimize the bootstrap plot for
      % many points.
      SBOOT{i} = [Out.BOOTST_EXT(II,5) round(Out.BOOTST_EXT(II,6:end))];
    end
  else % 3D, 4D
    for i=1:size(GRID,1)
      x = GRID(i,COLGRID_EFF(1));
      y = GRID(i,COLGRID_EFF(2));
      z = GRID(i,COLGRID_EFF(3));
      if ndim_eff == 3
        II = Out.BOOTST_EXT(:,COLGRID_EFF_old(1)) == x ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(2)) == y ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(3)) == z;
      else
        t = GRID(i,COLGRID_EFF(4));
        II = Out.BOOTST_EXT(:,COLGRID_EFF_old(1)) == x ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(2)) == y ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(3)) == z ...
          & Out.BOOTST_EXT(:,COLGRID_EFF_old(4)) == t;
      end
      % TR and PL are rounded to optimize the bootstrap plot for
      % many points.
      SBOOT{i} = [Out.BOOTST_EXT(II,5) round(Out.BOOTST_EXT(II,6:end))];
    end
  end
  
  for i=1:size(GRID,1)
    PHI_BOOT{i,1} = SBOOT{i,1}(:,1);
    R_BOOT{i,1}   = 1 - SBOOT{i,1}(:,1);
    TR_BOOT{i,1}  = SBOOT{i,1}(:,2:2:6);
    PL_BOOT{i,1}  = SBOOT{i,1}(:,3:2:7);
  end
end

% Finally, replace if necessary the GRID & Out.GRID with arbitrary points.
if ~isempty(GRID_ARB)
  if size(GRID_ARB) == size(GRID)
    GRID = GRID_ARB;
    Out.GRID = GRID_ARB;
  else
    warning('MSATSI_PLOT:MatrixNoMatch',['Arbitrary grid must match the original grid size [' num2str(size(GRID)) ']']);
  end
end

% Prepare the grid points.
if ndim_total == 2
  if ndim_eff == 2
    X0 = GRID(:,COLGRID_EFF(1));
    Y0 = GRID(:,COLGRID_EFF(2));
  elseif ndim_eff == 1
    X0 = GRID(:,COLGRID_EFF(1));
    Y0 = zeros(length(X0),1);
    if sum(X0) == 0 % I guess it does not make a lot of sense.
      Y0 = GRID(:,COLGRID_EFF(2));
    end
  elseif ndim_eff == 0;
    X0 = 0;
    Y0 = 0;
  end
end




%==== Plot generation routines ====================================================================

% Generate output handle structure.
Hs = []; % Handles of stress-related plots.
Hr = []; % Handles of r-related plots.

%==== A) Stress orientation with 'STEREONET' 0D, 1D case ==========================================
if strcmpi(type, 'stereonet')
  
  %---- A.1: Plot by data points from bootstrap resampling
  if strcmpi(confidence_intervals,'bootstraps')
    [XX_BEST,YY_BEST] = project(AZM_BEST, TKO_BEST);
    [dummy1,dummy2,Hs] = make_bootstrap_dots(type,GRID,SBOOT,XX_BEST,YY_BEST,X0,Y0, ...
      zeros(length(X0),1),view_stereonet, view_grid,stext, ...
      title_text,symbol,Hs); %#ok<ASGLU>
    
    % Generate R plot.
    if strcmpi(R_plot,'on')
      for i=1:size(GRID,1)
        R = 1 - SBOOT{i}(:,1);
        X = 0:0.1:1;
        Hr = newrfigure(Hr);
        [n,xout] = hist(R,X);
        h = bar(xout,n);
        set(h,'FaceColor','w','EdgeColor','k')
        best_value = (round(R_BEST(i)*100))/100;
        text(0.75,max(n)-0.5,['R best = ' num2str(best_value)],'Fontsize',12)
        xlim([0 1]);
        grid on; box on;
        xlabel('R value');
        ylabel('Number of cases');
        % **Save figure
        saveimage('-stereonet-R',GRID,i);
      end
    end
    
  %---- A.2 Representation using data points from interval patches
  elseif strcmpi(confidence_intervals, 'intervals') || strcmpi(confidence_intervals, 'off')
    
    if ndim_eff == 0 || ndim_eff == 1
      
      [XX_BEST, YY_BEST] = project(AZM_BEST, TKO_BEST);
      
      %---- Plot stereonet with trace of stress directions
      %Hs = newsfigure(Hs);
      % PUT THIS IN ORDER AGAIN!
      Hs = gcf;
      % plot_stereonet(0, 0, view_stereonet, view_grid) UNCOMMENT THIS
      hold on;
      makepatch(type,TR_MIN,TR_MAX,PL_MIN,PL_MAX,X0,Y0,PATCHCOLORS, ...
        AxesEnabled,confidence_intervals)
      for i = AXESI
        line(XX_BEST(:,i)',YY_BEST(:,i)','Color', AxesColor{i}, 'LineWidth', 1);
        plot(XX_BEST(:,i),YY_BEST(:,i),symbol,'Color',AxesColor{i},'MarkerSize',8);
      end
      hold off;
      axis equal; box on; axis off;
      if strcmpi(stext,'on')
        text(XX_BEST(1,1) - 0.05,YY_BEST(1,1) + 0.05,'\sigma_1','FontSize',20,'Color',AxesColor{1});
        text(XX_BEST(1,2) - 0.05,YY_BEST(1,2) + 0.05,'\sigma_2','FontSize',20,'Color',AxesColor{2});
        text(XX_BEST(1,3) - 0.05,YY_BEST(1,3) + 0.05,'\sigma_3','FontSize',20,'Color',AxesColor{3});
      end
      plot_title(title_text,Out.Caption);
      %---- END: Plot stereonet with trace of stress directions
      % **Save figure
      saveimage('-stereonet');
    else
      error(['You selected ''stereonet'' plot type, but the effective ' ...
        'dimension of data grid is not 0D or 1D (ndim_eff = ' num2str(ndim_eff) '). '  ' ']);
    end
    %plot_title(title_text,Out.Caption); % Moved up to if block above.
    
    %---- R plot.
    if strcmpi(R_plot,'on')
      Hr = newrfigure(Hr);
      hold on;
      for k = 1:size(GRID,1)
        if strcmpi(confidence_intervals,'off') == 0
          xline = R_MIN(k):0.01:R_MAX(k);
          yline = ones(length(xline),1)*k;
          d2 = line(yline,xline);
          set(d2,'Color','k','LineWidth',2);
        end
      end
      d1 = plot(1:k,R_BEST,'sr');
      set(d1,'MarkerSize',5,'MarkerFaceColor','r');
      grid on; box on;
      xlabel('X Y');
      ylabel('R values');
      set(gca,'XTick',1:k);
      set(gca,'XTickLabel',{num2str(GRID)});
      %---- R plot.
      % **Save figure
      saveimage('-stereonet-R',GRID,i);
    end
  end
  
  
  
%==== B) Plot data as "PROFILE" ===================================================================
elseif strcmpi(type, 'profile')
  
  if ndim_eff == 1 % data must be 1D
    
    XGRID = Out.GRID(:, COLGRID_EFF(1)); % TODO:  modify in arbitrary dimension
    % GK: The following two lines should be realized inside
    % profile_correct.
    if strcmpi(confidence_intervals,'intervals'); TR_BOOT = 0; end % dummy
    if strcmpi(confidence_intervals,'bootstraps'); TR_MIN = 0; TR_MAX = 0; PL_MIN = 0; PL_MAX = 0; end % dummy
    
    % Patricia's optimization scheme.
    %[TR_BEST,TR_BOOT,PL_BEST,TR_MIN,TR_MAX,PL_MIN,PL_MAX] = profile_correct(TR_BEST,TR_BOOT,PL_BEST,TR_MIN,TR_MAX,PL_MIN,PL_MAX);
    
    % Grzegorz's optimization scheme.
    [TR_BEST,TR_BOOT,TR_MIN,TR_MAX, JE, JD, TR_MAX_E, TR_MIN_E] = ...
      profile_correct2(TR_BEST,TR_BOOT,TR_MIN,TR_MAX);
    
    %---- B.1 Plot: Generate profile plot of Sx trends.
    Hs = newsfigure(Hs);
    k=0;
    for i=AXESI
      k = k + 1;
      subplot(naxes,1,k);
      if k == 1
        plot_title(title_text,Out.Caption);
      end
      hold on;
      if strcmpi(confidence_intervals,'intervals')
        errorbar(XGRID,TR_BEST(:,i),TR_BEST(:,i) - TR_MIN(:,i),TR_MAX(:,i) - TR_BEST(:,i),'--','Color',AxesColor{i});
        
        % Grzegorz's optimization scheme >>>
        if sum(JE(:,i)>0)
          errorbar(XGRID(JE(:,i)),zeros(size(XGRID(JE(:,i)))),zeros(size(XGRID(JE(:,i)))),TR_MAX_E(JE(:,i),i),'Color',AxesColor{i});
        end
        if sum(JD(:,i)>0)
          errorbar(XGRID(JD(:,i)),180*ones(size(XGRID(JD(:,i)))),TR_MIN_E(JD(:,i),i),zeros(size(XGRID(JD(:,i)))),'Color',AxesColor{i});
        end
        % <<< Grzegorz's optimization scheme
      elseif strcmpi(confidence_intervals,'bootstraps')
        for j = 1:size(GRID,1)
          TROPT = unique(TR_BOOT{j}(:,i));
          plot(ones(length(TROPT),1)*XGRID(j),TROPT,'LineStyle','none','Marker','.','Color',[0.6 0.6 0.6],'MarkerSize',6);
        end
      end
      plot(XGRID,TR_BEST(:,i),'-','Marker','.','Color',AxesColor{i});
      hold off;
      grid on; box on;
      xlim([min(XGRID) max(XGRID)]);
      ylim([0 180]) % for micheles plot
      set(gca,'YTick',0:45:180);
      ylabel('Trend');
      if k==naxes
        xlabel(['Grid # (' GridLabel{COLGRID_EFF} ' coordinate)']); % TODO: must be modified in arbitrary dimension
      end
    end
    %---- END: Generate profile plot of Sx trends.
    % **Save figure
    saveimage('-profile-trends');
    
    %---- B.2 Plot: Generate profile plot of Sx plunges.
    Hs = newsfigure(Hs);
    k=0;
    for i=AXESI
      k = k + 1;
      subplot(naxes,1,k);
      if k == 1
        plot_title(title_text,Out.Caption);
      end
      hold on;
      if strcmpi(confidence_intervals,'intervals')
        errorbar(XGRID,PL_BEST(:,i),PL_BEST(:,i) - PL_MIN(:,i),PL_MAX(:,i) - PL_BEST(:,i),'--','Color',AxesColor{i});
      elseif strcmpi(confidence_intervals,'bootstraps')
        for j = 1:size(GRID,1)
          PLOPT = unique(PL_BOOT{j}(:,i));
          plot(ones(length(PLOPT),1)*XGRID(j),PLOPT,'LineStyle','none','Marker','.','Color',[0.6 0.6 0.6],'MarkerSize',6);
        end
      end
      plot(XGRID,PL_BEST(:,i),'-','Marker','.','Color',AxesColor{i});
      hold off;
      grid on; box on;
      ylabel('Plunge');
      xlim([min(XGRID) max(XGRID)]);
      ylim([0 90]) % for micheles plot
      set(gca,'YTick',0:30:90);
      if k==naxes
        xlabel(['Grid # (' GridLabel{COLGRID_EFF} ' coordinate)']); % TODO: modify in arbitrary dimension
      end
    end
    %---- END: Generate profile plot of Sx plunges.
    % **Save figure
    saveimage('-profile-plunges');
    
    %---- B.3: Generate R plot.
    if strcmpi(R_plot,'on')
      Hr = newrfigure(Hr);
      hold on;
      if strcmpi(confidence_intervals,'intervals')
        errorbar(XGRID,R_BEST(:,1),R_BEST(:,1) - R_MIN(:,1),R_MAX(:,1) - R_BEST(:,1),'--','Color','k');
      elseif strcmpi(confidence_intervals,'bootstraps')
        for j = 1:size(GRID,1)
          plot(ones(length(R_BOOT{j}(:,1)),1)*XGRID(j),R_BOOT{j}(:,1),'LineStyle','none','Marker','.','Color',[0.6 0.6 0.6],'MarkerSize',6);
        end
      end
      plot(XGRID,R_BEST,'-','Marker','.','Color','k');
      plot_title(title_text,Out.Caption);
      hold off; grid on; box on;
      xlabel(['Grid # (' GridLabel{COLGRID_EFF} ' coordinate)']);
      ylabel('R values');
    end
    %---- END: Generate R plot.
    % **Save figure
    saveimage('-profile-R');
  else
    error(['You selected ''profile'' plot type, but the effective ' ...
      'dimension of data grid is not 1D (ndim_eff = ' num2str(ndim_eff) '). ' '']);
  end
  
  
  
%=== C) Plot data as World Stress Map project =====================================================
elseif strcmpi(type, 'wsm')
  if ndim_total == 2 && ndim_eff == 2   % Plot valid ONLY for 2D
    
    % Lund&Townend extension
    [XBEST,YBEST,ZBEST] = sph2cart(TR_BEST*pi/180,PL_BEST*pi/180,ones(length(GRID),3));
    TR_BEST_SHmax = zeros(length(GRID),1);
    AZM_BEST_SHmax = zeros(length(GRID),1);
    for i = 1:length(GRID)
      alpha = SH([XBEST(i,1),YBEST(i,1),ZBEST(i,1)],[XBEST(i,2),YBEST(i,2),ZBEST(i,2)],[XBEST(i,3),YBEST(i,3),ZBEST(i,3)],R_BEST(i));
      AZM_BEST_SHmax(i) = alpha;
      TR_BEST_SHmax(i) = alpha*180/pi; % in degrees
    end
    % end extension
    
    % Determine SHMax axes for each grid (Zoback's Criteria, 1992)
    AxesSHMax = zeros(length(X0),1);  color = zeros(length(X0),3);
    for i = 1:length(X0)
      if PL_BEST(i,1) > 52 && PL_BEST(i,3) < 35 % NF
        AxesSHMax(i) = 2;  color(i,:) = [1 0 0];  % 'r'
      elseif PL_BEST(i,1) < 20 && PL_BEST(i,2) > 45 && PL_BEST(i,3) < 40 % SS
        AxesSHMax(i) = 1; color(i,:) = [0 1 0];  % 'g';
      elseif PL_BEST(i,1) < 35 && PL_BEST(i,3) > 52 % TF
        AxesSHMax(i) = 1; color(i,:) = [0 0 1];  % 'b';
      else
        AxesSHMax(i) = 0; color(i,:) = [0 0 0];  % 'k';
      end
    end
    
    %---- C.1: Generate WSM figure.
    if strcmpi(confidence_intervals,'intervals') || strcmpi(confidence_intervals,'off')
      Hs = newsfigure(Hs);
      hold on;
      TKO_BEST_SHmax = ones(size(GRID,1),1)*pi/2; % We rewrite TKO here as pi/2
      [XX_BEST_SHmax, YY_BEST_SHmax] = project(AZM_BEST_SHmax, TKO_BEST_SHmax);
    else
      error('Bootstraps are not available for ''wsm'' plot. Try ''intervals'' or ''off'' ');
    end
    
    if strcmpi(confidence_intervals,'off')
      X = XX_BEST_SHmax';  Y = YY_BEST_SHmax';
      if  strcmp(intermediate,'on'); I = AxesSHMax ~= 5;
      elseif strcmp(intermediate,'off'); I = AxesSHMax ~= 0; end
      if sum(I) ~= 0
        Xline = [X(I);zeros(size(X(I)))] + [X0(I) X0(I)]';
        Yline = [Y(I);zeros(size(Y(I)))] + [Y0(I) Y0(I)]';
        Xline2 = [-X(I);zeros(size(X(I)))] + [X0(I) X0(I)]';
        Yline2 = [-Y(I);zeros(size(Y(I)))] + [Y0(I) Y0(I)]';
        linecolor = color(I,:)';
        for k = 1:size(Yline,2);
          line(Xline(:,k),Yline(:,k),'Color',linecolor(:,k),'LineWidth',2,'LineStyle','-');
          line(Xline2(:,k),Yline2(:,k),'Color',linecolor(:,k),'LineWidth',2,'LineStyle','-');
        end
      end
      
    elseif strcmpi(confidence_intervals,'intervals')
      
      [XMIN,YMIN,ZMIN] = sph2cart(TR_MIN.*pi/180,PL_MIN*pi/180,ones(length(PL_MIN),3));
      [XMAX,YMAX,ZMAX] = sph2cart(TR_MAX.*pi/180,PL_MAX*pi/180,ones(length(PL_MAX),3));
      
      for j = 1:size(AxesSHMax,1)
        if   (AxesSHMax(j) ~=0) || (AxesSHMax(j) == 0 && strcmp(intermediate,'on'))
          %   extension lund and Townend
          alpha_min = SH([XMIN(j,1),YMIN(j,1),ZMIN(j,1)],[XMIN(j,2),YMIN(j,2),ZMIN(j,2)],[XMIN(j,3),YMIN(j,3),ZMIN(j,3)],R_MIN(j));
          alpha_max = SH([XMAX(j,1),YMAX(j,1),ZMAX(j,1)],[XMAX(j,2),YMAX(j,2),ZMAX(j,2)],[XMAX(j,3),YMAX(j,3),ZMAX(j,3)],R_MAX(j));
          tr_min_shmax = alpha_min * 180/pi;
          tr_max_shmax = alpha_max * 180/pi;
          %   end extension
          ARANGE = tr_min_shmax:tr_max_shmax; %#ok<NASGU>
          
          % Display bootstrap grid points.
          TKO = pi/2 - SBOOT{j}(:,3+(AxesSHMax(j)-1)*2) * pi / 180;
          % TKO = pi/2 - []
          
          [XBOOT,YBOOT,ZBOOT] = sph2cart([SBOOT{j}(:,2),SBOOT{j}(:,4),SBOOT{j}(:,6)].*pi/180,[SBOOT{j}(:,3),SBOOT{j}(:,5),SBOOT{j}(:,7)].*pi/180,ones(length(SBOOT{j}(:,3)),3));
          AZM_BOOT_SHMAX = zeros(length(SBOOT{j}(:,1)),1);
          TR_BOOT_SHMAX = zeros(length(SBOOT{j}(:,1)),1);
          for k = 1:length(SBOOT{j}(:,1))
            alpha_boot = SH([XBOOT(k,1),YBOOT(k,1),ZBOOT(k,1)],[XBOOT(k,2),YBOOT(k,2),ZBOOT(k,2)],[XBOOT(k,3),YBOOT(k,3),ZBOOT(k,3)],SBOOT{j}(k,1));
            AZM_BOOT_SHMAX(k) = alpha_boot;
            TR_BOOT_SHMAX(k) = alpha_boot*180/pi;
          end
          AZM = AZM_BOOT_SHMAX.*180/pi;  % replacement for corrected L&T
          
          TKO = [TKO; TKO]; %#ok<AGROW> % mirror
          AZM = [AZM; AZM+180]; %#ok<AGROW>
          
          [XX,YY] = project(AZM*pi/180, TKO);
          RR = sqrt(XX.^2+YY.^2);
          AZM = atan2(XX,YY)*180/pi+360;
          OO = AZM>=360; AZM(OO) = AZM(OO)-360;
          
          st = 2;
          ARANGE2 = 0:st:(360-st);
          kk = 1;
          DR = zeros(length(ARANGE2)-1,1);
          DR2 = zeros(length(ARANGE2)-1,1);
          DA = zeros(length(ARANGE2)-1,1);
          for a = ARANGE2
            J = AZM >= a & AZM < a+st;
            if sum(J)
              DR(kk) = max(TKO(J));
              DR2(kk) = max(RR(J));
              DA(kk) = a + st/2;
            end
            kk = kk + 1;
          end
          XX = DR2.*sin(DA*pi/180);
          YY = DR2.*cos(DA*pi/180);
          
          Fvc.VERTICES = [0 0; XX YY] + repmat([X0(j) Y0(j)],length(XX)+1,1);
          Fvc.FACES = [];
          for a = 1: (length(XX)-1)
            if DR(a) > 0 && DR(a+1) > 0
              Fvc.FACES = [Fvc.FACES; 1 a+1 a+2];
            end
          end
          patch(Fvc,'FaceColor',color(j,:),'EdgeColor','none')
          
        end
      end
    end
    
    set(gca,'XTick',0:max(X0),'YTick',0:max(Y0),'Color','none');
    hold off; axis equal; box on; grid on;
    
    plot_title(title_text,Out.Caption);
    plot_grids(xgridlabel,ygridlabel,zgridlabel);
    plot_labels(x_label,y_label,z_label);
    %---- END: C.1: Generate WSM figure.
    % **Save figure
    saveimage('-wsm');
    
  else
    error(['You selected ''wsm'' plot type, but the effective ' ...
      'dimension of data grid is not 2D (ndim_eff = ' num2str(ndim_eff) ')']);
  end
  
  
  
%==== D) Map with S1,2,3 on "STEREOMAP" 2D ========================================================
elseif strcmpi(type, 'stereomap')  %(2D)
  
  if ndim_total == 2
    %---- D.1: Stereomap plot.
    if strcmpi(confidence_intervals,'intervals') || strcmpi(confidence_intervals,'off')
      Hs = newsfigure(Hs);
      plot_stereonet(X0, Y0, view_stereonet, view_grid)
      hold on;
    end
    [XX_BEST,YY_BEST] = project(AZM_BEST, TKO_BEST);
    for i = AXESI
      if strcmpi(confidence_intervals,'intervals') || strcmpi(confidence_intervals,'off')
        X = XX_BEST(:,i)';
        Y = YY_BEST(:,i)';
        line([X;zeros(size(X))] + [X0 X0]',[Y;zeros(size(Y))] + [Y0 Y0]', 'Color',AxesColor{i}, 'LineWidth', 3);
        plot(XX_BEST(:,i)+X0,YY_BEST(:,i) + Y0,'.','Color',AxesColor{i});
      end
    end
    if strcmpi(confidence_intervals,'intervals')
      makepatch(type,TR_MIN,TR_MAX,PL_MIN,PL_MAX,X0,Y0,PATCHCOLORS,AxesEnabled,confidence_intervals)
      txt = {'\sigma_1','\sigma_2','\sigma_3'};
      if strcmpi(stext,'on')
        for n = 1: size(GRID,1)
          for m = 1:3
            if AxesEnabled{m}
              text(XX_BEST(n,m)  + X0(n) - 0.05,YY_BEST(n,m)  + Y0(n) + 0.05,txt{m},'FontSize',20,'Color',AxesColor{m});
            end
          end
        end
      end
      
    elseif strcmpi(confidence_intervals,'bootstraps')
      Z0 = zeros(length(X0),1);  % dummy variables to fill
      [dummy1,dummy2,Hs] = make_bootstrap_dots(type,GRID,SBOOT,XX_BEST,YY_BEST,X0,Y0,Z0,view_stereonet,view_grid,stext,title_text,symbol,Hs); %#ok<ASGLU>
    end
    plot_title(title_text,Out.Caption);
    plot_grids(xgridlabel,ygridlabel,zgridlabel);
    plot_labels(x_label,y_label,z_label);
    hold off; axis equal;  box on;
    if isempty(GRID_ARB)
      set(gca,'XTick',0:max(X0),'YTick',0:max(Y0),'Color','none');
    else
      set(gca,'Color','none');
    end
    %---- END: D.1: Stereomap plot.
    % **Save figure
    saveimage('-stereomap');
    
    if strcmpi(confidence_intervals,'intervals') && strcmpi(R_plot,'on')
      %---- D.2: R values plot
      Hr = newrfigure(Hr);
      hold on;
      d1 = scatter(X0,Y0,90,R_BEST,'filled');
      set(d1,'Marker','s');
      axis([(min(X0)-0.5) (max(X0)+0.5) (min(Y0)-0.5) (max(Y0)+0.5)])
      colormap jet;
      xlabel('X grid number');  ylabel('Y grid number');
      h = colorbar;   ylabel(h,'R values');
      if isempty(GRID_ARB)
        set(gca,'XTick',min(X0):max(X0),'YTick',min(Y0):max(Y0));
      end
      hold off;
      grid on;  box on;
      plot_grids(xgridlabel,ygridlabel,zgridlabel);
      for i = 1: length(X0)
        text(X0(i),Y0(i) - 0.02,['[' num2str(R_MIN(i)) ' , ' num2str(R_MAX(i)) ']'],'HorizontalAlignment','center','VerticalAlignment','top');
      end
      %---- END: D.2: R values plot
      % **Save figure
      saveimage('-stereomap-R');
    elseif strcmpi(confidence_intervals,'bootstraps') && strcmpi(R_plot,'on')
      %---- D.2: R values plot
      Hr = newrfigure(Hr);
      for j = 1:size(GRID,1)
        R = R_BOOT{j};
        left = X0(j)/(max(X0) + 1) + 0.01;
        bottom = Y0(j)/(max(Y0) + 1) + 0.05;
        width = 1/length(unique(X0)) - 0.1;
        height = 1/length(unique(Y0)) - 0.1;
        subplot('Position',[left bottom width height]);
        x = 0:0.1:1;
        [n,xout]= hist(R,x);
        h = bar(xout,n);
        set(h,'FaceColor','w','EdgeColor','k')
        best_value = (round(R_BEST(j)*100))/100;
        title(['R best = ' num2str(best_value)])
        grid on; box on;
        xlim([min(R) max(R)]);
        xlabel('R'); ylabel('No. cases');
      end
      %---- END: D.2: R values plot
      % **Save figure
      saveimage('-stereomap-R');
    end
  else
    error(['You selected ''stereomap'' plot type, but the total ' ...
      'dimension of data grid is not 2D (ndim_total = ' num2str(ndim_total) ')']);
  end
  
  
%==== E: Plot data as "STEREOVOLUME" ==============================================================
elseif strcmpi(type, 'stereovolume')
  
  if (ndim_total == 4 && ndim_eff == 3) || (ndim_total == 2 && ndim_eff == 2)  % data must be 3D
    
    % Prepare grid points.
    X0 = Out.GRID(:,COLGRID_EFF(1));  Y0 = Out.GRID(:,COLGRID_EFF(2));
    if ndim_eff == 2
      Z0 = ones(size(X0));
    else
      Z0 = Out.GRID(:,COLGRID_EFF(3));
    end
    
    %---- E.1: Stereovolume plot of stress axis directions
    if strcmpi(confidence_intervals,'intervals')
      Hs = newsfigure(Hs);
      hold on;
      for i=AXESI
        for j=1:length(AZM_MIN)
          [AG,TG] = meshgrid(linspace(AZM_MIN(j,i),AZM_MAX(j,i),10), ...
            linspace(TKO_MIN(j,i),TKO_MAX(j,i),10));
          [XV, YV, ZV] = sph2cart( (pi/2 - AG), TG - pi/2, 0.5);
          surf(X0(j)+XV,Y0(j)+YV,Z0(j)+ZV,'FaceColor',PATCHCOLORS{i},'EdgeColor','k','FaceAlpha',0.4,'EdgeAlpha',0.2);
        end
      end
      
    elseif strcmpi(confidence_intervals,'bootstraps')
      [XX_BEST,YY_BEST] = project(AZM_BEST, TKO_BEST);
      [dummy1,dummy2,Hs] = make_bootstrap_dots(type,GRID,SBOOT,XX_BEST,YY_BEST,X0,Y0,Z0,view_stereonet,view_grid,stext,title_text,symbol,Hs); %#ok<ASGLU>
    end
    
    txt = {'\sigma_1','\sigma_2','\sigma_3'};
    for m=AXESI
      [XV, YV, ZV] = sph2cart( (pi/2 - AZM_BEST), TKO_BEST - pi/2, 0.5);
      line([X0'; X0'+XV(:,m)'],[Y0'; Y0'+YV(:,m)'],[Z0'; Z0'+ZV(:,m)'],'Color',AxesColor{m},'LineWidth',2);
      if strcmpi(stext,'on')
        text(XV(:,m) + X0, YV(:,m) + Y0, ZV(:,m) + Z0,txt{m},'FontSize',20,'Color',AxesColor{m});
      end
    end
    hold off; axis equal;  box on; grid on;
    plot_title(title_text,Out.Caption);
    plot_grids(xgridlabel,ygridlabel,zgridlabel);
    plot_labels(x_label,y_label,z_label);
    set(gca,'XTick',0:max(X0),'YTick',0:max(Y0),'ZTick',0:max(Z0)); % 'Color','none'
    % xlabel('X grid'); ylabel('Y grid'); zlabel('Z grid');
    axis([min(X0)-1 max(X0)+1 min(Y0)-1 max(Y0)-1 min(Z0)-1 max(Z0)+1]); %
    % GK: ^ Pati, for some reason this above does not work very good in some
    % cases...
    axis equal;
    %---- END: E.1: Stereovolume plot of stress axis directions
    % **Save figure
    saveimage('-stereovolume');
    
    if strcmpi(R_plot,'on')  % R plot(s)
      %---- E.2: R Plot for stereovolume
      Hr = newrfigure(Hr);
      hold on;
      d1 = scatter3(X0,Y0,Z0,90,R_BEST,'filled');
      set(d1,'Marker','s');
      axis([(min(X0)-0.5) (max(X0)+0.5) (min(Y0)-0.5) (max(Y0)+0.5) (min(Z0)-0.5) (max(Z0)+0.5)])
      colormap jet;
      xlabel('X grid number');
      ylabel('Y grid number');
      zlabel('Z grid number');
      h = colorbar;
      ylabel(h,'R values');
      set(gca,'XTick',min(X0):max(X0),'YTick',min(Y0):max(Y0),'ZTick',min(Z0):max(Z0));
      hold off;
      grid on; box on;
      plot_grids(xgridlabel,ygridlabel,zgridlabel);
      plot_labels(x_label,y_label,z_label);
      if strcmpi('ConfidenceIntervals','intervals')
        for i = 1: length(X0)
          text(X0(i),Y0(i) - 0.02,Z0(i),['[' num2str(R_MIN(i)) ' , ' num2str(R_MAX(i)) ']'],'HorizontalAlignment','center','VerticalAlignment','top');
        end
      end
      view([-40 30]);
      %---- END: E.2: R Plot for stereovolume
      % **Save figure
      saveimage('-stereovolume-R');
    end
    
  else
    error(['You selected ''stereovolume'' plot type, but the total ' ...
      'dimension of data grid is not 3D (ndim_total = ' num2str(ndim_total) ')']);
  end
  
end % across all type of plots.

%=========================================================================
%=========================================================================
%=========================================================================
function h = newrfigure(h)

global silent_mode;

if ~silent_mode
  h = [h figure];
else
  h = [h figure('Visible','off')];
end

%-------------------------------------------------------------------------
function h = newsfigure(h)

global silent_mode;

if ~silent_mode
  h = [h figure];
else
  h = [h figure('Visible','off')];
end

%-------------------------------------------------------------------------
function saveimage(filename,varargin)

global save_image;
global project_name;

if save_image
  
  % Force smoothing of various objects (may not work correctly!)
%   h = findobj(gca,'Type','line'); set(h,'LineSmoothing','on');
%   h = findobj(gca,'Type','patch'); set(h,'LineSmoothing','on');
%   h = findobj(gca,'Type','surf'); set(h,'LineSmoothing','on');
%   h = findobj(gca,'Type','mesh'); set(h,'LineSmoothing','on');
  
  % Printinf driver properties.
  driver = '-dpng';
  resolution = '-r150';
  extension = '.png';
  
  if nargin == 3
    % Provide additional grid information if necessary.
    GRID = varargin{1};
    i = varargin{2};
    if size(GRID,1) == 1
      % If 0D problem, skip grid information.
      print(driver,resolution,[project_name filename extension]);
    else
      % Provide more or less extensive description of grid points depending on number of dimentions.
      if size(GRID,2) == 2
        gridstr = sprintf('-x%02dy%02d',GRID(i,:));
      elseif size(GRID,2) == 4
        gridstr = sprintf('-x%02dy%02dz%02dt%02d',GRID(i,:));
      else
        error('Unknown number of dimensions');
      end
      print(driver,resolution,[project_name filename gridstr extension]);
    end
  else
    % Ordinary saving with no grid information.
    print(driver,resolution,[project_name filename extension]);
  end
end

%-------------------------------------------------------------------------
function plot_grids(xgridlabel,ygridlabel,zgridlabel)

if ~isempty(xgridlabel)
  set(gca,'XTickLabel',xgridlabel);
end
if ~isempty(ygridlabel)
  set(gca,'YTickLabel',ygridlabel);
end
if ~isempty(zgridlabel)
  set(gca,'ZTickLabel',zgridlabel);
end

%-------------------------------------------------------------------------
function plot_labels(x_label,y_label,z_label)

if ~isempty(x_label)
  xlabel(x_label)
end
if ~isempty(y_label)
  ylabel(y_label);
end
if ~isempty(z_label)
  zlabel(z_label)
end

%-------------------------------------------------------------------------
function makepatch(type,TR_MIN,TR_MAX,PL_MIN,PL_MAX,X0,Y0,PATCHCOLORS,AxesEnabled,confidence_intervals)

% Design and plot patches.
for i = 1:size(TR_MIN,1)
  for j=1:3
    if ~AxesEnabled{j}, continue; end
    if strcmpi(confidence_intervals,'off'), continue; end
    % Procedure for plotting the confidence intervals as patches.
    [P,A] = meshgrid(linspace(PL_MIN(i,j),PL_MAX(i,j),20),...
      linspace(TR_MIN(i,j),TR_MAX(i,j),20));
    TKO = (pi/2 - P*pi/180);
    AZM = A*pi/180;
    [X,Y] = project(AZM,TKO);
    [IDX, C] = kmeans([X(:) Y(:)],2);
    dist = sqrt(sum((C(2,:)-C(1,:)).^2));
    
    if dist > 0.6 % (45°)   % 0.707
      % if distance > 0.707, we have potentially two patches of
      % confidence intervals. Don't know if 0.707 is the appropriate
      % value. In this case, plot confidence intervals
      % as two separate patches.
      if strcmpi(type,'stereonet') % In X0 = Y0 = 0
        plotpatch(X(IDX==1),Y(IDX==1), PATCHCOLORS{j});
        plotpatch(X(IDX==2),Y(IDX==2), PATCHCOLORS{j});
      else
        plotpatch(X(IDX==1) + X0(i),Y(IDX==1) + Y0(i), PATCHCOLORS{j});
        plotpatch(X(IDX==2) + X0(i),Y(IDX==2) + Y0(i), PATCHCOLORS{j});
      end
    else
      % one patch is plotted.
      if strcmpi(type,'stereonet')
        plotpatch(X,Y,PATCHCOLORS{j});
      else
        plotpatch(X + X0(i),Y + Y0(i),PATCHCOLORS{j});
      end
    end
  end
end

%------------------------------------------------------------------------------
function plotpatch(X,Y,COLOR)

X = X(:);
Y = Y(:);
K = convhull(X,Y);
patch(X(K),Y(K),0,'FaceColor',COLOR,'EdgeColor',COLOR*0.9,'FaceAlpha',0.5);

%------------------------------------------------------------------------------
function [XX_BEST,YY_BEST,Hs] = make_bootstrap_dots(type,GRID,SBOOT,XX_BEST,YY_BEST,X0,Y0,Z0,view_stereonet,view_grid,stext,title_text,symbol,Hs)
global Out;
global AxesColor;
global AXESI;
global confidence_intervals;

for i=1:size(GRID,1)
  BOOT = SBOOT{i};
  XX = zeros(length(BOOT),3);  YY = zeros(length(BOOT),3);
  XV = zeros(length(BOOT),3);  YV = zeros(length(BOOT),3); ZV = zeros(length(BOOT),3);
  
  TR1 = BOOT(:,2); PL1 = BOOT(:,3);
  TR2 = BOOT(:,4); PL2 = BOOT(:,5);
  TR3 = BOOT(:,6); PL3 = BOOT(:,7);
  
  ID = false(size(BOOT,1),3); % defined for bootstrap plot optimization
  
  for j = 1:3
    [NULL,M]= unique(BOOT(:,(2*j):(2*j+1)),'rows'); %#ok<ASGLU>
    ID(M,j) = true;
  end
  
  
  if i == 1 || (strcmpi(confidence_intervals,'bootstraps') && strcmpi(type,'stereonet'))
    Hs = newsfigure(Hs);
    if strcmpi(type,'stereomap')
      plot_stereonet(X0, Y0, view_stereonet, view_grid);
    end
  end
  
  if strcmpi(type,'stereonet')
    if strcmpi(type,'stereovolume') == 0
      plot_stereonet(X0(i), Y0(i), view_stereonet, view_grid);
    end
  end
  
  TRs = [TR1 TR2 TR3];  PLs = [PL1 PL2 PL3];
  
  for j = 1:3
    TKO_PP = (pi/2 - PLs(:,j)*pi/180);
    AZM_PP = TRs(:,j)*pi/180;
    if strcmpi(type,'stereovolume')
      [XV(:,j), YV(:,j), ZV(:,j)] = sph2cart( (pi/2 - AZM_PP), TKO_PP - pi/2, 0.5);
    else
      [XX(:,j),YY(:,j)] = project(AZM_PP,TKO_PP);
    end
  end
  
  hold on;
  
  if strcmpi(type,'stereovolume') == 0
    if size(X0,1) == size(GRID,1)
      for m = fliplr(AXESI)
        plot(XX(ID(:,m),m) + X0(i),YY(ID(:,m),m) + Y0(i),'.','Color',AxesColor{m},'MarkerSize',15,'LineWidth',2);
        
      end
      if strcmpi(type,'stereomap'); marker_size = 10; else marker_size = 15; end
      for m = fliplr(AXESI)
        plot(XX_BEST(i,m) + X0(i),YY_BEST(i,m) + Y0(i),'k','MarkerSize',marker_size,'LineWidth',2,'Marker',symbol);
      end
    else
      for m = fliplr(AXESI)
        plot(XX(ID(:,m),m) + X0,YY(ID(:,m),m) + Y0,'.','Color',AxesColor{m},'MarkerSize',15,'LineWidth',2);
      end
      for m = AXESI
        plot(XX_BEST(i,m) + X0,YY_BEST(i,m) + Y0,'k','MarkerSize',18,'LineWidth',2,'Marker',symbol);
      end
    end
    
    txt = {'\sigma_1','\sigma_2','\sigma_3'};
    if strcmpi(stext,'on')
      for m = AXESI
        text(XX(1,m)  + X0(i) - 0.05,YY(1,m)  + Y0(i) + 0.05,txt{m},'FontSize',20,'Color',AxesColor{m});
      end
    end
  else % Stereovolume
    for m = fliplr(AXESI)
      plot3(XV(ID(:,m),m) + X0(i),YV(ID(:,m),m) + Y0(i),ZV(ID(:,m),m) + Z0(i),'Color',AxesColor{m},'LineWidth',2,'Marker','.','LineStyle','none');
    end
    
    txt = {'\sigma_1','\sigma_2','\sigma_3'};
    if strcmpi(stext,'on')
      for m = AXESI
        text(XV(1,m)  + X0(i) - 0.05,YV(1,m)  + Y0(i) + 0.05,ZV(1,m)  + Z0(i) + 0.05,txt{m},'FontSize',20,'Color',AxesColor{m});
      end
    end
  end
  
  plot_title(title_text,Out.Caption);
  hold off; axis equal; grid on
  if strcmpi(type,'stereonet');  axis off; end
  
  % **Save figure
  if strcmpi(type,'stereomap');  continue; end
  saveimage('-stereonet',GRID,i);

end

%------------------------------------------------------------------------------
function [X, Y, R] = project(AZM, TKO)

global scale_factor;
global stereo_projection;

I = TKO > pi/2;
AZM(I) = AZM(I) + pi;
TKO(I) = pi - TKO(I);
if strcmpi(stereo_projection, 'schmidt')
  R = sqrt(2)*sin(TKO/2); % Schmidt (Lambert, equal-area)
else
  R = tan(TKO/2);  % Wulff projection (Stereographic, equal-angle)
end
X = scale_factor*R.*sin(AZM);
Y = scale_factor*R.*cos(AZM);

%------------------------------------------------------------------------------
function plot_stereonet(X0, Y0, view_stereonet, view_grid)

global scale_factor;
global grid_color;
global grid_step_azimuth;
global grid_step_plunge;

if strcmpi(view_stereonet, 'off')
  return;
end

for i=1:length(X0)
  x0 = X0(i);
  y0 = Y0(i);
  
  a = (0:4:360)'*pi/180;
  
  % Plot background as patch object.
  Fvc.VERTICES = [scale_factor*cos(a) + x0 scale_factor*sin(a) + y0];
  Fvc.FACES = 1:size(Fvc.VERTICES,1);
  patch(Fvc,'FaceColor','w','EdgeColor','none','LineWidth',2);
  
  % Plot a circle.
  %line(Fvc.VERTICES(:,1),Fvc.VERTICES(:,2),'Color','k','LineWidth',2)
  
  % Plot grid lines if necessary according to the projection.
  if strcmpi(view_grid,'on')
    for GAZM = (0:grid_step_azimuth:360)*pi/180
      GTKO = [0 90] * pi / 180;
      [GX,GY] = project(GAZM, GTKO);
      line(GX+ x0,GY+ y0,'Color',grid_color);
    end
    for GTKO = (0:grid_step_plunge:90)*pi/180
      GAZM = (0:2:360) * pi / 180;
      [GX,GY] = project(GAZM, GTKO);
      line(GX'+ x0,GY'+ y0,'Color',grid_color);
    end
    %[GAZM, GTKO] = meshgrid( (0:grid_step_azimuth:360)*pi/180, (0:grid_step_plunge:90)*pi/180);
    %[GX,GY] = project(GAZM, GTKO);
    %line(GX'+ x0,GY'+ y0,'Color',grid_color);
    %line(GX+ x0,GY+ y0,'Color',grid_color);
  end
  
  patch(Fvc,'FaceColor','none','EdgeColor','k','LineWidth',2);
  set(gcf,'Color','w');
end
% ------------------------------------------------------------------------
function plot_title(title_text,ctext)

if strcmpi(title_text,'Caption');
  title(ctext,'Interpreter','none');
elseif strcmpi(title_text,'off') == 0;
  title(title_text,'Interpreter','none');
end

%-------------------------------------------------------------------------
function [TR_BEST, TR_BOOT, TR_MIN, TR_MAX, JE, JD, TR_MAX_E, TR_MIN_E] = profile_correct2(TR_BEST, TR_BOOT, TR_MIN, TR_MAX)

% Wrap angles to 0-180 degrees.
if ~isscalar(TR_MIN)
  TR_MAX_E = NaN*ones(size(TR_MAX));
  TR_MIN_E = NaN*ones(size(TR_MIN));
  JE = false(size(TR_MAX));
  JD = true(size(TR_MIN));
else
  JE = [];
  JD = [];
  TR_MAX_E = [];
  TR_MIN_E = [];
end

for i=1:3
  I = TR_BEST(:,i) < 0;
  TR_BEST(I,i) = 180 + TR_BEST(I,i);
  if ~isscalar(TR_MIN)
    TR_MIN(I,i) = 180 + TR_MIN(I,i);
    TR_MAX(I,i) = 180 + TR_MAX(I,i);
    
    % Finds right confidence intervals exceeding 180.
    JE(:,i) = TR_MAX(:,i) > 180;
    TR_MAX_E(JE(:,i),i) = TR_MAX(JE(:,i),i) - 180;
    TR_MAX(JE(:,i),i) = 180;
    
    % Find left confidence interval lower than 0.
    JD(:,i) = TR_MIN(:,i) < 0;
    TR_MIN_E(JD(:,i),i) = -TR_MIN(JD(:,i),i); % for the purpos of errorplot
    TR_MIN(JD(:,i),i) = 0;
    
    % Remove intervals that fills the whole space (0-180 deg)
    J = TR_MIN(:,i) == 0 & TR_MAX(:,i) == 180;
    TR_MIN_E(J,i) = NaN;
    TR_MAX_E(J,i) = NaN;
    JD(J,i) = false;
    JE(J,i) = false;
  else
    for j=1:size(TR_BOOT,1)
      I = TR_BOOT{j}(:,i) < 0;
      TR_BOOT{j}(I,i) = 180 + TR_BOOT{j}(I,i);
    end
  end
end

%pause;
%-------------------------------------------------------------------------
function [TR_BEST,TR_BOOT,PL_BEST,TR_MIN,TR_MAX,PL_MIN,PL_MAX] = profile_correct(TR_BEST,TR_BOOT,PL_BEST,TR_MIN,TR_MAX,PL_MIN,PL_MAX)

% Constrain trend values between 0 and 360: No harm to plunge

if size(TR_MIN) ~= [1,1]
  TR_MIN(TR_BEST < 0) = TR_MIN(TR_BEST < 0) + 180;
  TR_MAX(TR_BEST < 0) = TR_MAX(TR_BEST < 0) + 180;
  
  TR_MIN(TR_BEST >= 180) = TR_MIN(TR_BEST >= 180) - 180;
  TR_MAX(TR_BEST >= 180) = TR_MAX(TR_BEST >= 180) - 180;
end

TR_BEST(TR_BEST < 0) = TR_BEST(TR_BEST < 0) + 180;
TR_BEST(TR_BEST >= 180) = TR_BEST(TR_BEST >= 180) - 180;

if iscell(TR_BOOT)
  for m = 1:length(TR_BOOT)
    TR_BOOT{m}(TR_BOOT{m} < 0) = TR_BOOT{m}(TR_BOOT{m} < 0) + 180;
    TR_BOOT{m}(TR_BOOT{m} >= 180) = TR_BOOT{m}(TR_BOOT{m} >= 180) - 180;
  end
end

% Problem A: The vertical axis
for b = 1:3
  for a = 2:size(TR_BEST,1)
    dfA = TR_BEST(a,b) - TR_BEST(a-1,b);
    if dfA > 270 && PL_BEST(a,b) > 60
      TR_BEST(a,b) = TR_BEST(a,b) - 360;
      if size(TR_MIN) ~= [1,1]
        TR_MIN(a,b) = TR_MIN(a,b) - 360;
        TR_MAX(a,b) = TR_MAX(a,b) - 360;
      end
    elseif dfA < -270 && PL_BEST(a,b) > 60
      TR_BEST(a,b) = TR_BEST(a,b) + 360;
      if size(TR_MIN) ~= [1,1]
        TR_MIN(a,b) = TR_MIN(a,b) + 360;
        TR_MAX(a,b) = TR_MAX(a,b) + 360;
      end
    end
    
    if iscell(TR_BOOT)
      % % In case the first bootstrap is opposite to best solution
      if sign(TR_BEST(a,b)) == 1 && sign(TR_BOOT{a}(1,b)) == -1
        TR_BOOT{a}(1,b) = TR_BOOT{a}(1,b) + 360;
      elseif sign(TR_BEST(a,b)) == -1 && sign(TR_BOOT{a}(1,b)) == 1
        TR_BOOT{a}(1,b) = TR_BOOT{a}(1,b) - 360;
      end
      
      for c = 2:length(TR_BOOT{a}) % Type A: bootstrap
        df_boot = TR_BOOT{a}(c,b) - TR_BOOT{a}(c-1,b);
        if df_boot >= 270 && PL_BEST(a,b) > 60
          TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) - 360;
        elseif df_boot <= -270 && PL_BEST(a,b) > 60
          TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) + 360;
        end
      end
    end
    
  end
  
end

% Problem C: The axes on shallow plunges
for b=1:3
  if iscell(TR_BOOT)
    a = 1; % to correct bootstrap from 1st grid
    for c = 2:length(TR_BOOT{a}) % Type c: bootstrap
      df_boot = TR_BOOT{a}(c,b) - TR_BOOT{a}(c-1,b);
      if df_boot >= 150 && PL_BEST(a,b) < 20
        TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) - 180;
      elseif df_boot <= -150 && PL_BEST(a,b) < 20
        TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) + 180;
      end
    end
  end
  for a=2:size(TR_BEST,1)
    df = TR_BEST(a,b) - TR_BEST(a-1,b);
    
    if df >= 174 && PL_BEST(a,b) < 20   % Problem type C: best
      TR_BEST(a,b) = TR_BEST(a,b) - 180;
      if size(TR_MIN) ~= [1,1]
        TR_MIN(a,b) = TR_MIN(a,b) - 180;
        TR_MAX(a,b) = TR_MAX(a,b) - 180;
      end
    elseif df <= -174 && PL_BEST(a,b) < 20 % Problem type C
      TR_BEST(a,b) = TR_BEST(a,b) + 180;
      if size(TR_MIN) ~= [1,1]
        TR_MIN(a,b) = TR_MIN(a,b) + 180;
        TR_MAX(a,b) = TR_MAX(a,b) + 180;
      end
    end
    
    if iscell(TR_BOOT)
      % % In case the first bootstrap is opposite to best solution
      if sign(TR_BEST(a,b)) == 1 && sign(TR_BOOT{a}(1,b)) == -1
        TR_BOOT{a}(1,b) = TR_BOOT{a}(1,b) + 180;
      elseif sign(TR_BEST(a,b)) == -1 && sign(TR_BOOT{a}(1,b)) == 1
        TR_BOOT{a}(1,b) = TR_BOOT{a}(1,b) - 180;
      end
      
      
      for c = 2:length(TR_BOOT{a}) % Type c: bootstrap
        df_boot = TR_BOOT{a}(c,b) - TR_BOOT{a}(c-1,b);
        if df_boot >= 150 && PL_BEST(a,b) < 20
          TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) - 180;
        elseif df_boot <= -150 && PL_BEST(a,b) < 20
          TR_BOOT{a}(c,b) = TR_BOOT{a}(c,b) + 180;
        end
      end
      
    end
  end
end
% project between 0 and 180
TR_BEST(TR_BEST >= 180) = TR_BEST(TR_BEST >= 180) - 180;
%         TR_BEST(TR_BEST < 0) = abs(TR_BEST(TR_BEST < 0));

if size(TR_MIN) ~= [1,1]
  TR_MIN(TR_BEST >= 180) = TR_MIN(TR_BEST >= 180) - 180;
  TR_MAX(TR_BEST >= 180) = TR_MAX(TR_BEST >= 180) - 180;
  
  %              TR_MIN(TR_BEST < 0) = abs(TR_MIN(TR_BEST < 0));
  %              TR_MAX(TR_BEST < 0) = abs(TR_MAX(TR_BEST < 0));
end

if iscell(TR_BOOT)
  for m = 1:length(TR_BOOT)
    TR_BOOT{m}(TR_BOOT{m} >= 180) = TR_BOOT{m}(TR_BOOT{m} >= 180) - 180;
    TR_BOOT{m}(TR_BOOT{m} < 180) = abs(TR_BOOT{m}(TR_BOOT{m} < 180));
  end
end
% %         % END TRIAL FIX AXES
%-------------------------------------------------------------------------
function alpha = SH(S1,S2,S3,R)

% SH() Matlab function for the calculation of the direction of maximum
%      horizontal stress.
%
%   Copyright (C) 1998-2007 Bjorn Lund
%
% This program is free software; you have the permission to use, copy, and
% distribute, either verbatim or with modifications, either without cost or for a
% fee, as long as the copyright notice above is distributed with the software
% and you clearly state modifications made to the original code.
% If you use this algorithm for work that is to be published, please cite
% the paper:
% Lund and Townend, (2007). Calculating horizontal stress orientations
%    with full or partial knowledge of the tectonic stress tensor,
%    Geophys. J. Int., 170, 1328-1335, doi: 10.1111/j.1365-246X.2007.03468.x.
%
%---------------------------------------------------------------------------
%   SH()
%      Calculate the direction of maximum horizontal stress using only the
%      directions of the principal stress and R = (S1 - S2)/(S1 - S3).
%      This function uses the methodology described in
%      Lund and Townend, (2007). Calculating horizontal stress orientations
%      with full or partial knowledge of the tectonic stress tensor,
%      Geophys. J. Int., 170, 1328-1335, doi: 10.1111/j.1365-246X.2007.03468.x.
%      In particular, SH() uses Equations 11 and 10 from the paper.
%
%   Input:
%      S1, S2, S3 are the principal stress orientations.
%                 The variables holds the coordinates in the North, East
%                 and Down geographical coordinate system as vectors
%                 of three floats, e.g. S1 = [s1N s1E s1D]
%      R is a single number describing the relative magnitude of S2 with
%        respect to S1 and S3. R = (S1 - S2)/(S1 - S3)
%
%   Returns:
%      The direction of SH from North, angle in radians
%

EPS = 1.0e-8;
UNDEF = -999;

% Calculate the direction of maximum horizontal stress using Eq. 11 in
% Lund and Townend (2007)
Y = 2.0*(S1(1)*S1(2) + (1.0 - R)*(S2(1)*S2(2)));
X = S1(1)*S1(1) - S1(2)*S1(2) + (1.0 - R)*(S2(1)*S2(1) - S2(2)*S2(2));

% Is the denominator (X here) from Eq. 11 in Lund and Townend (2007) zero?
if abs(X) < EPS
  % If so, the first term in Eq. 10 is zero and we are left with the
  % second term, which must equal zero for a stationary point.
  % The second term is zero either if
  % s1Ns1E + (1-R)*s2Ns2E = 0   (A)
  % or if
  % cos(2*alpha) = 0            (B)
  % If (A) holds, the direction of SH is undefined since Eq. 10 is zero
  % irrespective of the value of alpha. We therefore check for (A) first.
  % If (A) holds, R = 1 + s1Ns1E/s2Ns2E unless s2Ns2E = 0, in which case
  % s1Ns1E also has to be zero for (A) to hold.
  %
  if abs(S2(1)*S2(2)) < EPS
    % s2Ns2E = 0
    if abs(S1(1)*S1(2)) < EPS
      alpha = UNDEF;
      return
    else
      alpha = pi/4.0;
    end
  else
    if abs(R - (1.0 + S1(1)*S1(2)/S2(1)*S2(2))) < EPS
      alpha = UNDEF;
      return
    else
      alpha = pi/4.0;
    end
  end
  
  % The denominator is non-zero
else
  alpha = atan(Y/X)/2.0;
end

% Have we found a minimum or maximum? Use 2nd derivative to find out.
% A negative 2nd derivative indicates a maximum, which is what we want.
dev2 = -2.0*X*cos(2.0*alpha) - 2.0*Y*sin(2.0*alpha);
if dev2 > 0
  % We found a minimum. Add 90 degrees to get the maximum.
  alpha = alpha + pi/2.0;
end

% The resulting direction of SH is given as [0,180[ degrees.
if alpha < 0
  alpha = alpha + pi;
end
%if alpha > pi or abs(alpha - pi) < EPS
if alpha > pi || abs(alpha - pi) < EPS
  alpha = alpha - pi;
end
