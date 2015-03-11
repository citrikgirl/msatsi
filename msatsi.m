function [OUT] = msatsi(projectname, TABLE, varargin)
%MSATSI Stress tensor inversion and uncertainty assesment.
%   Calculate stress tensor orientation using focal mechanisms. For
%   full documentation of msatsi.m routine, see http://induced.pl/msatsi 
%   or documentation files provided in MSATSI package.
%

%   Copyright 2013-2014 Patricia Martínez-Garzón <patricia@gfz-potsdam.de>
%                       Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%   $Revision: 1.0.8 $  $Date: 2015.02.04 $ 
% 
% If you use MSATSI in your research, please refer to the following papers:
%
% Martínez-Garzón et al. (2014). Seismol. Res. Lett., 85, 4, doi: 10.1785/0220130189
% Hardebeck and Michael (2006). J. Geophys. Res. Solid Earth 111, B11310, doi 10.1029/2005JB004144.
% Lund and Townend,(2007). Geophys. J. Int., 170, 1328-1335, doi: 10.1111/j.1365-246X.2007.03468.x.
%
%     1.0.8 SilentMode and SaveImage options added. Small corrections to existing code.
%           Correction to the best solution of satsi2d.
%

%
% Interpretation of input parameters.
p = inputParser;
p.addRequired('projectname', @(x) ischar(x) && length(x) < 20);
p.addRequired('TABLE', @(x) isnumeric(x) && (size(x,2) == 5 || size(x,2) == 7)); 
p.addParamValue('Damping', 'on', @(x)any(strcmpi(x,{'on','off'})));  
p.addParamValue('DampingCoeff', 0, @(x) isscalar(x));  
p.addParamValue('ConfidenceLevel', 95, @(x) isscalar(x) && x > 0 && x < 100);  
p.addParamValue('FractionValidFaultPlanes', 0.5, @(x) isscalar(x) && x > 0 && x < 1);  
p.addParamValue('MinEventsNode',20, @(x) isscalar(x) && x > 0);  
p.addParamValue('BootstrapResamplings', 2000, @(x) isscalar(x) && x > 0);  
p.addParamValue('Caption', '', @(x) ischar(x));  
p.addParamValue('TimeSpaceDampingRatio', 1, @(x) isnumeric(x));
p.addParamValue('PTPlots', 'on', @(x)any(strcmpi(x,{'on','off'})));

% Parse input parameters.
p.parse(projectname,TABLE,varargin{:});

damping_coeff = p.Results.DampingCoeff;
% Interpret parsing.
if strcmp(p.Results.Damping,'on')
  damping = true;
else
  damping = false;
  pause(5);
end
PTplots = p.Results.PTPlots;

if strcmp(PTplots,'on')
  PT = true;
else
  PT = false;
end

confidence_level = p.Results.ConfidenceLevel;
fraction_corr_picker = 1 - p.Results.FractionValidFaultPlanes;
min_events_per_node = p.Results.MinEventsNode;
n_bootstrap_resamplings = p.Results.BootstrapResamplings;
ts_damp_ratio = p.Results.TimeSpaceDampingRatio;

% Define if inversion is 2D or 4D
if size(TABLE,2) == 5
  is_2D = true;
  n = 0;
else
  is_2D = false;
  n = 2;
  Z = TABLE(:,3);
  T = TABLE(:,4);
end

% Get the absolute path of msatsi.m and keep the path information only.
msatsi_path = mfilename('fullpath');
msatsi_path = fileparts(msatsi_path);

% Determine platform (Windows / Linux is supported)
archstr = computer('arch');
if strcmp(archstr,'win32') || strcmp(archstr,'win64')
  %  win = true;
  if is_2D
    exe_satsi = [msatsi_path '\satsi2d.exe'];
    exe_tradeoff = [msatsi_path '\satsi2d_tradeoff.exe'];
    exe_bootmech = [msatsi_path '\bootmech2d.exe'];
    exe_bootuncert = [msatsi_path '\bootuncert.exe'];
  else
    exe_satsi = [msatsi_path '\satsi4d.exe'];
    exe_tradeoff = [msatsi_path '\satsi4d_tradeoff.exe'];
    exe_bootmech = [msatsi_path '\bootmech4d.exe'];
    exe_bootuncert = [msatsi_path '\bootuncert.exe'];
  end
elseif strcmp(archstr,'glnx86') || strcmp(archstr,'glnxa64') || strcmp(archstr,'maci64')
  %  win = false;
  if is_2D
    exe_satsi = [msatsi_path '/satsi_2D'];
    exe_tradeoff = [msatsi_path '/satsi_2D_tradeoff'];
    exe_bootmech = [msatsi_path '/bootmech_2D'];
    exe_bootuncert = [msatsi_path '/boot_uncert'];
  else
    exe_satsi = [msatsi_path '/satsi_4D'];
    exe_tradeoff = [msatsi_path '/satsi_4D_tradeoff'];
    exe_bootmech = [msatsi_path '/bootmech_4D'];
    exe_bootuncert = [msatsi_path '/boot_uncert'];
  end
else
  error('Platform is not supported.');
end

% Round numbers.
TABLE = round(TABLE);

% Extract desired information from TABLE.
X = TABLE(:,1);
Y = TABLE(:,2);

% Detect cases of only 1 grid (0D)
single = false;
switch is_2D
  case true
    GRIDS = unique([X Y], 'rows');
    if size(GRIDS,1) == 1
      single = true;
      if GRIDS(1) ~= 0 || GRIDS(2) ~= 0;
        error('For 0D inversion please fill X and Y with 0,0');
      end
      new_len = size(TABLE,1);
      TABLE(end+1:new_len*2,:) = [GRIDS(1,1)*ones(new_len,1) (GRIDS(1,2)+1)*ones(new_len,1) TABLE(:,3:5)];
      X = TABLE(:,1);
      Y = TABLE(:,2);
      if damping == true || (damping == false && damping_coeff ~= 0)
        warning('MSATSI:badParameter','Damping cannot be used with single grid point: Setting damping_coeff = 0');
        damping = false;
        damping_coeff = 0;
      end
    end
  case false
    GRIDS = unique([X Y, Z, T], 'rows');
    if size(GRIDS,1) == 1
      error('4D case does not work with only one grid: Use 2D version');
    end
end

DIP_DIRECTION = TABLE(:,n + 3);
DIP_ANGLE = TABLE(:,n + 4);
RAKE = TABLE(:,n + 5);
comment = 'default';
    
if ~exist(projectname,'dir')
  if ~exist(projectname,'file')
    mkdir(projectname);
  else
    error(['''' projectname ''' folder cannot be created since there is a file with the same name.']);
  end
else
  reply = input(['Folder ' projectname ' already exist. Delete and continue [y/n]?'], 's');
  if strcmp(reply,'y')
    [status, message] = rmdir(projectname,'s');
    if status ~= 1
      error(['There is a problem to delete folder: ' projectname ' [' message ']. Please try another project name or close MATLAB and delete the folder manually.']);
    end
    mkdir(projectname);
  else
    return;
  end
end
% Save initial .sat file
sat_input_file = [projectname '.sat']; 

switch is_2D
  case true
    ap = [];
  case false
    ap = [Z T];
end

savesat(sat_input_file, 'w', comment,[X Y ap DIP_DIRECTION DIP_ANGLE RAKE],is_2D,single);

caption = [p.Results.Caption ' (' projectname ') '];

%==== Plot pictures with P/T axes. ============================================
if PT
    plotaxes(TABLE, projectname, caption,is_2D,single)
end


%======= Calculate damping parameter if necessary. ============================
if damping
    [damping_coeff] = tradeoff(projectname,caption,is_2D,ts_damp_ratio,exe_tradeoff);
end

GRIDS = unique([X Y ap], 'rows');
N_EVENTS = zeros(size(GRIDS,1),1);

for i=1:size(GRIDS,1)
  switch is_2D
    case true
      N_EVENTS(i) = sum(X == GRIDS(i,1) & Y == GRIDS(i,2));
    case false
      N_EVENTS(i) = sum(X == GRIDS(i,1) & Y == GRIDS(i,2) & Z == GRIDS(i,3) & T == GRIDS(i,4));
  end
end

%==== Perform SATSI inversion. ===========================================
sat_output_file = [projectname '/' projectname '.out'];
satsi_cmstr = [sat_input_file ' ' sat_output_file ' ' num2str(damping_coeff)];

disp(['Executing ' upper(exe_satsi)]);
switch is_2D
  case true
    command = [exe_satsi ' ' satsi_cmstr];
  case false
    command = [exe_satsi ' ' satsi_cmstr ' ' num2str(ts_damp_ratio)];
end
disp(command); [status] = system(command);
disp(['Exit status = ' num2str(status) ]);
%==== Read SATSI out file and keep best solutions ========================

 BEST_TENSOR = read_out(projectname,GRIDS,is_2D);
 BEST_TRPL = get_trpl(BEST_TENSOR,is_2D);

%==== Run BOOTMECH (bootstrap resampling) ================================
boot_cmstr = [sat_input_file ' ' num2str(n_bootstrap_resamplings) ' ' num2str(fraction_corr_picker) ' ' num2str(damping_coeff)];
disp(['Executing ' upper(exe_bootmech)]);
switch is_2D
  case true
    callline = [exe_bootmech ' ' boot_cmstr];
  case false
    callline = [exe_bootmech ' ' boot_cmstr ' ' num2str(ts_damp_ratio)];
end

disp(callline); [status] = system(callline);
disp(['Exit status = ' num2str(status)]);

%==== Prepare and Call BOOTUNCERT.EXE ==================================== 
bootstrap_file_temp = [projectname '.sat.slboot'];
boot_uncertainty = [projectname '.summary'];
grid_uncertainty = [projectname '.grid'];
boot_uncertainty_ext = [projectname '.summary_ext'];

if exist(boot_uncertainty,'file'); delete(boot_uncertainty); end  
if exist(grid_uncertainty,'file'); delete(grid_uncertainty); end

fid = fopen(bootstrap_file_temp,'r');
try
tline = fgetl(fid);
n_lines_slboot = n_bootstrap_resamplings * size(GRIDS,1); 
if is_2D
  n_dim_add = 0;
  str_a = '%d %d %f %f %f %f %f %f';
  str_b = '%d %d %f %f %f %f %f %f %f';
else
  n_dim_add = 2;
  str_a = '%d %d %d %d %f %f %f %f %f %f';
  str_b = '%d %d %d %d %f %f %f %f %f %f %f';
end

SLBOOT_TENSOR = NaN*ones(n_lines_slboot,8+n_dim_add); 
SLBOOT_TRPL = NaN*ones(n_lines_slboot,9+n_dim_add);
j = 1;
k = 1;
disp('Reading global .slboot file');
while ischar(tline)
  if mod(j,2) == 1
    A = sscanf(tline,str_a);
  else
    B = sscanf(tline,str_b);
    if single
      if A(2) == 1   % remove the second grid
      else
        SLBOOT_TENSOR(k,:) =  A;
        SLBOOT_TRPL(k,:) = B;
        k = k + 1;
      end
    else
      SLBOOT_TENSOR(k,:) =  A;
      SLBOOT_TRPL(k,:) = B;
      k = k + 1;
    end  
  end
  j = j + 1;
  tline = fgetl(fid);
end
fclose(fid);
catch Me
  fclose('all');
  disp(Me.message);
  error(['Problem with interpreting .slboot file: ' bootstrap_file_temp]);
end

disp(['Executing ' upper(exe_bootuncert)]);
fid3 = fopen(grid_uncertainty,'a'); 
GRID = NaN*ones(size(GRIDS,1),3+n_dim_add);
GRID_REJ = [];
for i=1:size(GRIDS,1)
  
  % Skip processing if not enough events.
  if N_EVENTS(i) < min_events_per_node
    GRID_REJ = [GRID_REJ; GRIDS(i,:)];  %#ok<AGROW>
    continue;
  end
  
  % Store grid points.
  x = GRIDS(i,1);
  y = GRIDS(i,2);
  if(is_2D)
    fprintf(fid3,'%d %d %d\r\n',x,y,N_EVENTS(i));
    GRID(i,:) = [x y N_EVENTS(i)]; 
    sbootfile = sprintf('%d_%d.slboot',x,y);
    I_SEL = SLBOOT_TENSOR(:,1) == x & SLBOOT_TENSOR(:,2) == y;
    I_BEST = BEST_TENSOR(:,1) == x & BEST_TENSOR(:,2) == y;
  else
    z = GRIDS(i,3);
    t = GRIDS(i,4);
    fprintf(fid3,'%d %d %d %d %d\r\n',x,y,z,t,N_EVENTS(i));
    GRID(i,:) = [x y z t N_EVENTS(i)];
    sbootfile = sprintf('%d_%d_%d_%d.slboot',x,y,z,t);
    I_SEL = SLBOOT_TENSOR(:,1) == x & SLBOOT_TENSOR(:,2) == y & SLBOOT_TENSOR(:,3) == z & SLBOOT_TENSOR(:,4) == t;
    I_BEST = BEST_TENSOR(:,1) == x & BEST_TENSOR(:,2) == y & BEST_TENSOR(:,3) == z & BEST_TENSOR(:,4) == t; 
  end
  
  % Generate .slboot file with single grid set.
  if exist(sbootfile,'file')
    delete(sbootfile);
  end
  fid2 = fopen(sbootfile,'w');
  % Print best solutions as first lines
  fprintf(fid2, '%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f\r\n%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f\r\n', ...
    [BEST_TENSOR(I_BEST,(3+n_dim_add):end) BEST_TRPL(I_BEST,(3+n_dim_add):end)]');
  % Print the bootstrap resamplings
  fprintf(fid2, '%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f\r\n%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f\r\n', ...
    [SLBOOT_TENSOR(I_SEL,(3+n_dim_add):end) SLBOOT_TRPL(I_SEL,(3+n_dim_add):end)]');
  fclose(fid2);  
  
  switch is_2D
    case true
      callline = [exe_bootuncert ' ' sbootfile ' ' boot_uncertainty ' ' num2str(confidence_level) ' ' boot_uncertainty_ext ' ' num2str(x) ' ' num2str(y) ' ' num2str(0) ' ' num2str(0)];
    case false
      callline = [exe_bootuncert ' ' sbootfile ' ' boot_uncertainty ' ' num2str(confidence_level) ' ' boot_uncertainty_ext ' ' num2str(x) ' ' num2str(y) ' ' num2str(z) ' ' num2str(t)];
  end
  disp(callline);
  [status] = system(callline);
  disp(['Exit status = ' num2str(status) ]);
  delete(sbootfile) 
end
fclose(fid3);

% Remove rejected grids from SLBOOT_TENSOR and SLBOOT_TRPL.
if ~isempty(GRID_REJ)
  for i=1:size(GRID_REJ,1)
    if is_2D
      x = GRID_REJ(i,1);
      y = GRID_REJ(i,2);
      I_SEL1 = SLBOOT_TENSOR(:,1) == x & SLBOOT_TENSOR(:,2) == y;
      I_SEL2 = SLBOOT_TRPL(:,1) == x & SLBOOT_TRPL(:,2) == y;
    else
      x = GRID_REJ(i,1);
      y = GRID_REJ(i,2);
      z = GRID_REJ(i,3);
      t = GRID_REJ(i,4);
      I_SEL1 = SLBOOT_TENSOR(:,1) == x & SLBOOT_TENSOR(:,2) == y & SLBOOT_TENSOR(:,3) == z & SLBOOT_TENSOR(:,4) == t;
      I_SEL2 = SLBOOT_TRPL(:,1) == x & SLBOOT_TRPL(:,2) == y & SLBOOT_TRPL(:,3) == z & SLBOOT_TRPL(:,4) == t;
    end
    SLBOOT_TENSOR = SLBOOT_TENSOR(~I_SEL1,:);
    SLBOOT_TRPL = SLBOOT_TRPL(~I_SEL2,:);
  end
end

% Modify SLBOOT_TENSOR and SLBOOT_PLTL in case of 0D inversion.
if single
  I_SEL = SLBOOT_TENSOR(:,1) == 0 & SLBOOT_TENSOR(:,2) == 0;
  SLBOOT_TENSOR = SLBOOT_TENSOR(I_SEL,:);
  I_SEL = SLBOOT_TRPL(:,1) == 0 & SLBOOT_TRPL(:,2) == 0;
  SLBOOT_TRPL = SLBOOT_TRPL(I_SEL,:);
end
I_SEL = isnan(SLBOOT_TENSOR(:,1));
SLBOOT_TENSOR = SLBOOT_TENSOR(~I_SEL,:);
I_SEL = isnan(SLBOOT_TRPL(:,1));
SLBOOT_TRPL = SLBOOT_TRPL(~I_SEL,:);

% Save .slboot matrixes in corresponding output files:
fid7 = fopen([projectname '.slboot_tensor'],'w');
fid8 = fopen([projectname '.slboot_trpl'],'w');
switch is_2D
  case true
    fprintf(fid7,'%s %s %s %s %s %s %s %s\n','X','Y','See','Sen','Seu','Snn','Snu','Suu');
    fprintf(fid7,'%d %d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n',SLBOOT_TENSOR');
    fprintf(fid8,'%s %s %s %s %s %s %s %s %s\n','X','Y','Phi','Tr1','Pl1','Tr2','Pl2','Tr3','Pl3');
    fprintf(fid8,'%d %d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n',SLBOOT_TRPL');
  case false
    fprintf(fid7,'%s %s %s %s %s %s %s %s %s %s\n','X','Y','Z','T','See','Sen','Seu','Snn','Snu','Suu');
    fprintf(fid7,'%d %d %d %d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n',SLBOOT_TENSOR');
    fprintf(fid8,'%s %s %s %s %s %s %s %s %s %s %s\n','X','Y','Z','T','Phi','Tr1','Pl1','Tr2','Pl2','Tr3','Pl3');
    fprintf(fid8,'%d %d %d %d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n',SLBOOT_TRPL');
end
fclose(fid7);
fclose(fid8);


%==== Creation of output files =========================
switch is_2D
    case true
        I = GRID(:,3) >= min_events_per_node;
    case false
        I = GRID(:,5) >= min_events_per_node;
end
GRID = GRID(I,:);

% Take the result from .summary file
fid4 = fopen(boot_uncertainty,'r');

SUMMARY_TAB = textscan(fid4,'%s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f'); 

fid5 = fopen([projectname '\' projectname '.summary'],'w');
fprintf(fid5,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', ...
    'PhiBest','PhiMin','PhiMax','Tr1Best','Tr1Min','Tr1Max','Pl1Best','Pl1Min','Pl1Max', ...
    'Tr2Best','Tr2Min','Tr2Max','Pl2Best','Pl2Min','Pl2Max','Tr3Best','Tr3Min','Tr3Max','Pl3Best','Pl3Min','Pl3Max');

SUMMARY = [SUMMARY_TAB{1,2:4} SUMMARY_TAB{1,6:8} SUMMARY_TAB{1,10:12} SUMMARY_TAB{1,14:16} SUMMARY_TAB{1,18:20} SUMMARY_TAB{1,22:24} SUMMARY_TAB{1,26:28}];
if single
    fprintf(fid5,'%1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f\n',SUMMARY(1,:)');
else
    fprintf(fid5,'%1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f\n',SUMMARY');
end
fclose(fid4);
fclose(fid5);
  
% ==================== CLEAN UP AND COPY FILES. ==========================
delete(sat_output_file);
delete(grid_uncertainty);
copyfile(bootstrap_file_temp,[projectname '/']);
delete([projectname '.summary']);
movefile([projectname '.summary_ext'],[projectname '/']);
delete([projectname '.sat.slboot'],[projectname '/' projectname '.sat.slboot']);
movefile([projectname '.slboot_tensor'],[projectname '/']);
movefile([projectname '.slboot_trpl'],[projectname '/']);
if exist(grid_uncertainty,'file'); delete(grid_uncertainty); end
if exist(bootstrap_file_temp,'file'); delete(bootstrap_file_temp); end
for i=1:size(GRIDS,1)
    sbootfile = sprintf('%d_%d.slboot',x,y); 
    if exist(sbootfile,'file')     
        delete(sbootfile);
    end
end

%==== Creation of OUTPUT structure ========================================
OUT = struct;
OUT.Damping = damping;
OUT.DampingCoeff = damping_coeff;
OUT.ConfidenceLevel = confidence_level;
OUT.FractionValidFaultPlanes = fraction_corr_picker;
OUT.MinEventsNode = min_events_per_node;
OUT.BootstrapResamplings = n_bootstrap_resamplings;
OUT.Caption = caption;
OUT.TimeSpaceDampingRatio = ts_damp_ratio;
OUT.PTPlot = PT;
OUT.SLBOOT_TENSOR = SLBOOT_TENSOR;
OUT.SLBOOT_TRPL = SLBOOT_TRPL;
OUT.BOOTST_EXT = load([projectname '/' projectname '.summary_ext']);
if single
    dim = size(TABLE,1)/2;  
    OUT.INPUT_TABLE = TABLE(1:dim,:);
    OUT.SUMMARY_TABLE = SUMMARY(1,:);
    OUT.GRID = GRID(1,:);
    OUT.BEST_TENSOR = BEST_TENSOR(1,:);
    OUT.BEST_TRPL = BEST_TRPL(1,:);
else
    OUT.INPUT_TABLE = TABLE;
    OUT.SUMMARY_TABLE = SUMMARY;
    OUT.GRID = GRID;
    OUT.BEST_TENSOR = BEST_TENSOR;
    OUT.BEST_TRPL = BEST_TRPL;
end
% Save output structure in the directory.
save([projectname '/' projectname '_OUT.mat'],'OUT');
close all;

%=========================================================================
%================== Auxiliary functions ==================================
%=========================================================================

%=========================================================================
% SAVESAT : Save the .sat file
%=========================================================================
function savesat(filename, mode, comment, TABLE,is_2D,single,varargin)
    fid = fopen(filename, mode);
    if nargin == 7 && strcmp(varargin{1},'nohead')
    else
        fprintf(fid,'%s\n', comment);
    end
    switch is_2D
        case true
        switch single
            case false
                fprintf(fid,'%d %d %d %d %d\n',TABLE');
                folder = filename(1:end-4);
                fid2 = fopen([folder '\' filename],'w');
                fprintf(fid2,'%d %d %d %d %d\n',TABLE');
                fclose(fid2);
            case true
                fprintf(fid,'%d %d %d %d %d\n',TABLE');
                dim = size(TABLE,1)/2;  
                folder = filename(1:end-4);
                fid2 = fopen([folder '\' filename],'w');
                fprintf(fid2,'%d %d %d %d %d\n',TABLE(1:dim,:)');
                fclose(fid2);
        end
        case false
        switch single
            case false
                fprintf(fid,'%d %d %d %d %d %d %d\n',TABLE');
                folder = filename(1:end-4);
                fid2 = fopen([folder '\' filename],'w');
                fprintf(fid2,'%d %d %d %d %d\n',TABLE');
                fclose(fid2);
            case true
                error('Single inversion not permitted using 4D configuration');
        end
    end
    fclose(fid);
%=========================================================================
% DEFINE_SECOND_PLANE : Compute second nodal plane dip when are given strike, 
% dip and rake for the first nodal plane with AKI & RICHARD's convention. 
% Angles are in degrees.
% Rewritten from C code to MATLAB by PM. 
% Author: Genevieve Patau
% Source code: GMT package (http://gmt.soest.hawaii.edu/)
%=========================================================================
function [str2, dip2, rake2] = define_second_plane(str1,dip1,rake1)

[str2] = computed_strike1(str1,dip1,rake1);
[dip2] = computed_dip1(str1,dip1,rake1);
[rake2] = computed_rake1(str1,dip1,rake1);
%=========================================================================
% COMPUTED_DIP1: Compute rake in the second nodal plane when strike,dip and rake are 
% given for the first nodal plane with AKI & RICHARD's convention.
% Angles are in degrees.
% Rewritten from C code to MATLAB by PM. 
% Author: Genevieve Patau
% Source code: GMT package (http://gmt.soest.hawaii.edu/)
%=========================================================================
function [dip2] = computed_dip1(str1,dip1,rake1)
   
str1 = str1 * pi / 180; %#ok<NASGU>
dip1 = dip1 * pi / 180;
rake1 = rake1 * pi / 180;

if rake1 == 0
  am = 1.0;
else
  am = rake1 / abs(rake1);
end

dip2 = (acos(am * sin(rake1) * sin(dip1)));
dip2 = dip2 * 180/pi;
%=========================================================================
% COMPUTED_RAKE1: Compute rake in the second nodal plane when strike,dip and rake are 
% given for the first nodal plane with AKI & RICHARD's convention.
% Angles are in degrees.
% Rewritten from C code to MATLAB by PM. 
% Author: Genevieve Patau
% Source code: GMT package (http://gmt.soest.hawaii.edu/)
%=========================================================================
function [rake2] = computed_rake1(str1,dip1,rake1)

EPSIL = 0.0001; % Tolerance index

[str2] = computed_strike1(str1,dip1,rake1);
[dip2] = computed_dip1(str1,dip1,rake1);

str1 = str1 * pi / 180;
dip1 = dip1 * pi / 180;
rake1 = rake1 * pi / 180;

str2 = str2 * pi / 180;
dip2 = dip2 * pi / 180;

if rake1 == 0
  am = 1.0;
else
  am = rake1 / abs(rake1);
end

sd = sin(dip1); cd = cos(dip1);
cd2 = cos(dip2);
ss = sin(str1 - str2); cs = cos(str1 - str2);

if(abs(dip2 - pi/2) < EPSIL)
  sinrake2 = am * cd;
else
  sinrake2 = -am * sd * cs / cd2;  		% cd2 [cos(DIP2)] must be used not cd [cos(DIP1)] */
end
rake2 = atan2(sinrake2, -am * sd * ss);
rake2 = rake2 * 180/pi;
%=========================================================================
% COMPUTED_STRIKE1: Compute strike in the second nodal plane when strike,dip and rake are 
% given for the first nodal plane with AKI & RICHARD's convention.
% Angles are in degrees.
% Rewritten from C code to MATLAB by PM. 
% Author: Genevieve Patau
% Source code: GMT package (http://gmt.soest.hawaii.edu/)
%=========================================================================
function [str2] = computed_strike1(str1,dip1,rake1)
  
  EPSIL = 0.0001; % Tolerance index
  
  str1 = str1 * pi / 180;
  dip1 = dip1 * pi / 180;
  rake1 = rake1 * pi / 180;
  
  cd1 = cos(dip1);
  
  if rake1 == 0
      am = 1.0;
  else
      am = rake1/abs(rake1);
  end
  
  sr = sin(rake1); cr = cos(rake1);
  ss = sin(str1); cs = cos(str1);
  
  if cd1 < EPSIL && abs(cr) < EPSIL    
      str2 = str1 + pi;   
  else    
    temp = cr * cs;
    temp = temp + (sr * ss * cd1);
    sp2 = -am * temp;
    temp = ss * cr;
    temp = temp - sr *  cs * cd1;
    cp2 = am * temp;
    str2 = atan2(sp2, cp2);
    
    [str2] = zero_twopi(str2);
  end
  
  str2 = str2 * 180 / pi; 
%=========================================================================
% ZERO_TWOPI
%=========================================================================
function [str] = zero_twopi(str)
  if str >= 2 * pi
    str = str - 2 * pi;
  elseif str < 0
    str =str + 2 * pi;
  end
%=========================================================================
% DC2AXES: Convert strike/dip/rake into P/T axes directions.
% Source code: GMT package PSMECA (http://gmt.soest.hawaii.edu/)
% Converted from original C code to MATLAB by PM. 
%=========================================================================        
function [TS,TD,PS,PD] = dc2axes(S1,D1,R1)

[S2, D2, R2] = define_second_plane(S1,D1,R1);

pure_strike_slip = false;

if abs(sin(R1*pi/180)) > 0.0001
  IM = floor(R1./abs(R1));
elseif abs(sin(R2*pi/180)) > 0.0001
  IM = floor(R2./abs(R2));
else
  pure_strike_slip = true;
end

if pure_strike_slip
  if cos(R1*pi/180) < 0.0 
    PS = S1 + 45;
    TS = S1 - 45;
  else
    PS = S1 - 45;
    TS = S1 + 45;
  end
  
  I = PS >= 360;
  PS(I) = PS(I) - 360;
  I = TS >= 360;
  TS(I) = TS(I) - 360;
  
  PD = 0;
  TD = 0;
else 
  cd1 =  cos(D1 * pi/180) * sqrt(2);
  sd1 =  sin(D1 * pi/180) * sqrt(2);
  cd2 =  cos(D2 * pi/180) * sqrt(2);
  sd2 =  sin(D2 * pi/180) * sqrt(2);
  cp1 = -cos(S1 * pi/180) * sd1;
  sp1 =  sin(S1 * pi/180) * sd1;
  cp2 = -cos(S2 * pi/180) * sd2;
  sp2 =  sin(S2 * pi/180) * sd2;

  amz = - (cd1 + cd2);
  amx = - (sp1 + sp2);
  amy = cp1 + cp2;
  
  dx = atan2(sqrt(amx * amx + amy * amy), amz) - pi/2;
  px = atan2(amy, - amx);
  if px < 0.0
    px = px + 2*pi;
  end

  amz = cd1 - cd2;
  amx = sp1 - sp2;
  amy = - cp1 + cp2;
  dy = atan2(sqrt(amx * amx + amy * amy), - abs(amz)) - pi/2;
  py = atan2(amy, - amx);
  if amz > 0.0 
    py = py - pi;
  end
    
  if py < 0.0 
    py = py + 2*pi;
  end
  
  if IM == 1 
    PD = dy;
    PS = py;
    TD = dx;
    TS = px;
  else 
    PD = dx;
    PS = px;
    TD = dy;
    TS = py;
  end
    
  TS = TS * 180/pi;
  TD = TD * 180/pi;
  PS = PS * 180/pi;
  PD = PD * 180/pi;
end
%=========================================================================
% PLOTAXES: Plot the P and T axes of the given focal mechanisms.
%=========================================================================  
function plotaxes(TABLE,projectname,caption,is_2D,single)
switch is_2D
    case true
        switch single
            case false
             GRIDS = unique(TABLE(:,1:2), 'rows');
            case true
             GRIDS = [0,0];
        end
    case false
        GRIDS = unique(TABLE(:,1:4), 'rows');
end

for i = 1:size(GRIDS,1)
    
    x = GRIDS(i,1);
    y = GRIDS(i,2);
    
    switch is_2D
        case true
            I = TABLE(:,1) == x & TABLE(:,2) == y;
            n = 0;
        case false
            z = GRIDS(i,3);
            t = GRIDS(i,4);
            I = TABLE(:,1) == x & TABLE(:,2) == y & TABLE(:,3) == z & TABLE(:,4) == t;
            n = 2;
    end
    DIP_DIR = TABLE(I,n + 3);
    STRIKE = DIP_DIR - 90;
    DIP    = TABLE(I,n + 4);
    RAKE   = TABLE(I,n + 5);
    
    TTREND = zeros(length(STRIKE),1);
    TPLUNGE = zeros(length(STRIKE),1);
    PTREND = zeros(length(STRIKE),1);
    PPLUNGE = zeros(length(STRIKE),1);
    
    for k = 1:length(STRIKE)
        [TTREND(k),TPLUNGE(k),PTREND(k),PPLUNGE(k)] = dc2axes(STRIKE(k),DIP(k),RAKE(k));
    end
    
    TTREND = TTREND * pi / 180;
    TPLUNGE = TPLUNGE * pi / 180;
    PTREND = PTREND * pi / 180;
    PPLUNGE = PPLUNGE * pi / 180;
    
    [XT,YT] = project(TTREND,TPLUNGE);
    [XP,YP] = project(PTREND,PPLUNGE);
    
    % Prepare figure of P/T axes.
    figure('Visible','off');
    A = (0:5:360)'*pi/180;
    line(cos(A),sin(A),'Color','k','LineWidth',2)
    
    % Plot the indicative lines of the azimuths
    for A = (0:15:345)*pi/180;
        line([0,sin(A)],[0,cos(A)],'Color',[0.6 0.6 0.6],'LineWidth',1)
    end
    
    % Plot the indicative lines of the plunges
    A = (0:5:360)'*pi/180;
    for R = 1/6:1/6:(1-1/6)
        line(R*cos(A),R*sin(A),'Color',[0.6 0.6 0.6],'LineWidth',1)
    end
    
    hold on
    h1 = plot(XT,YT,'bo','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
    h2 = plot(XP,YP,'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
    hold off;
    axis equal
    set(gca,'Color','w');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'XColor','w');
    set(gca,'YColor','w');
    set(gcf,'Color','w')
    
    legend([h1,h2],{'Tension','Pressure'},'Location','SouthOutside')
    switch is_2D
        case true
            grid_info = sprintf('(x= %02.0f,y=%02.0f)',x,y);
            content = strcat(caption,grid_info);
            title([caption content]);
            fileout = [projectname '/' projectname grid_info  '.png'];
        case false
            grid_info = sprintf('(x= %02.0f,y=%02.0f,z=%2.0f,t=%2.0f)',x,y,z,t);
             content = strcat(caption,grid_info);
            title([caption content]);
            fileout = [projectname '/' projectname grid_info '.png'];
    end
    print('-r300','-dpng',fileout);
end
%=========================================================================
% PROJECT: Project according to Schmidt projection (lower hemisphere). Points with 
% Takeoff angles > 90 deg are reverted. 
% Attention!! AZM and TKO are only mathematically, not always have physical
% meaning except for case of polarities!!
%=========================================================================  
function [X,Y,R] = project(AZM,TKO)

  TKO = pi/2 - TKO; % Plunge is calculated in same way as stresses
  I = TKO > pi/2;
  AZM(I) = AZM(I) + pi;
  TKO(I) = pi - TKO(I);   %  TKO(I) = pi - TKO(I); 
  R = sqrt(2)*sin(TKO/2);   % schmidt radius    
  X = R.*sin(AZM);
  Y = R.*cos(AZM);
%=========================================================================
% TRADEOFF: Calculate optimal damping parameter.
%=========================================================================  
function [damping, status, result] = tradeoff(projectname, caption, is_2D, ts_damp_ratio,exe_tradeoff)

damp_parameters = [0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6];
% 0.2, 0.3, , 8, 10, 20, 50
OUTPUT = zeros(length(damp_parameters),3);
input_file = [projectname '.sat'];

disp('Calculating damping parameter.');
for i = 1:length(damp_parameters)
    output_file = [projectname '/' projectname '_' sprintf('%1.2f',damp_parameters(i)) '.trd'];
    switch is_2D
        case true
            [status, result] = system([exe_tradeoff ' ' input_file ' ' output_file ' ' num2str(damp_parameters(i))]);
        case false
            [status, result] = system([exe_tradeoff ' ' input_file ' ' output_file ' ' num2str(damp_parameters(i)) ' ' num2str(ts_damp_ratio)]);
    end
    disp(['Exit status = ' num2str(status) ]);
end

for i = 1:length(damp_parameters)
    output_file = [projectname '/' projectname '_' sprintf('%1.2f',damp_parameters(i)) '.trd'];
    TEMP = load(output_file);
    OUTPUT(i,:) = [damp_parameters(i) TEMP];
    delete(output_file);
end

% Selection of the most apropiate damping parameter based on the tradeoff
% curve.
X = OUTPUT(:,2); % Tradeoff curve points.
Y = OUTPUT(:,3);

figure('Visible','off');
hold on;
plot(X,Y,'ok');
text(X,Y,num2cell(OUTPUT(:,1)));

% Normalization
XX = (X-min(X))/(max(X)-min(X));
YY = (Y-min(Y))/(max(Y)-min(Y));
R = sqrt(YY.^2+XX.^2);

% Simple finding of inflection point (we believe it is a hyperbola...)
i = find(R==min(R),1,'first');
plot(X(i),Y(i),'+b','MarkerSize',20);

damping = OUTPUT(i,1);

hold off;
box on;
xlabel('Data misfit');
ylabel('Model length');
title({'Trade-off curve',caption});
print('-dpng','-r300',[projectname '/' projectname '_tradeoff.png']);
close all;
%=========================================================================
% BEST_TENSOR: Reads output file with the best stress tensor solutions
%=========================================================================
function BEST_TENSOR = read_out(projectname,GRIDS,is_2D)

archstr = computer('arch');
if strcmp(archstr,'win32') || strcmp(archstr,'win64')
    fid = fopen([projectname '\' projectname '.out']);
elseif strcmp(archstr,'glnx86') || strcmp(archstr,'glnxa64') || strcmp(archstr,'maci64')
    fid = fopen([projectname '/' projectname '.out']);
else
  error('Platform is not supported.');
end

tline = fgetl(fid);
if is_2D   % make it only str for fscanf
    BEST_TENSOR = zeros(size(GRIDS,1),8);
    k = 1;
    while ischar(tline)      % While we have string lines...
        if length(tline) <= 2 | length(tline) == 30 | length(tline) == 19 | length(tline) == 27
            % Headers of the file
        elseif length(tline) == 67
            % fscanf think about
            xb = str2double(strtrim(tline(1:3)));
            yb = str2double(strtrim(tline(5:7)));
            see = str2double(strtrim(tline(9:17)));
            sen = str2double(strtrim(tline(19:27)));
            seu = str2double(strtrim(tline(29:37)));
            snn = str2double(strtrim(tline(39:47)));
            snu = str2double(strtrim(tline(49:57)));
            suu = str2double(strtrim(tline(59:67)));
            
            BEST_TENSOR(k,:) = [xb yb see sen seu snn snu suu];
            k = k + 1;
            if k > size(GRIDS,1)
                break
            end
        end
        tline = fgetl(fid);
    end
else
   BEST_TENSOR = zeros(size(GRIDS,1),10);
    k = 1;
    while ischar(tline)      % While we have string lines...
        if length(tline) <= 2 | length(tline) == 30 | length(tline) == 19 | length(tline) == 31
            % Headers of the file
        elseif length(tline) == 75
            xb = str2double(strtrim(tline(1:3)));
            yb = str2double(strtrim(tline(5:7)));
            zb = str2double(strtrim(tline(9:11)));
            tb = str2double(strtrim(tline(13:15)));
            see = str2double(strtrim(tline(17:25)));
            sen = str2double(strtrim(tline(27:35)));
            seu = str2double(strtrim(tline(37:45)));
            snn = str2double(strtrim(tline(47:55)));
            snu = str2double(strtrim(tline(57:65)));
            suu = str2double(strtrim(tline(67:75)));
            
            BEST_TENSOR(k,:) = [xb yb zb tb see sen seu snn snu suu];
            k = k + 1;
            if k > size(GRIDS,1)
                break
            end
        end
        tline = fgetl(fid);
    end
end
fclose(fid);
%=========================================================================
% BEST_TRPL: Calculate trend and plunges of the best solution stress tensors
%=========================================================================
function BEST_TRPL = get_trpl(BEST_TENSOR,is_2D)
% Select one stress tensor for each grid
if is_2D
    I =0;
else
    I = 2;
end    
BEST_TRPL = zeros(size(BEST_TENSOR,1),9+I);
for i = 1:size(BEST_TENSOR,1)
    STRESS_LINE = BEST_TENSOR(i,3+I:end);
    STRESS = [STRESS_LINE(1) STRESS_LINE(2) STRESS_LINE(3); ...
        STRESS_LINE(2) STRESS_LINE(4) STRESS_LINE(5); ...
        STRESS_LINE(3) STRESS_LINE(5) STRESS_LINE(6)];
    
    % Calculate eigenvalues and phi
    [VECS,D] = eig(STRESS);
    LAM = diag(D);
    [LAM,IX] = sort(LAM); % Sort from smaller to largest, since tension >0 in SATSI
    VECS = VECS(:,IX);
    phi = (LAM(2) - LAM(3))/(LAM(1) - LAM(3));
    % Define trend and plunge of each eigenvector
    E = VECS(1,:); N = VECS(2,:); U = VECS(3,:);
    Z = sqrt(E.^2 + N.^2);
    PPLG = atan2(-U,Z) * 180/pi;
    for j = 1:3
        if PPLG(j) < 0
            PPLG(j) = -1*PPLG(j);
            E(j) = -1*E(j);
            N(j) = -1*N(j);
        end
    end
    PDIR = atan2(E,N) * 180/pi;
    
    BEST_TRPL(i,:) = [BEST_TENSOR(i,1:2+I)...
        phi PDIR(1) PPLG(1) PDIR(2) PPLG(2) PDIR(3) PPLG(3)];
end






