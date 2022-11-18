function NetFlags = settings_conn_default(settings_fun)

try

NetFlags.lookup_file = '';
NetFlags.subdatafield = 'hasRS';

NetFlags.zscoreTS = true;
NetFlags.detrendTS = [1, 1]; % detrend the time series and detrend the confounds

NetFlags.parpool = 0; % how many parallel processes to use: depends on the server used

NetFlags.doPredictions = false;
NetFlags.Group = ''; % this will be used as directory name to save results

NetFlags.usegpu = false;
NetFlags.project    = 'results'; % a sub directory inside results
NetFlags.Analysis   = {'Rew_NoGMmask_NoPCA_lWMCSF_HP'};
NetFlags.GroupNotIsPat      = [];  % use [] if the groups are already coded by isPat
NetFlags.GroupPrintNames    = {'young','old'};

% applied to whole dataset: just ignore these
NetFlags.Movement2Correct   = ''; %'DVARS';
NetFlags.Volume2Correct     = []; %'TBV'; 
NetFlags.SDthreshold        = [Inf Inf];
NetFlags.MinkowskiExponent  = 2;

% how to subsample the ROI and WB
NetFlags.subsampleWB = 1;
NetFlags.subsampleROI = 1;

NetFlags.conn_metrics = {};

NetFlags.ConfoundsTS = 'WM+CSF+GS'; % time-series confounds taken from Confound_sub.mat
NetFlags.Confounds = {};  % removes from connectivity over all subjects can be 'age' and/or 'gender'
NetFlags.CorrelationMasking    = 'UN' ; % First letter: N one / U ncorrected / C orrected between group differences; Second letter: N one / E ffectsize
NetFlags.EffectSizeThreshold   = -Inf;
NetFlags.Correlation.Type      = 'Pearson'; % Spearman / Pearson

% settings for prediction
NetFlags.pred.Method = 'RVM';
NetFlags.pred.CVarg = [];
NetFlags.pred.CVarg.NumFolds = 10;
NetFlags.pred.CVarg.NumRepeats = 10;
NetFlags.pred.Methodarg = []; % define any hyperparam grid and inner-loop settings here, not needed for RVM
NetFlags.pred.Confound = [];
NetFlags.pred.Preprocarg = [];
NetFlags.pred.Preprocarg.preprocX = 'zscore';
NetFlags.pred.Preprocarg.preprocy = 'zscore';
NetFlags.pred.replaceNaNwithZeroX  = true;

NetFlags.ConInfo.ForRegressionAll  = [2]; % these are 'Age'
NetFlags.ConInfo.ForRegressionPat  = []; % Those in ForRegressionAll are also used
NetFlags.ConInfo.ForRegressionCon  = []; % Those in ForRegressionAll are also used
NetFlags.ConInfo.FilterOutliers    = true; 
NetFlags.ConInfo.SDthreshold       = 8;
NetFlags.ConInfo.ExForNoCov        = false;

% Kaustubh: some additional options for flexibility
NetFlags.data = {''};
NetFlags.global_mask_to_use = 'CAT12_02'; % mask to subset voxels for all subjects ('MNI152GM', 'CanLabGM', 'FSL04', 'None')
NetFlags.global_mask_to_use_vbm = 'CAT12_02_vbm';
NetFlags.regress_confounds  = true; % whether to regress out the confounds (true or false)
NetFlags.filter_freq        = true; % whether to filter the time-series in frequency (true or false)
NetFlags.filterFrequency    = [0.01, 0.1]; % use from Amico and Goni, 2018 Sci Reports
NetFlags.sphere_radius      = 5;    % radius of the sphere
NetFlags.sphere_boundary    = '<='; % how to get the sphere boundary ('<=' or '<')
NetFlags.preproc_data       = ''; % preprocessed data to use ('Smooth_5mm', 'Normalised', 'Norm_sub')

% for prediction analysis
NetFlags.slope_adjustment   = false; % whether to use slope adjustment (true or false)
NetFlags.scale_cov          = true; % whether to scale the dependent variable to mean 0 and std 1 
NetFlags.datareduction      = 'mean'; % either 'ev' for eigenvariate or 'mean' for mean
NetFlags.connectome_cov     = 'COV';

% calculate connectome using ROIs in a nifti file
NetFlags.connectome_nii     = {};
NetFlags.niiref = ''; % reference nifti file

%NetFlags.conndim_dist = [10,50,90];
NetFlags.maxSubConnectome  = inf; % how many first subjects to use
NetFlags.conndim_randWB = 0;
%NetFlags.centrality_coord = 'VOIs/unused/Power2013.txt';
%NetFlags.smoothness = 100;

% nifti file for vbm extraction
NetFlags.vbm_nii = '';
%NetFlags.vbm_nii = '/data/BnB2/TOOLS/masks/Power/Power_5mm_181x217x181.nii';

NetFlags.datadir = '';
NetFlags.confoundsdir = 'data/group/appliedml/data/HCP1200_Confounds'; % empty assumes confounds are in the datadir

usersettings = nan;

if isstruct(settings_fun)
    usersettings = settings_fun;
else
    fprintf('getting user settings: %s\n', settings_fun)
    [~,~,ext] = fileparts(settings_fun);
    switch ext
        case {'','.m'}
            settings_func = str2func(settings_fun);
            usersettings = settings_func();
        case '.mat'
            usersettings = load(settings_fun);
        otherwise
            error('cannot deal with user settings!')
    end
end

required_settings = {'lookup_file', 'datadir', 'spmdir', 'data', 'preproc_data', 'conn_metrics', ...
    'subdatafield', 'connectome_nii', 'savesubfile'};

if isstruct(usersettings)
  fs = fields(usersettings);
  for f=1:length(required_settings)
      assert(ismember(required_settings{f}, fs), ['missing user_setting: ' required_settings{f}])
  end

  for f=1:length(fs)
    NetFlags.(fs{f}) = usersettings.(fs{f});
  end
end

fprintf('data dir: %s\n', NetFlags.datadir)

assert(isfolder(NetFlags.datadir), ['datadir is not a folder: ' NetFlags.datadir])
assert(isfile(NetFlags.lookup_file), ['lookup_file is not a file: ' NetFlags.lookup_file])
assert(isfolder(NetFlags.spmdir), ['spmdir is not a folder: ' NetFlags.spmdir])
if ~isdeployed
    addpath(NetFlags.spmdir)
    addpath(fullfile(NetFlags.spmdir, 'compat'))
end

if ischar(NetFlags.data)
    NetFlags.data = {NetFlags.data};
end

if ~isfield(NetFlags, 'Group') || isempty(NetFlags.Group)
    [~, NetFlags.Group, ~] = fileparts(NetFlags.lookup_file);
end

if ~isfield(NetFlags, 'outputdir') || isempty(NetFlags.outputdir)
    NetFlags.outputdir = fullfile(pwd, 'output', NetFlags.Group);
else
    NetFlags.outputdir = fullfile(NetFlags.outputdir, NetFlags.Group);
end
[~,~,~] = mkdir(NetFlags.outputdir);
fprintf('output dir: %s\n', NetFlags.outputdir)
groupwrite = true;
if isfield(NetFlags, 'groupwrite') && ~isempty(NetFlags.groupwrite)
     groupwrite = NetFlags.groupwrite;
end
if groupwrite
     [status,msg,msgID] = fileattrib(NetFlags.outputdir,'+w', 'g');
     if status==0
         fprintf('WARN: groupwrite attribute could not be set for %s\n', NetFlags.outputdir)
     end
end


catch ME

getReport(ME)
throw(ME)

end
