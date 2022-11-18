% connectivity dimensionality given a nifti file
% saves data in a per subject structure
% data is saved after each subject is processed
% Kaustubh 23 Oct 2017: changed to handle 4D data
function niiprop = kp_conn_glocal_sub_noparfor(nii,niiwb,NetFlags,subsampleWB,subsampleROI,subs,createnii)
% nii  : cell array with paths to nii files
% niiwb: path to niiwb file

%assert(exist('laplace_pca','file')==2)

switch nargin
    case {0,1,2}
        error('give me at least 3 arguments')
    case 3
        subsampleWB=1; subsampleROI=1; subs=[]; createnii=true;
    case 4
        subsampleROI=1; subs=[]; createnii=true;
    case 5
        subs=[]; createnii=true;
    case 6
        createnii=true;
end

cwd=pwd;

groupwrite = true;
if isfield(NetFlags, 'groupwrite') && ~isempty(NetFlags.groupwrite)
     groupwrite = NetFlags.groupwrite;
end

[NetFlags, des] = analysis_settings(NetFlags.Analysis{1}, NetFlags);

if ischar(nii)
    nii = {nii};
end

if isfield(NetFlags, 'subsampleWB')
    subsampleWB = NetFlags.subsampleWB;
end
if isfield(NetFlags, 'subsampleROI')
    subsampleROI = NetFlags.subsampleROI;
end
fprintf('sumsample ROI:%d WB:%d\n', subsampleROI, subsampleWB)

usegpu = false;
if isfield(NetFlags, 'usegpu')
    usegpu = NetFlags.usegpu;
end
if usegpu
  fprintf(' (gpu) ')
end

savesubfile = false;
if isfield(NetFlags, 'savesubfile') && ~isempty(NetFlags.savesubfile)
    savesubfile = NetFlags.savesubfile;
end

if isfield(NetFlags, 'outputdir') && ~isempty(NetFlags.outputdir)
    outputdir = NetFlags.outputdir;
else % for backward compatibility reasons
    outputdir = fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'));
end
if ~isfolder(outputdir)
    [~,~,~] = mkdir(outputdir);
end

conn_metrics = {'pairwise', 'reho', 'reho10', 'rehots', 'fcnesswithin', 'fcness', 'alff', 'graph', 'tsnr', 'sampen', 'hacf'};
% alff also calculates falf% alff also calculates falfff
% hacf = high-amplitude co-fluctuations
if isfield(NetFlags, 'conn_metrics') && ~isempty(NetFlags.conn_metrics)
  conn_metrics = NetFlags.conn_metrics;
end
assert(~isempty(conn_metrics), 'conn_metrics must not be empty')
if ischar(conn_metrics)
    conn_metrics = {conn_metrics};
end
conn_metrics = lower(conn_metrics)

nsub = min(length(NetFlags.sub), NetFlags.maxSubConnectome);

tmpdir = tempdir;
if isfield(NetFlags, 'tempdir') && ~isempty(NetFlags.tempdir)
    tmpdir = NetFlags.tempdir;
end

% assuming that the mask and subject data in same space
if isempty(niiwb)
    assert(isfield(NetFlags, 'global_mask_to_use') && ~isempty(NetFlags.global_mask_to_use), 'provide global_mask_to_use')
    niiwb = NetFlags.global_mask_to_use;
end
[wbgm, wbgmfile, Vgm] = get_mask(niiwb);

Vref = [];
Dref = [];
if ~isempty(Vgm)
    Vref = Vgm;
    Vref = Vref(1);
    Dref = Vref.dim;
elseif isfield(NetFlags, 'niiref') && ~isempty(NetFlags.niiref)
    Vref = spm_vol(NetFlags.niiref);
    Vref = Vref(1);
    Dref = Vref.dim;
else % get reference from the atlas files
    for ii=1:length(nii)
        [~,~,ext] = fileparts(nii{ii});
        switch ext
            case {'.nii', '.gz'}
                Vref = spm_vol(nii{ii});
                Vref = Vref(1);
                Dref = Vref.dim;
                break
            otherwise
                Vref = [];
                Dref = [];
        end
    end
end
assert(~isempty(Vref), 'Vref turned out empty')
assert(~isempty(Dref), 'Dref turned out empty')
assert(length(Dref)==3, 'Dref must be 3D volume')
fprintf('ref nifti dim: %d %d %d\n', Dref(1), Dref(2), Dref(3))

% get ROI GM voxels
gm = zeros(Dref);
gm(1:subsampleROI:end,1:subsampleROI:end,1:subsampleROI:end) = 1;
if ~isempty(wbgm)
    gm = gm .* wbgm;
end
IndROI = find(gm>0);

if isfield(NetFlags, 'rmvox')
    rmvox = NetFlags.rmvox;
    mask = ones(Dref);
    rmD = spm_read_vols(spm_vol(rmvox{1}));
    torm = ismember(rmD(:), rmvox{2});
    fprintf('%d voxels removed!\n', sum(torm))
    mask(torm) = 0;
    gm = gm.*mask;
end

% get WB GM voxels
gmWB = zeros(Dref);
gmWB(1:subsampleWB:end,1:subsampleWB:end,1:subsampleWB:end) = 1;
if ~isempty(wbgm)
    gmWB = gmWB .* wbgm;
end
IndWB = find(gmWB>0)';

hemi = [];
if ~isempty(Vref)
    hemi = hemisphere_mask(Vref);
    IndLeft = find(hemi==1);
    IndRight = find(hemi==3);
    IndWBLeft = intersect(IndWB, IndLeft);
    IndWBRight = intersect(IndWB, IndRight);
end

XYZ = [];
[XYZ(:,1), XYZ(:,2), XYZ(:,3)] = ind2sub(Dref, find(ones(Dref)));
XYZ = XYZ';
XYZmm = Vref.mat * [XYZ; ones(1,size(XYZ,2))];
XYZmm = XYZmm(1:3,:);

% get properties of each nifti
niiprop = [];
niinames = {};
IndAll = IndWB; % for all niftis and WB

confdir='';
datalad=false;
if isfield(NetFlags, 'datalad') && ~isempty(NetFlags.datalad)
    datalad = NetFlags.datalad;
end

if datalad
    fprintf('datalad will be used\n')
end

datasetdir = NetFlags.datadir;
settings = kp_get_settings(NetFlags);

logdir = fullfile(outputdir, 'log');
if ~isfolder(logdir)
    [~,~,~] = mkdir(logdir);
end
if groupwrite
    [status,msg,msgID] = fileattrib(outputdir,'+w', 'g');
    if status==0
        fprintf('WARN: groupwrite attribute could not be set for %s\n', logdir)
    end
end

save(fullfile(logdir, ['settings' datestr(now,'ddmmyyyy_HHMMSS') '.mat']), '-struct', 'NetFlags')

% use the gray matter Volume as reference
% we rely on the user to provide a proper mask file
V = Vref; % Vref is the GM mask volume if provided by user
for i=1:length(nii)
    [~,niiname,niiext] = fileparts(nii{i});
    niiname = matlab.lang.makeValidName(niiname);
    niinames{i} = niiname;
    niiprop.(niiname) = [];
    niiprop.(niiname).name = niiname;
    niiprop.(niiname).file = nii{i};
    fprintf('connectome from nifti file: %s\n', niiname);
    [NetFlags, des] = analysis_settings(NetFlags.Analysis{1}, NetFlags);
    destdir = fullfile(outputdir, niiname);
    niiprop.(niiname).destdir = destdir;
    zfile = fullfile(destdir, ['conn_' settings '.mat']);
    niiprop.(niiname).zfile = zfile;

    switch niiext
        case {'.nii', '.gz'}
            [Indnii, VOInii] = nii2voiind(nii{i}, NetFlags, V, 0);
        case '.txt'
            [Indnii, VOInii] = coord2sphereind(nii{i}, NetFlags, V, 0);
	    otherwise
             error('unknown nii file type: expected .nii or .txt')
    end
    ii = ismember(Indnii, IndROI);
    assert(length(ii)==length(Indnii))
    fprintf('keeping %d GM voxels out of %d total\n', sum(ii), length(ii))
    niiprop.(niiname).Ind = Indnii(ii);
    niiprop.(niiname).VOI = VOInii(ii);
    if isfield(NetFlags, 'conndim_randVOI') && NetFlags.conndim_randVOI==true
        fprintf('nifti %s: randomizing VOI!!!\n', niiname);
        niiprop.(niiname).zfile = fullfile(destdir, [niiname '_conn_randVOI_' des '.mat']);
        niiprop.(niiname).VOI = randsample(niiprop.(niiname).VOI, length(niiprop.(niiname).VOI));
    end

    [~, ~, ~] = mkdir(destdir);
    if groupwrite
        [status,msg,msgID] = fileattrib(destdir,'+w', 'g');
        if status==0
            fprintf('WARN: groupwrite attribute could not be set for %s\n', destdir)
        end
    end

    nROI = numel(unique(niiprop.(niiname).VOI));
    IndAll = [IndAll, niiprop.(niiname).Ind];

    niiprop.(niiname).('VOIhemi') = nan(1, nROI);
    for ana=1:nROI
        xQ = niiprop.(niiname).Ind(niiprop.(niiname).VOI==ana);
        niiprop.(niiname).VOIhemi(ana) = mode(hemi(xQ));
    end
end
IndAll(isnan(IndAll)) = []; % remove NaNs that might come from the niftis
[IndAll,ii] = sort(unique(IndAll)); % unique also sorts but sort anyway in case MATLAB changes its mind
XYZAll = XYZ(:,IndAll);
XYZAllmm = XYZmm(:,IndAll);

session = '';
poolsize = 0;
ppool = gcp('nocreate'); % If no pool, do not create new one.
if ~isempty(ppool)
    poolsize = ppool.NumWorkers;
end

forcecompute = false;
if isfield(NetFlags, 'forcecompute') && ~isempty(NetFlags.forcecompute)
    forcecompute = NetFlags.forcecompute;
end

substouse = true(1, nsub);
conn_nii = containers.Map();
substrings = {};
% keep subject as the outer loop so we read in the data only once
% this will save quite some time, especially for data sets with many time-points like HCP
for s=1:nsub
    subid = NetFlags.sub{s};
    substr = num2str(subid);
    if ~isvarname(substr)
        substr = matlab.lang.makeValidName(substr);
    end
    substrings = [substrings, substr];
    % check if any nifti is processed
    nifti_processed = false(1, length(nii));
    for inii=1:length(nii) % for each nifti with VOIs
        % set properties for this nii
        niiname = niinames{inii};
        if savesubfile
            zfile = fullfile(niiprop.(niiname).destdir, [substr '_conn_' settings '.mat']);
            if isfile(zfile)
                nifti_processed(inii) = true;
            end
        else
            zfile = niiprop.(niiname).zfile;
            conn = struct();
            if ~isKey(conn_nii, niiname)
                conn_nii(niiname) = conn;
                if exist(zfile, 'file')==2
                    fprintf('reading conn file: %s\n', zfile)
                    conn_nii(niiname) = load(zfile);
                end
            end
            % check if we have this subject
            if isfield(conn_nii(niiname), substr)
                nifti_processed(inii) = true;
            end
        end % savesubfile
    end

    if all(nifti_processed)
        fprintf('sub %s %d/%d all nifti processed, continuing to next one\n', substr, s , nsub)
        substouse(s) = false;
        continue
    end

    subfiles = {};
    if isfield(NetFlags, 'subfile') && ~isempty(NetFlags.subfile)
        subfiles = NetFlags.subfile{s};
        if ischar(subfiles)
            subfiles = {subfiles};
        end
    end
    if ~datalad
        %Vols = check_data_sub_4D_concat(NetFlags.sub{s}, NetFlags.SubDir{s}, NetFlags, subfiles);
        %if isempty(Vols)
        %    substouse(s) = false;
        %    continue
        %end
    end
end % for nsub

for s=1:nsub
    if ~substouse(s)
        if ~forcecompute
            fprintf('skipping %s: output files exists, forcecompute=%d (set forcecompute=true to override)\n', ...
            substrings{s}, forcecompute)
            continue
        end
    end
    substr = substrings{s};
    timerSub = tic;

    % get subject specific GM voxels
    myIndAll = IndAll;
    if isfield(NetFlags, 'subject_gm_mask') && NetFlags.subject_gm_mask>0
        subInd = subject_gm_voxels(NetFlags, s, nii{1}, NetFlags.subject_gm_mask);
        myIndAll = intersect(myIndAll, subInd);
        fprintf('(subGM:%d/%d) ', length(myIndAll), length(IndAll))
    end

    Vols = [];
    TSAll = [];
    fprintf('reading 4D data: subject %s %d/%d ', substr, s , nsub)
    % data from different niftis is concatenated
    % read in all volumes to minimize io
    % also process the time-series
    % do this only once per subject and use the data for all nii
    try
        % read_data_sub_4D_concat also performs confound removal
	    subfiles = {};
        if isfield(NetFlags, 'subfile') && ~isempty(NetFlags.subfile)
            subfiles = NetFlags.subfile{s};
	        if ischar(subfiles)
                subfiles = {subfiles};
            end
        end
        [Vols, TSAll, TR, TSAllunfilt] = read_data_sub_4D_concat(NetFlags.sub{s}, NetFlags.SubDir{s}, XYZ(:,myIndAll), ...
                                         NetFlags, Vref, subfiles);
        datafield = matlab.lang.makeValidName(NetFlags.data{1});
        V1 = Vols.(datafield);
        V1 = V1(1);
        if ~all(V1.dim==Dref)
            fprintf('got dimensions')
            V1.dim
            fprintf('expected dimensions')
            Dref
            error('nifti dimensions mismatch')

        end
    catch ME
        getReport(ME)
        errstr = sprintf('ERR: %s data could not be read %s... skipping!\n', ME.message, NetFlags.sub{s});
        fprintf('ERR: %s!\n', errstr)
        if nsub==1
            rethrow(ME)
        else
            continue
        end
    end

    if ~ismember('alff', conn_metrics)
        TSAllunfilt = [];
    end

    elapsed = toc(timerSub);
    fprintf('%g seconds\n', elapsed)

    if isfield(NetFlags, 'subsampleTS')
        subTS = NetFlags.subsampleTS;
        if ~isnan(subTS) && subTS>0
            fprintf(' subsampling TS at %d ', subTS)
            its = 1:subTS:size(TSAll,1);
            TSAll = TSAll(its,:);
	    if ~isempty(TSAllunfilt)
	        TSAllunfilt = TSAllunfilt(its,:);
	    end
        end
    end

    [nTimePoints, nTimeSeries] = size(TSAll);
    haveunfilt = ~isempty(TSAllunfilt);

    fprintf('done\n')
    fprintf('conn: processing %d nifti\n', length(nii))
    for inii=1:length(nii) % for each nifti with VOIs
        timerNii = tic;
        % set properties for this nii
        niiname = niinames{inii};
		if savesubfile
			zfile = fullfile(niiprop.(niiname).destdir, [substr '_conn_' settings '.mat']);
            if isfile(zfile)
	            nifti_processed(inii) = true;
            end
		else
        	zfile = niiprop.(niiname).zfile;
		end
        Ind = niiprop.(niiname).Ind;
        VOI = niiprop.(niiname).VOI;
        uROI = unique(VOI);
        assert(~any(isnan(uROI)))
        assert(~any(uROI==0))
        nROI = numel(uROI);
        assert(length(VOI)==length(Ind), 'VOI and Ind lengths must match');

        fprintf('conn: subject %s %d/%d nifti %s ', NetFlags.sub{s}, s , nsub, niiname)

        if savesubfile
              conn = struct();
        else
            conn = conn_nii(niiname);
            if isfield(conn, substr)
                fprintf('nii processed...skipping!\n')
                continue
            end
        end

        timerPrep = tic;
        % adjust Ind and VOI to myIndAll
        [~,myInd] = ismember(Ind, myIndAll);
        %assert(sum(myInd==0)==0 && sum(isnan(myInd))==0) % this assumes all of Ind is in myIndAll which might not be the case with subject masks
        %assert(sum(isnan(myInd))==0)
        % remove 0s and NaNs
        assert(length(myInd)==length(Ind));
        VOI(myInd==0 | isnan(myInd)) = [];
        myInd(myInd==0 | isnan(myInd)) = [];
        % recreate VOI as myVOI
        myVOI = zeros(1, length(myIndAll));
        myVOI(myInd) = VOI;

        % get averaged time-series per ROI
        TSVOIAV = nan(nTimePoints, nROI);
        TSVOIEV = nan(nTimePoints, nROI);
        for ana=1:nROI
            roivox = find(myVOI==uROI(ana));
            if isempty(roivox)
            elseif length(roivox)>1
                TSVOIAV(:,ana) = nanmean(TSAll(:,roivox),2);
                TSVOIEV(:,ana) = ev(TSAll(:,roivox));
            else
                TSVOIAV(:,ana) = TSAll(:,roivox);
                TSVOIEV(:,ana) = TSAll(:,roivox);
            end
        end
        %TSVOIAV = TSVOIAV(:,var(TSVOIAV)>eps);

        % get WB time-series index
        [~,myIndWB] = ismember(IndWB, myIndAll);
        myIndWB(myIndWB==0) = [];
        myIndWB = myIndWB(var(TSAll(:,myIndWB))>eps);
        [~,myIndWBLeft] = ismember(IndWBLeft, myIndAll);
        myIndWBLeft(myIndWBLeft==0) = [];
        myIndWBLeft = myIndWBLeft(var(TSAll(:,myIndWBLeft))>eps);
        [~,myIndWBRight] = ismember(IndWBRight, myIndAll);
        myIndWBRight(myIndWBRight==0) = [];
        myIndWBRight = myIndWBRight(var(TSAll(:,myIndWBRight))>eps);

        % create anotehr layer of GM masking, by default it's inactive
        % but will be activated if individual gm mask is aksed for
        gm = ones(1, numel(myIndAll));
        if NetFlags.gotMask==2
            try
                Vgm = spm_vol(fullfile(NetFlags.datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},NetFlags.data{dt}, ['c1w' NetFlags.sub{s} '_meanEPI.nii']));
            catch
                Vgm = spm_vol(fullfile(cwd,'MaskenEtc','forRS',[spm_str_manip(NetFlags.Group,'rt') '.nii']));
            end
            % Kaustubh: not sure why the following is in mask space changing to data space
            % gm  = spm_get_data(Vgm,Vgm.mat \ VM.mat * [maskXYZ(:,Ind); ones(1,size(Ind,2))]);
            gm  = spm_get_data(Vgm,Vgm.mat \ V(1).mat * [XYZ(:,myIndAll); ones(1,length(myIndAll))]);
        end
        assert(length(gm)==length(myIndAll))

        fprintf(' (%d workers)', poolsize)
        fprintf(' (%d ROI)', nROI)
        fprintf(' (%d VOL)', nTimePoints)
        fprintf(' (%.3f sec)\n', toc(timerPrep))

        allroiprop = cell(1, nROI);

        % get dimensionality
        % using atanh seems to mess up the dimwithin which always ends up being size(y,2)-2
        % not using atanh for anything below
        for ana = 1:nROI
            roiprop = struct();
            xQ = find(myVOI==uROI(ana));
            xQ = xQ(gm(xQ)>=median(gm(xQ)));
            TSROI = TSAll(:, xQ);
            goodvar = var(TSROI)>eps;
            xQ = xQ(goodvar);
            TSROI = TSROI(:, goodvar);
            roiprop.('nvox') = length(xQ);
            roiprop.('ntime') = nTimePoints;
            roiprop.('nvoxWB') = length(myIndWB);
            roiXYZmm = XYZAllmm(:,xQ);
            if ~isempty(xQ)
                if ismember('fcnesswithin', conn_metrics)
                  roiprop.('eigvalWITHIN') = kp_eigval_mean0(TSROI);
                end

                if ismember('reho', conn_metrics)
                  roiprop.('reho') = KendallsW(TSROI);
                end

                if ismember('sampen', conn_metrics)
		            % using parameters from https://pubmed.ncbi.nlm.nih.gov/29676512/
		            mysampen = nan(1,size(TSROI,2));
		            for irep=1:size(TSROI,2)
		                mysampen(irep) = sampen(TSROI(:,irep), 2, 0.3, 'chebychev');
		            end
		            roiprop.('sampen') = nanmean(mysampen);
                end

		        if ismember('tsnr', conn_metrics)
                    roiprop.('tsnr') = nanmean(nanmean(TSROI)./nanstd(TSROI));
                end

                if ismember('reho10', conn_metrics) 
		            if size(TSROI,2)<=10
		                roiprop.('reho10') = NaN;
		            else
		                rhrep = nan(1,100);
		                for irep=1:100
		                    ivr = randsample(1:size(TSROI,2), 10);
		                    rhrep(irep) = KendallsW(TSROI(:,ivr));
		                end
		                roiprop.('reho10') = nanmean(rhrep);
		            end
                end % reho10

                if ismember('rehots', 'conn_metrics')
                   % time-shifted reho w.r.t. mean TS
                   tsmean = mean(TSROI, 2);
                   tsroi = TSROI;
                   crr = [];
                   for sf=1:10
                       tsm = tsmean;
                       tsm(end) = [];
                       tsroi(1,:) = [];
                       cr = [];
                       for ti=1:size(tsroi,2)
                           cr = [cr, corr(tsm, tsroi(:,ti))];
                       end
                       crr = [crr, mean(cr)];
                   end
                   roiprop.('rehots_mean') = nanmean(crr);
                   [crmx, crmxi] = nanmax(crr);
                   roiprop.('rehots_max') = crmx;
                   roiprop.('rehots_maxi') = crmxi;
                   clearvars tsroi
                end

	            % see: https://web.conn-toolbox.org/definitions/connectivity-measures/other
                % see: https://fcp-indi.github.io/docs/user/alff.html
                if ismember('alff', conn_metrics)
                  TSROIunfilt = TSAllunfilt(:, xQ);
                  alff = nan(1, size(TSROIunfilt,2));
                  falff = nan(1, size(TSROIunfilt,2));
                  for vi=1:size(TSROIunfilt,2)
                    alff(vi)  = nanstd(idealfilter_pass_col(TR, [0.01, 0.1], TSROIunfilt(:,vi)));
                    falff(vi) = alff(vi)/nanstd(TSROIunfilt(:,vi));
                  end
                  roiprop.('alff') = nanmean(alff);
                  roiprop.('falff') = nanmean(falff);
                  TSROIunfilt = [];
                end

                if isfield(NetFlags, 'conndim_saveTS') && NetFlags.conndim_saveTS
                    roiprop.('TS') = TSROI;
                end

                if ismember('fcd', conn_metrics) && ~ismember('fcness', conn_metrics) && ~isempty(myIndWB)
                    wbXYZmm = XYZAllmm(:,myIndWB);
                    z = corr(TSROI, TSAll(:,myIndWB));
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('FCdensity') = kp_FCdensity_corr(z, wbXYZmm, roiXYZmm);
                    end
                end

                if ~isempty(myIndWB) && ismember('fcness', conn_metrics)
                    wbXYZmm = XYZAllmm(:,myIndWB);
                    %z = atanh(corr_smooth_extreme(corr(TSAll(:,xQ), TSAll(:,myIndWB))));
                    z = corr(TSROI, TSAll(:,myIndWB));
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('eigvalWB') = kp_eigval_mean0(z, usegpu);
                        roiprop.('FCdensity') = kp_FCdensity_corr(z, wbXYZmm, roiXYZmm);
                        if isfield(NetFlags, 'conndim_randWB') && NetFlags.conndim_randWB
                            roiprop.('eigvalWB_rand') = [];
                            for rr=1:NetFlags.conndim_randWB
                                z = reshape(randsample(z(:), numel(z)), size(z));
                                roiprop.('eigvalWB_rand') = [roiprop.('eigvalWB_rand'), kp_eigval_mean0(z, usegpu)];
                            end
                        end
                        
                        if isfield(NetFlags, 'conndim_dist') && ~isempty(NetFlags.conndim_dist)
                            Dmm = pdist2(roiXYZmm', wbXYZmm', 'euclidean', 'Smallest', 1);
                            assert(length(Dmm)==size(z,2))
                            for di=1:length(NetFlags.conndim_dist)
                                distd = NetFlags.conndim_dist(di);
                                dists = num2str(distd);
                                imm = find(Dmm<=distd);
                                fprintf(' dist:%d #vox:%.3f ', distd, length(imm)/size(z,2))
                                zmm = z(:,imm);
                                roiprop.(['nvoxWB' dists 'mm']) = length(imm);
                                roiprop.(['eigvalWB' dists 'mm']) = kp_eigval_mean0(zmm, usegpu);
                                roiprop.(['FCdensity' dists 'mm']) = kp_FCdensity_corr(zmm, wbXYZmm(:,imm), roiXYZmm);
                                imm = find(Dmm>distd);
                                fprintf(' dist:%d #vox:%.3f ', distd, length(imm)/size(z,2))
                                zmm = z(:,imm);
                                roiprop.(['nvoxWB_' dists 'mm']) = length(imm);
                                roiprop.(['eigvalWB_' dists 'mm']) = kp_eigval_mean0(zmm, usegpu);
                                roiprop.(['FCdensity_' dists 'mm']) = kp_FCdensity_corr(zmm, wbXYZmm(:,imm), roiXYZmm);
                            end
                        end
                    end

                    % get hemispheric dimensionality
                    roiprop.('nvoxWBLeft') = length(myIndWBLeft);
                    z = corr(TSROI, TSAll(:,myIndWBLeft));
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('eigvalWBLeft') = kp_eigval_mean0(z, usegpu);
                    end

                    roiprop.('nvoxWBRight') = length(myIndWBRight);
                    z = corr(TSROI, TSAll(:,myIndWBRight));
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('eigvalWBRight') = kp_eigval_mean0(z, usegpu);
                    end
                end % myIndWB

                % corr with roi averages connections
                if size(TSVOIAV,2)
                    %z = atanh(corr_smooth_extreme(corr(TSAll(:, xQ), TSVOIAV)));
                    z = corr(TSROI, TSVOIAV);
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('eigvalROIAV') = kp_eigval_mean0(z);
                    end
                    z = corr(TSROI, TSVOIEV);
                    if ~isempty(z)
                        z(isnan(z)) = 0;
                        roiprop.('eigvalROIEV') = kp_eigval_mean0(z);
                    end
                end
            end % isempty xQ
            allroiprop{ana} = roiprop;
            fprintf('.')
        end % for ana=1:nROI

        %save info about this nifti
        conn.(substr).glocal = allroiprop;

        %if ~isempty(myIndWB)
        %    timerWBWB = tic;
        %    z = corr(TSAll(:,myIndWB), TSAll(:,myIndWB));
        %    z(isnan(z)) = 0;
        %    conn.(substr).eigvalWBWB = kp_eigval_mean0(z);
        %    fprintf('WBWB processed in %g seconds\n', toc(timerWBWB))
        %end

        % get pairwise connectivity
        if ismember('pairwise', conn_metrics)
            cr = corr(TSVOIAV);
            conn.(substr).roipair.mean = cr(tril(true(size(cr)),-1))';
            cr = corr(TSVOIEV);
            conn.(substr).roipair.ev   = cr(tril(true(size(cr)),-1))';
            if isfield(NetFlags, 'pairwise_saveTS') && NetFlags.pairwise_saveTS==true
                 fprintf('(pairwise_saveTS) ')
                 conn.(substr).roipair.meanTS = TSVOIAV;
                 conn.(substr).roipair.evTS = TSVOIEV;
            end
        end

        if ismember('pairwisespearman', conn_metrics)
            cr = corr(TSVOIAV, 'type', 'spearman');
            conn.(substr).roipair.mean_spearman = cr(tril(true(size(cr)),-1))';
            cr = corr(TSVOIEV, 'type', 'spearman');
            conn.(substr).roipair.ev_spearman   = cr(tril(true(size(cr)),-1))';
            if isfield(NetFlags, 'pairwise_saveTS') && NetFlags.pairwise_saveTS==true
                 fprintf('(pairwise_saveTS) ')
                 conn.(substr).roipair.meanTS = TSVOIAV;
                 conn.(substr).roipair.evTS = TSVOIEV;
            end
        end

        % partial
        if ismember('pairwisepartial', conn_metrics)
                fprintf('connectome: partial mean ')
                timepart = tic;
                cr = partialcorr(TSVOIAV);
                conn.(substr).roipair.meanpartial = cr(tril(true(size(cr)),-1))';
                fprintf(' %.3f sec\n', toc(timepart))
                fprintf('connectome: partial ev')
                timepart = tic;
                cr = partialcorr(TSVOIEV);
                conn.(substr).roipair.evpartial   = cr(tril(true(size(cr)),-1))';
                fprintf(' %.3f sec\n', toc(timepart))
        end

        if ismember('pairwiseLW', conn_metrics)
            try % in some cases cov1para returns all NaN
                cr = corrcov(TSVOIAV);
                conn.(substr).roipair.meanLW = cr(tril(true(size(cr)),-1));
                cr = corrcov(TSVOIEV);
                conn.(substr).roipair.evLW  = cr(tril(true(size(cr)),-1));
            catch MEME
            end
        end

        if ismember('edge', conn_metrics)
            ets = zscore(TSVOIAV);
            [ntime,nnodes] = size(ets);
            nedges = nnodes*(nnodes - 1)/2;
            [u,v] = find(triu(ones(nnodes),1));
            idx = (v - 1)*nnodes + u;
            ets = ets(:,u).*ets(:,v);
            cr = corr(ets);
            conn.(substr).edge.mean = cr(tril(true(size(cr)),-1))';
            % ev
            ets = zscore(TSVOIEV);
            [ntime,nnodes] = size(ets);
            nedges = nnodes*(nnodes - 1)/2;
            [u,v] = find(triu(ones(nnodes),1));
            idx = (v - 1)*nnodes + u;
            ets = ets(:,u).*ets(:,v);
            cr = corr(ets);
            conn.(substr).edge.ev = cr(tril(true(size(cr)),-1))';
        end

        if ismember('graph', conn_metrics)
            conn.(substr).roipair.mean_centrality5percent  = kp_conn2centrality(conn.(substr).roipair.mean, 0.05);
            conn.(substr).roipair.mean_centrality10percent = kp_conn2centrality(conn.(substr).roipair.mean, 0.1);
            conn.(substr).roipair.mean_centrality25percent = kp_conn2centrality(conn.(substr).roipair.mean, 0.25);
            conn.(substr).roipair.mean_centrality50percent = kp_conn2centrality(conn.(substr).roipair.mean, 0.50);
        end

        % https://www.pnas.org/content/early/2020/10/21/2005531117
        % https://github.com/brain-networks/edge-ts/blob/master/main.m
        if ismember('hacf', conn_metrics)
            conn.(substr).roipair.hacfmean = high_amplitude_cofluctuations(TSVOIAV);
            conn.(substr).roipair.hacfev = high_amplitude_cofluctuations(TSVOIEV);
        end % hafc

    fprintf('nifti %s processed in %.3f seconds\n', niiname, toc(timerNii))

    if savesubfile
        fprintf('saving file: %s\n', zfile)
        save(zfile, '-struct', 'conn');
    else
        conn_nii(niiname) = conn;
        % save regularly to avoid data loss
	    if s==1 || s==nsub || rand()<0.5
	  	    fprintf('saving file: %s\n', zfile)
            save(zfile, '-struct', 'conn');
	    end
    end % savesubfile
    end % nii
    fprintf('subject %s processed in %.3f seconds\n', substr, toc(timerSub))
end % sub

