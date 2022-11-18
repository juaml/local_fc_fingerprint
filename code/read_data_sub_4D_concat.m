function [Vols, TSAll, TR, TSAllunfilt, settings] = read_data_sub_4D_concat(sub, subdir, XYZ, NetFlags, Vref, subfiles)

switch nargin
    case {0,1,2,3,4}
        error('I take 5 or more arguments!')
    case 5
        subfiles = {};
end

cwd = pwd;

settings = kp_get_settings(NetFlags);

Filt = []; % default filter frequency: no filtering
if NetFlags.gotFilt==1 % BP
    Filt = [0.01, 0.08];
elseif NetFlags.gotFilt==2 % HP
    Filt = [0.01, inf];
end

% override Filt if provided directly with NetFlags
if isfield(NetFlags, 'filterFrequency') % this can be empty
    Filt = NetFlags.filterFrequency;
    assert(isempty(Filt) || length(Filt)==2, 'filterFrequency must be empty or have 2 element')
else
    error('must have NetFlags.filterFrequency!!!')
end

if length(Filt)==2
    assert(Filt(1)<Filt(2), 'filterFrequency(1) must be < filterFrequency(2)')
    fprintf(' (filt:%g-%g) ', Filt(1), Filt(2))
else
    fprintf(' (filt:none) ')
end

regress_confounds = true;
if isfield(NetFlags, 'regress_confounds') && ~isempty(NetFlags.regress_confounds)
    regress_confounds = NetFlags.regress_confounds;
end

confounds_extra = true;
if isfield(NetFlags, 'confounds_extra') && ~isempty(NetFlags.confounds_extra)
    confounds_extra = NetFlags.confounds_extra;
end

badTP = false;
if isfield(NetFlags, 'confound_badTP')
    badTP = NetFlags.confound_badTP;
end

useSubFile = false;
if length(subfiles)==length(NetFlags.data)
    useSubFile = true;
end

nvolrmstart = 0;
if isfield(NetFlags, 'nvolrmstart') && ~isempty(NetFlags.nvolrmstart)
    nvolrmstart = NetFlags.nvolrmstart;
end

nvolkeepstart = 0;
if isfield(NetFlags, 'nvolkeepstart') && ~isempty(NetFlags.nvolkeepstart)
    nvolkeepstart = NetFlags.nvolkeepstart;
end

volblocktime = [];
if isfield(NetFlags, 'volblocktime')
    volblocktime = NetFlags.volblocktime;
end

volblockindex = [];
if isfield(NetFlags, 'volblockindex')
    volblockindex = NetFlags.volblockindex;
end

confoundsdir = '';
if isfield(NetFlags, 'confoundsdir') && ~isempty(NetFlags.confoundsdir)
    confoundsdir = NetFlags.confoundsdir;
    assert(isfolder(confoundsdir), 'confoundsdir must be a folder')
end

datalad = false;
if isfield(NetFlags, 'datalad') && ~isempty(NetFlags.datalad)
    datalad = NetFlags.datalad;
end

dataladdrop = true;
if datalad && isfield(NetFlags, 'dataladdrop') && ~isempty(NetFlags.dataladdrop)
    dataladdrop = NetFlags.dataladdrop;
end

dataladclonetmp = true;
if datalad && isfield(NetFlags, 'datalad_clonetmp') && ~isempty(NetFlags.datalad_clonetmp)
    dataladclonetmp  = NetFlags.datalad_clonetmp;
end

expectconffile = true;
if isfield(NetFlags, 'expectconffile')
    expectconffile = NetFlags.expectconffile;
end

files4d = {};
Vols = struct();
TSAll = [];
TSAllunfilt = [];
try
    for dt=1:numel(NetFlags.data)
        datafield = matlab.lang.makeValidName(NetFlags.data{dt});
        % load the subject volume and confounds
        if useSubFile
            sub = fullfile(NetFlags.datadir, subfiles{dt});
            fprintf('using settings provided subject file: %s\n', sub)
        end
        [files4d{dt}, conffile, fextra, confdir, jsonfile] = sub_files(NetFlags.preproc_data, sub, subdir, NetFlags.data{dt} , confoundsdir);
        files4d{dt}
        if ~isempty(conffile) && isstr(conffile)
            conffile = {conffile};
        end

        % read in data
        ladfiles = {};
        % explicitly add '.gz' file
        ladfiles = [files4d{dt}, [files4d{dt} '.gz'], fextra];

        % check if the confound file is part of the datalad dataset
        conffileindatlad = false;
        ladfilesconf = [];
        if datalad && expectconffile && ~isempty(conffile)
            conffilefull = fullfile(confdir, conffile);
            if isstr(conffilefull)
                conffilefull = {conffilefull};
            end
            for icf=1:length(conffilefull)
                [conffiledir,conffilename,conffileext] = fileparts(conffilefull{icf});
                niidir = fileparts(files4d{dt});
                if strcmp(conffiledir, niidir) && datalad % also get conffile
                    conffilelad = strrep(conffilefull, NetFlags.datadir, '');
                    conffilelad = regexprep(conffilelad, ['^' filesep], '');
                    ladfilesconf = [ladfilesconf, conffilelad];
                    conffileindatlad = true;
                    conffile{icf} = [conffilename, conffileext];
                end
            end
        end

        if datalad
            if conffileindatlad
                ladfiles = [ladfiles, ladfilesconf];
            end
            % get additional datalad options
            datalad_args = struct();
            datalad_args.fail = false;
            if isfield(NetFlags, 'datalad_verbose')
                datalad_args.verbose = NetFlags.datalad_verbose;
            end
            if isfield(NetFlags, 'datalad_loglevel')
                datalad_args.loglevel = NetFlags.datalad_loglevel;
            end
            if isfield(NetFlags, 'datalad_ntries')
                datalad_args.ntries = NetFlags.datalad_ntries;
            end
            if isfield(NetFlags, 'datalad_fail')
                datalad_args.fail = NetFlags.datalad_fail;
            end
            if isfield(NetFlags, 'datalad_getsh')
                datalad_args.getsh = NetFlags.datalad_getsh;
            end
            [datasetdir, timetaken, status, cmdout] = datalad_get(ladfiles, NetFlags.datadir, datalad_args);
            fprintf('datalad: get %.4f seconds\n', timetaken)
        else
            datasetdir = NetFlags.datadir;
        end

        file4dfull = fullfile(datasetdir, files4d{dt});
        if conffileindatlad
            confdir = fileparts(file4dfull);
        end

        if ~isfile(file4dfull)
            file4dfull = [file4dfull '.gz'];
        end
        assert(isfile(file4dfull), ['nifti file (.nii/.nii.gz) does not exist: ' file4dfull])
        Vols.(datafield) = spm_vol(file4dfull);


        % we have all the volumes ready here
        if ~isempty(Vref)
            V1 = Vols.(datafield);
            V1 = V1(1);
            if ~isequal(Vref.mat, V1.mat)
                fprintf('ERROR: the reference affine different than subject affine\n')
                fprintf('reference affine\n')
                Vref.mat
                fprintf('subject affine')
                V1.mat
                error('ERROR: the reference affine different than subject affine')
            end
        end

        if nargout>1
            [y,itime,ncorrupt] = spm_get_data_corrupt(Vols.(datafield),XYZ);
            if ncorrupt>0
                fprintf('(%d corrupted vol) ', ncorrupt);
            end
            Qreg = [];
            Nreg = [];
            conf  = [];
            regextra = [];
            if regress_confounds
                % regress out the confounds from all time-series
                % reg is loaded from confounds above
                if isempty(confdir)
                    confdir = datasetdir;
                end

                if iscell(conffile) % check which one exists and take it
                    conffilefound = false;
                    for icf=1:length(conffile)
                        confdir
                        conffile{icf}
                        conffilefull = fullfile(confdir, conffile{icf})
                        fprintf('checking conffile: %s', conffilefull)
                        if isfile(conffilefull)
                            conffile = conffile{icf};
                            fprintf(' exists!\n')
                            conffilefound = true;
                            break
                        else
                            fprintf(' does not exist!\n')
                        end
                    end
                    assert(conffilefound==true, 'confound file not found')
                end

                conffilefull = fullfile(confdir, conffile);
                fprintf('confounds file: %s\n', conffilefull)
                [~,~,conffileext] = fileparts(conffile);

                if strcmpi(conffileext, '.mat')
                    conf = load(conffilefull);
                    if isfield(NetFlags, 'ConfoundsTS')
                        if ~isempty(NetFlags.ConfoundsTS)
                            if ~isfield(conf, 'regnames')
                                fprintf('your confound mat file has no regnames: %s (going old school) ', conffile)
                                %error('confound mat file has no regnames')
                                Qreg = confound_indices(NetFlags.gotPCA, NetFlags.gotGlob, NetFlags.data{dt});
                            else
                                [Qreg, Nreg] = confound_indices_TS(NetFlags.ConfoundsTS, conf.regnames);
                                if isfield(NetFlags, 'ConfoundsTSderiv') && ~isempty(NetFlags.ConfoundsTSderiv)
                                    fprintf('derivatives: ')
                                    [conf, Qreg, Nreg, nDerivs] = confound_add_derivatives(conf, Qreg, Nreg, NetFlags.ConfoundsTSderiv);
				                    fprintf('%d added\n', nDerivs)
                                end
                            end
                            fprintf(' (ConfoundsTS:%s (%d)) ', NetFlags.ConfoundsTS, length(Qreg));
                        end
                    end
                elseif strcmpi(conffileext, '.tsv') || strcmpi(conffileext, '.csv') % must be a file with headers
                        confTable = readtable(fullfile(confdir, conffile), 'FileType', 'text');
                        conf = [];
                        vars = confTable.Properties.VariableNames;
                        % make all columns double
                        for ii=1:length(vars)
                            if iscell(confTable.(vars{ii}))
                                confTable.(vars{ii}) = str2double(confTable.(vars{ii}));
                            end
                        end
                        % workaround: create a conf structure that we expect later
                        conf.regnames = confTable.Properties.VariableNames;
                        conf.reg = table2array(confTable);
                        assert(length(conf.regnames) == size(conf.reg,2), 'all confounds must me named')
                        if isfield(NetFlags, 'ConfoundsTS')
                            [Qreg, Nreg] = confound_indices_TS(NetFlags.ConfoundsTS, conf.regnames);
                            if isfield(NetFlags, 'ConfoundsTSderiv') && ~isempty(NetFlags.ConfoundsTSderiv)
                                fprintf('derivatives: ')
                                [conf, Qreg, Nreg, nDerivs] = confound_add_derivatives(conf, Qreg, Nreg, NetFlags.ConfoundsTSderiv);
                                fprintf('%d added\n', nDerivs)
                            end
                        end
                elseif strcmpi(conffileext, '.txt')
                    conf.reg = textread(fullfile(confdir, conffile));
                    Qreg = 1:size(conf.reg,2);
                else
                    if expectconffile
                        error(['unknown confound file type: ', conffileext])
                    else
                        conf = [];
                        Qreg = [];
                    end
                end

                regextra = [];
                if ~isempty(fextra) && confounds_extra
                    regextra = read_reg(fextra, datasetdir);
                    fprintf(' (%d extra confounds) ', size(regextra,2))
                end
            end % regress_confounds

            detrendTS = false;
            if isfield(NetFlags, 'detrendTS') && ~isempty(NetFlags.detrendTS)
                detrendTS = NetFlags.detrendTS;
            end

            detrendConfTS = false;
            if isfield(NetFlags, 'detrendConfTS') && ~isempty(NetFlags.detrendConfTS)
                detrendConfTS = NetFlags.detrendConfTS;
            end

            detrends = [detrendTS, detrendConfTS];
            if any(detrends>0)
                fprintf(' (detrend %g %g)', detrends(1), detrends(2))
            end

            % get the TR, user provided value takes precedence
            if isfield(NetFlags, 'TR')
                conf.TR = NetFlags.TR;
                assert(conf.TR > 0, 'user TR must be > 0')
                % get the TR from jsonfile
            elseif ~isempty(jsonfile)
                fprintf('\nTR from json: %s\n', fullfile(confdir, jsonfile))
                fid = fopen(fullfile(confdir, jsonfile));
                raw = fread(fid,inf);
                str = char(raw');
                fclose(fid);
                json = jsondecode(str);
                conf.TR = json.RepetitionTime;
                assert(conf.TR > 0, 'json TR must be > 0')
            else
                fprintf(' (TR?)')
            end

            [y, TR, yunfilt] = time_series_preprocess(y,itime,conf,Qreg,badTP,Filt,NetFlags.data{dt},regextra, ...
                                detrends, nvolrmstart,nvolkeepstart, volblocktime, volblockindex);
            if isempty(y)
                error(' (could not process time-series) ')
            end

            if isfield(NetFlags, 'zscoreTS') && ~isempty(NetFlags.zscoreTS) && NetFlags.zscoreTS
                y = zscore(y);
                %yunfilt = zscore(yunfilt); % do not zscore unfilt as its used for alff/falff
            end

	        if nargout<=3
		        yunfilt = [];
	        end

            TSAll = [TSAll; y];
            TSAllunfilt = [TSAllunfilt; yunfilt];
            clearvars y, yunfilt;
        end
    end
catch ME
    getReport(ME)
    [~, fname] = fileparts(mfilename('fullpath'));
    fprintf('ERR: %s %s data could not be read... skipping!\n', fname, sub)
    rethrow(ME)
end

% remove datalad files
if datalad && dataladdrop
    [datasetdir, timetaken, status, cmdout] = datalad_drop(ladfiles, datasetdir);
    fprintf('datalad: drop %.4f minutes\n', timetaken/60)
end


