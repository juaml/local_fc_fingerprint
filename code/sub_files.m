% Kaustubh
% a function to get "relative" path of nifti and confound files
function [file4d, conffile, fextra, confdir, jsonfile] = sub_files(preproc_data, sub, SubDir, ...
                                                         RSdir, confdir, verbose)

if nargin<4
    error('give me 4 or more arguments!')
elseif nargin<5
    confdir=''; verbose = false;
elseif nargin<6
    verbose = false;
end

if exist(sub,'file')==2
    preproc_data = 'subisfile';
end

conffile = '';
jsonfile = '';

    switch lower(preproc_data)
        case 'subisfile'
	        file4d = sub;
        case {'rs_fix'}
            subsess = strsplit(sub, '/');
            assert(length(subsess)==2)
            file4d = fullfile(SubDir, subsess{1}, subsess{2}, ['w' subsess{1} '_' subsess{2} '.nii']);
        case {'normalised', 'norm', 'hcp_normalised'}
            file4d = fullfile(SubDir, sub, RSdir, 'Normalised',['w' sub '.nii']);
        case 'smooth_5mm'
            file4d = fullfile(SubDir, sub, RSdir, 'Smooth_5mm',['sw' sub '.nii']);
        case {'normalised_sub', 'norm_sub'} % directly under sub dir
            file4d = fullfile(RSdir, sub, ['w' sub '.nii']);
	    case {'cobre'}
	        file4d = fullfile(SubDir, sub, RSdir, ['w' sub '.nii']);
        case {'ucla'}
            file4d = fullfile(SubDir, RSdir, sub, ['w' sub '.nii']);
        case {'ucla_smooth'}
            file4d = fullfile(SubDir, RSdir, sub, ['sw' sub '.nii']);
        case 'bids'
            sess = strsplit(sub, '/');
            file4d = fullfile(SubDir, RSdir, sub, ['w' sess{1} '_' sess{2} '.nii']);
        case 'bids_smooth'
            sess = strsplit(sub, '/');
            file4d = fullfile(SubDir, RSdir, sub, ['sw' sess{1} '_' sess{2} '.nii']);
        case 'bids_cat'
            sess = strsplit(sub, '/');
            file4d = fullfile(SubDir, RSdir, sub, 'mri', ['m0wp1' sess{1} '_' sess{2} '.nii']);
        case {'hbn'}
            file4d = fullfile(SubDir, RSdir, sub, ['w' sub '.nii']);
	    case {'enki', 'enkibreath','enkichecker'}
            file4d = fullfile(SubDir, upper(RSdir), sub, 'ses-DS2', ['w' sub '_ses-DS2' '.nii']);
        case '1000rs'
            file4d = fullfile(SubDir, RSdir, sub, 'ses-1', ['w' sub '_ses-1' '.nii']);
        case 'cat'
            file4d = fullfile(sub, 'mri', ['m0wp1' sub '.nii']);
	    case 'cathcp'
            %% not implemented 
	    case {'hcp1', 'hcp1lr', 'hcp1lrfix'}
            scan    = 'rfMRI_REST1_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '_hp2000_clean.nii']);
        case {'hcp1nofix', 'hcp1lrnofix'}
            scan    = 'rfMRI_REST1_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
	    case 'hcp1lrs5'
            scan    = 'rfMRI_REST1_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, ['s5_' scan '_hp2000_clean.nii']);
        case {'hcp1rl', 'hcp1rlfix'}
            scan    = 'rfMRI_REST1_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '_hp2000_clean.nii']);
        case 'hcp1rlnofix'
            scan    = 'rfMRI_REST1_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
        case {'hcp2', 'hcp2lr', 'hcp2lrfix'}
            scan    = 'rfMRI_REST2_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '_hp2000_clean.nii']);
        case {'hcp2rl', 'hcp2rlfix'}
            scan    = 'rfMRI_REST2_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '_hp2000_clean.nii']);
        case {'hcp2nofix', 'hcp2lrnofix'}
            scan    = 'rfMRI_REST2_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
        case 'hcp2rlnofix'
            scan    = 'rfMRI_REST2_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
        case {'hcpsocial', 'hcpsociallr'}
            scan = 'tfMRI_SOCIAL_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
        case {'hcpwm', 'hcpwmlr'}
            scan = 'tfMRI_WM_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/', scan, [scan '.nii']);
        case {'hcpemotion', 'hcpemotionlr'}
            scan = 'tfMRI_EMOTION_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR.nii');
	    case {'hcprelational', 'hcprelationallr'}
            scan = 'tfMRI_RELATIONAL_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_RELATIONAL_LR/tfMRI_RELATIONAL_LR.nii');
	    case {'hcpmotor', 'hcpmotorlr'}
            scan = 'tfMRI_MOTOR_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR.nii');
	    case {'hcplanguage','hcplanguagelr'}
            scan = 'tfMRI_LANGUAGE_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR.nii');
        case {'hcpgambling', 'hcpgamblinglr'}
            scan = 'tfMRI_GAMBLING_LR';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_GAMBLING_LR/tfMRI_GAMBLING_LR.nii');
        case 'hcpsocialrl'
            scan = 'tfMRI_SOCIAL_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_SOCIAL_RL/tfMRI_SOCIAL_RL.nii');
        case 'hcpwmrl'
            scan = 'tfMRI_WM_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_WM_RL/tfMRI_WM_RL.nii');
        case 'hcpemotionrl'
            scan = 'tfMRI_EMOTION_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL.nii');
        case 'hcprelationalrl'
            scan = 'tfMRI_RELATIONAL_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_RELATIONAL_RL/tfMRI_RELATIONAL_RL.nii');
        case 'hcpmotorrl'
            scan = 'tfMRI_MOTOR_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_MOTOR_RL/tfMRI_MOTOR_RL.nii');
        case 'hcplanguagerl'
            scan = 'tfMRI_LANGUAGE_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_LANGUAGE_RL/tfMRI_LANGUAGE_RL.nii');
        case 'hcpgamblingrl'
            scan = 'tfMRI_GAMBLING_RL';
            file4d = fullfile(SubDir, sub, 'MNINonLinear/Results/tfMRI_GAMBLING_RL/tfMRI_GAMBLING_RL.nii');
        case 'macaque'
            file4d = fullfile(SubDir, sub, 'wmfunc.nii');
        case {'fmriprep'}
	        subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-preproc_bold.nii']);
        case {'fmriprep_openneuro_res2'}
	        subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii']);
        case {'fmriprepucla'} % func folder after sub and no desc at end, e.g. UCLA data from OpenNeuro
	        subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin2009cAsym_preproc.nii']);
        case {'fmripreparoma'}
            subsess = strsplit(sub, '/');
            assert(length(subsess) > 1)
            if length(subsess) == 2
                subsess{3} = '';
            end
            if length(subsess{3}) > 1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii']);
        case {'fmriprep2'} % version2 for fmriprep assumes func folder after session folder
	        subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, subsess{2}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-preproc_bold.nii']);
        case {'fmripreparoma2'}
            subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, subsess{2}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii']);
        case {'fmriprep-mni152nlin2009casym'}
	        subsess = strsplit(sub, '/');
            assert(length(subsess)>1)
            if length(subsess)==2
                subsess{3} = '';
            end
            if length(subsess{3})>1
                subsess{3} = ['_' subsess{3}];
            end
            file4d = fullfile(SubDir, subsess{1}, subsess{2}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii']);
        case 'msc' % midnight scan club preprocessed file
            subsess = strsplit(sub, '/');
            assert(length(subsess)==2)
            %derivatives/volume_pipeline/sub-MSC01/processed_restingstate_timecourses/ses-func01/talaraich/sub-MSC01_ses-func01_task-rest_bold_talaraich.nii.gz
            file4d = fullfile(SubDir,'volume_pipeline', subsess{1}, 'processed_restingstate_timecourses', subsess{2},'talaraich', [subsess{1} '_' subsess{2} '_task-rest_bold_talaraich.nii']);
        case 'hcp_ep'
            subsess = strsplit(sub, '/');
            sub = subsess{1};
            scan = subsess{2};
            assert(length(subsess)==2)
            file4d = fullfile(SubDir, sub, 'MNINonLinear', 'Results', scan, [scan '.nii.gz']);
        case 'hcp_epfix'
            subsess = strsplit(sub, '/');
            sub = subsess{1};
            scan = subsess{2};
            assert(length(subsess)==2)
            file4d = fullfile(SubDir, sub, 'MNINonLinear', 'Results', scan, [scan '_hp0_clean.nii.gz']);
        otherwise
            error(['Unknown data preproc: ' preproc_data])
        end % switch


    if verbose
        fprintf('file4d: %s\n', file4d);
    end

    % get the jsonfile
    switch lower(preproc_data)
        case {'fmriprep', 'fmripreparoma'}
            jsonfile = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-preproc_bold.json']);
        case {'fmriprep2', 'fmripreparoma2'}
            jsonfile = fullfile(SubDir, subsess{1}, subsess{2}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin6Asym_desc-preproc_bold.json']);
        case {'fmriprep-mni152nlin2009casym'}
            jsonfile = fullfile(SubDir, subsess{1}, subsess{2}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin2009cAsym_desc-preproc_bold.json']);
        case {'fmriprep_openneuro_res2'}
            jsonfile = [subsess{1} '_' subsess{2} subsess{3} '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.json'];
        otherwise
            jsonfile = '';
    end


    if regexp(lower(preproc_data), '^hcp')
        if regexp(lower(preproc_data), 'nofix$')
            preproc_data = 'hcpdatanofix';
        else
            preproc_data = 'hcpdata';
        end
    end

    fextra = {};

    if ~isempty(confdir)
        assert(isfolder(confdir))
        fprintf('(user confdir) ')
        % generate path relative to confdir
        switch lower(preproc_data)
            case {'hcpdata', 'hcpdatanofix'}
                scandir = fullfile(sub, 'MNINonLinear', 'Results', scan);
                if strcmp(lower(preproc_data), 'hcpdatanofix')
                    conffile = fullfile(scandir, ['Confounds_' sub '_noFIX.mat']);
                    if ~exist(fullfile(confdir, conffile), 'file')
                        conffile = fullfile(scandir, ['Confounds_' sub '_noFIX.tsv']);
                    end
                else
                    conffile = fullfile(scandir, ['Confounds_' sub '.mat']);
                    if ~exist(fullfile(confdir, conffile), 'file')
                        conffile = fullfile(scandir, ['Confounds_' sub '.tsv']);
                    end
                end
                fextra = {fullfile(SubDir, scandir, 'Movement_Regressors_dt.txt')};
            otherwise
                error(['read_data: ' preproc_data ' not implemented with user provided confoundsdir!'])
        end % switch

    else % confdir

      switch lower(preproc_data)
         case {'rs_fix'}
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' subsess{1} '_' subsess{2} '.mat']);
         case {'normalised', 'norm','smooth_5mm', 'hcp_normalised'}
             confdir = fileparts(file4d);
             confdir = fileparts(confdir);
             conffile = fullfile(confdir, ['Counfounds_' sub '.mat']);
         case {'enki', 'enkibreath','enkichecker'}
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sub '_ses-DS2' '.mat']);
         case '1000rs'
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sub '_ses-1' '.mat']);
         case {'cobre', 'ucla', 'ucla_smooth'}
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sub '.mat']);
         case {'cat','cathcp','bids_cat', 'msc'} % no confound file
             confdir = '';
             conffile = '';
         case {'bids', 'bids_smooth'}
             sess = strsplit(sub, '/');
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sess{1} '_' sess{2} '.mat']);
         case {'hcpdata'}
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sub '.mat']);
		     fextra = {fullfile(confdir, 'Movement_Regressors_dt.txt')};
        case {'hcpdatanofix'}
            confdir = fileparts(file4d);
            conffile = fullfile(confdir, ['Confounds_' sub '_noFIX.mat']);
            fextra = {fullfile(confdir, 'Movement_Regressors_dt.txt')};
        case 'macaque'
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Confounds_' sub '.mat']);
	        % create conf with regnames if it doesnt exist
	        if exist(conffile, 'file')~=2
                xx = load('/data/BnB_USER/xliu/Macaque_data//Macaque_preprocessing_Yerkes19/site_oxford/0001/Confounds_0001.mat');
                tmpconf = load(fullfile(confdir, ['Counfounds_' sub '.mat']));
                tmpconf.regnames = xx.regnames;
                fprintf('conf with regnames: %s\n', conffile) 
                save(conffile, '-struct', 'tmpconf');
             end
         case {'fmriprep', 'fmripreparoma', 'fmriprep2', 'fmripreparoma2', 'fmriprep-mni152nlin2009casym', 'fmriprep_openneuro_res2'}
             confdir = fileparts(file4d);
             conffile = {};
             conffile{1} = fullfile(confdir, [subsess{1} '_' subsess{2} subsess{3} '_desc-confounds_regressors.tsv']);
             conffile{2} = fullfile(confdir, [subsess{1} '_' subsess{2} subsess{3} '_desc-confounds_timeseries.tsv']);
             %conffile{3} = fullfile(confdir, [subsess{1} '_' subsess{2} '_' subsess{3} '_desc-confounds_regressors.tsv']);
             %conffile{4} = fullfile(confdir, [subsess{1} '_' subsess{2} '_' subsess{3} '_desc-confounds_timeseries.tsv']);
         case {'fmriprepucla'}
             confdir = fileparts(file4d);
             conffile = {};
             conffile{1} = fullfile(SubDir, subsess{1}, 'func', [subsess{1} '_' subsess{2} subsess{3} '_confounds.tsv']);
         otherwise
             confdir = fileparts(file4d);
             conffile = fullfile(confdir, ['Counfounds_' sub '.mat']);
      end % switch
      confdir = ''; % this will case use of datasetdir as the conffile now contains path path relative to datasetdir

    end % confdir empty


    if ~isempty(conffile) && isstr(conffile)
        conffile = {conffile};
    end

