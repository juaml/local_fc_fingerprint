addpath code

%% get the fingerprints

% create settings
settings = load('settings_fmriprep.mat');
settings.preproc_data = 'fmriprep_openneuro_res2'
settings.datadir = fullfile(pwd, 'data');
settings.spmdir = fullfile(pwd, 'spm12');
settings.lookup_file = fullfile(pwd, 'Lookups', 'ds000115_0back.txt');
settings.global_mask_to_use = fullfile(pwd, 'CAT12_IXI555_MNI152_TPM_GS_GMprob0.2_ds000115-fmriprep.nii');
settings.connectome_nii = {fullfile(pwd,'Power_3mm_ds000115-fmriprep.nii')};
settings.conn_metrics = {'pairwise', 'reho', 'alff'};
settings.outputdir = fullfile(pwd, 'output');
settings.datalad = 1;
settings.Group = 'ds000115_0back';
save('settings_fmriprep_0back.mat', '-struct', 'settings')

settings.lookup_file = fullfile(pwd, 'Lookups', 'ds000115_2back.txt');
settings.Group = 'ds000115_2back';
save('settings_fmriprep_2back.mat', '-struct', 'settings')

cd('code')
run_conn('../settings_fmriprep_0back.mat', 1:10)
run_conn('../settings_fmriprep_2back.mat', 1:10)

