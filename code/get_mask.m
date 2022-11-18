function [mask, maskfile, VM] = get_mask(maskname, maskdir, V)

switch nargin
    case 0
        error('need maskname')
    case 1
        maskdir = fullfile('.','MaskenEtc'); V = [];
    case 2
        V = []; % reference volume
end

mask = [];
VM = [];

fprintf('checking mask: %s\n', maskname)

if isempty(maskdir)
  maskdir = '.';
end

havefile = false;
[maskd,maskf,maske] = fileparts(maskname);
if ismember(maske, {'.nii', '.nii.gz', '.mat'}) && isfile(maskname)
    havefile = true;
    maskfile = maskname;
    fprintf('found mask as a file!\n')
end

if ~havefile
  % the default files are 2x2x2 considering RS/task data
  switch maskname
    case 'MNI152GM'
        maskfile = fullfile(maskdir, 'MNI152GM.nii');
    case 'CanLabGM'
        maskfile = fullfile(maskdir, 'CanLabGM.nii');
    case 'FSL04'
        maskfile = fullfile(maskdir, 'FSL_MNI152_GM04.nii');
    case 'FSL025'
        maskfile = fullfile(maskdir, 'FSL_MNI152_GM025.nii');
    case 'CAT12_02'
        if isempty(V)
            maskfile = fullfile(maskdir, 'CAT12_IXI555_MNI152_TPM_GS_GMprob0.2_clean_2x2x2.nii');
        else
            foundSubDim = false;
            fs = dir(fullfile(maskdir,'CAT12_*prob0.2*.nii'));
            for ii=1:length(fs)
                maskfile = fullfile(maskdir,fs(ii).name);
                VM = spm_vol(maskfile);
                if all(V.dim==VM.dim)
                    foundSubDim = true;
                    break
                end
            end
            if ~foundSubDim
                error('could not find CAT12_02 mask matching subject dimentions')
            end
        end
    case 'CAT12_02_vbm'
        maskfile = fullfile(maskdir, 'CAT12_IXI555_MNI152_TPM_GS_GMprob0.2_clean.nii');
    case 'MacaqueYerkes99'
        maskfile = fullfile(maskdir,'forRS', 'Macaque_Yerkes19_Template_GM99.nii');
    case {'none', 'None', 'NONE'}
        maskfile = 'none';
    otherwise
        error(['unknown mask: ' maskname])
  end
end % havefile

if ~isempty(maskfile) && ~strcmpi(maskfile, 'none')
    [~,~,maske] = fileparts(maskfile);
    switch maske
        case '.mat'
            load(maskfile)
            mask = zeros(V(1).dim(1:3));
            ind = sub2ind(V(1).dim(1:3), maskXYZ(1,:), maskXYZ(2,:), maskXYZ(3,:));
            mask(ind) = 1;
        case {'.nii', '.nii.gz'}
            VM = spm_vol(maskfile);
            if ~isempty(V)
                assert(all(VM.dim==V.dim))
            end
	        mask = spm_read_vols(VM);
	        mask(mask>0) = 1;
        otherwise
            error(['unknown maskfile: ' maskfile])
      end
end

