function [Ind, CATnii, CATniiV, th] = subject_gm_voxels(NetFlags, sub, V, th)
% NetFlags: struct with settings
% sub: which sub to use as reference, defaults to 1
% V: spm vol as reference, if provides takes preferences over sub

switch nargin
    case {0,1}
        error('give me arguments!')
    case 2
        V = []; th=0;
    case 3
        th=0;
end

assert(isfield(NetFlags, 'CATdir') && exist(NetFlags.CATdir, 'dir')==7)

if isempty(V)
    try
      [V,~,~,~,file4d] = read_data_sub_4D(NetFlags.preproc_data, NetFlags.sub{sub}, NetFlags.SubDir{sub}, NetFlags.datadir, NetFlags.data{1});
      V = V(1);
    catch ME
        fprintf('reading sub %d failed, wont continue as we want to use nativespace\n', sub)
        rethrow(ME)
    end
else
   if ischar(V)
      file4d = V;
      V = spm_vol(file4d);
   else
      file4d = V.fname;
   end
end

assert(exist(file4d, 'file')==2)

% get the CAT nifti

CATnii = fullfile(NetFlags.CATdir, NetFlags.sub{sub}, 'mri', ['wp1' NetFlags.sub{sub} '.nii']);
CATniiV = spm_vol(CATnii);

if ~all(V.dim==CATniiV.dim)
    %tmpnii = [tempname '.nii'];
    %voxsize = abs(diag(V.mat));
    %voxsize = voxsize(1:3)';
    %reslice_nii(CATnii, tmpnii, voxsize);
    tmpnii = reslice_spm(file4d, CATnii);
    CATnii = tmpnii;
    CATniiV = spm_vol(CATnii);
    D = spm_read_vols(CATniiV);
    delete(tmpnii);
else
    D = spm_read_vols(CATniiV);  
end

if ~isempty(V)
    assert(all(V.dim==size(D)))
end

Ind = find(D>th);


