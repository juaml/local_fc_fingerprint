function [Ind, VOI, XYZ, XYZmm, Vnii] = nii2voiind(nii, NetFlags, V, sub)
% nii: nifti file location as string
% NetFlags: struct with settings
% sub: which sub to use as reference, defaults to 1
% V: spm vol as reference, if provides takes preferences over sub
% returned XYZ and XYZmm are for whole 3D volume and not only for roi

switch nargin
    case {0,1,2}
        error('give me arguments!')
    case 3
        sub = 1;
end

nativespace = false;
if isfield(NetFlags, 'nativespace')
    nativespace = NetFlags.nativespace;
end
if isnan(nativespace) || isempty(nativespace)
    nativespace = false;
end

if isempty(V) % get reference nifti Volume for one subject
    try
      %[V,~,~,~,file4d] = read_data_sub_4D(NetFlags.preproc_data, NetFlags.sub{sub}, NetFlags.SubDir{sub}, NetFlags.datadir, NetFlags.data{1});
      file4d = sub_files(NetFlags.preproc_data, NetFlags.sub{sub}, NetFlags.SubDir{sub}, NetFlags.data{1})
      V = spm_vol(fullfile(NetFlags.datadir, file4d));
      V = V(1);
    catch ME
      if ~nativespace
        fprintf('reading sub %d failed\n', sub)
      else
        fprintf('reading sub %d failed, wont continue as we want to use nativespace (not implemented)\n', sub)
        rethrow(ME)
      end
    end
end

% get the ROI nifti
if ~isfile(nii)
    fprintf('expected a file but this one is not: %s\n', nii)
    assert(isfile(nii))
end
VniiV = spm_vol(nii);
Vnii = spm_read_vols(VniiV);
Vniisize = size(Vnii);
Vniid = length(Vniisize);
Vniisize3D = Vniisize(1:3);
Vnii(isnan(Vnii)) = 0;

if nativespace
   error('nativespace not implemented yet :(')
else
    if ~(all(V.dim==Vniisize3D)) % make sure that the mask and subjects are in the same space
        if ~exist('file4d', 'var') || isempty(file4d)
            file4d = 'argument_provided';
        end
        fprintf('subject file: %s\n', file4d)
        V.dim
        fprintf('ROI file: %s\n', nii)
        size(Vnii)
        error('ERROR: reference and nii dimensions mismatch')
    end
    % check the affine
    if ~isequal(V.mat, VniiV.mat)
        fprintf('ERROR: reference and nii affine mismatch\n')
        fprintf('reference affine\n')
        V.mat
        fprintf('nii affine\n')
        Vnii.mat
        error('ERROR: reference and nii affine mismatch')
    end
end

mask = [];
if isfield(NetFlags, 'global_mask_to_use') && ~isempty(NetFlags.global_mask_to_use)
  [mask, maskfile] = get_mask(NetFlags.global_mask_to_use);
  fprintf('using global mask: %s\n', maskfile)
end

% get all indices
clearvars XYZ;
[XYZ(:,1), XYZ(:,2), XYZ(:,3)] = ind2sub(V(1).dim(1:3), find(ones(V(1).dim(1:3))));
XYZ = XYZ';
XYZmm = V(1).mat * [XYZ; ones(1,size(XYZ,2))];
XYZmm = XYZmm(1:3,:);

Ind = []; VOI = [];
if Vniid==3
    uROI = unique(Vnii(:));
    uROI(uROI==0 | isnan(uROI)) = [];
    nROI = length(uROI);
    fprintf('%d ROIs in file %s\n', nROI, nii)
    if ~isempty(mask)
        Vnii = mask .* Vnii;
    end
    for rr=1:nROI
        roi = uROI(rr);
        roi_vox = find(Vnii==roi);
        if isempty(roi_vox)
            roi_vox = [NaN];
        end
        Ind = [Ind roi_vox'];
        VOI = [VOI repmat(rr,1,length(roi_vox))];
    end
elseif Vniid==4
    nROI = size(Vnii,4);
    fprintf('%d ROIs in file %s\n', nROI, nii)
    for rr=1:nROI
        if ~isempty(mask)
            Vnii(:,:,:,rr) = mask .* Vnii(:,:,:,rr);
        end
        roi_vox = find(Vnii(:,:,:,rr)>0);
        if isempty(roi_vox)
            roi_vox = [NaN];
        end
        Ind = [Ind roi_vox'];
        VOI = [VOI repmat(rr,1,length(roi_vox))];
    end
else
    error(sprintf('the nii volume must be either 3D or 4D, got %dD', Vniid))
end

assert(length(Ind)==length(VOI))
assert(size(XYZ,2)>=max(Ind))

if length(unique(VOI))<nROI
    fprintf('Warning: some ROIs ended up empty: %d remaining out of %d\n', length(unique(VOI)),nROI)
end

