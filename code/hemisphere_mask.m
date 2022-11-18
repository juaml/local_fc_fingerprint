function [hemi] = hemisphere_mask(V, nifti)
% V: spm_vol
% nifti: string with file name
% returns a mask with three regions 1:left, 2:exact center and 3:right

switch nargin
    case 0
        error('give me a som_vol')
    case 1
        nifti = '';
end

hemi = ones(V.dim(1:3));
[XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(size(hemi),find(hemi));
XYZ = XYZ';
XYZmm = V.mat * [XYZ; ones(1,size(XYZ,2))];
XYZmm = XYZmm(1:3,:);

hemi(:) = 0;
hemi(XYZmm(1,:)<0)  = 1;
hemi(XYZmm(1,:)==0) = 2;
hemi(XYZmm(1,:)>0)  = 3;

if ~isempty(nifti)
    V.fname = nifti;
    V.provate.fname = nifti;
    [~] = spm_write_vol(V, hemi);
end


