function rtmpf = reslice_spm(refimg, targetimg, interp, clean)
switch nargin
   case {0,1}
       error('give me ref and target!')
   case 2
       interp = 4; clean=true;
   case 3
       clean=true;
end


[tmpd, tmpn] = fileparts(tempname);
tmpf = fullfile(tmpd, [tmpn '.nii']);
rtmpf = fullfile(tmpd, ['r' tmpn '.nii']);

try
copyfile(targetimg, tmpf);
catch ME
  fprintf('copying failed: %s -> %s\n', targetimg, tmpf)
  rethrow(ME);
end

targetimg = tmpf;

matlabbatch{1}.spm.spatial.coreg.write.ref = {[refimg,',1']};
matlabbatch{1}.spm.spatial.coreg.write.source = {[targetimg, ',1']};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

if clean
    delete(tmpf);
end

