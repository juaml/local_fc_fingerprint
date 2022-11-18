function [datasetdir, timetaken, status, dropresult] = datalad_drop(files, datasetdir, params)
% returning datasetdir for uniform call to datalad_ functions
switch nargin
    case {0,1}
        error('need files and datasetdir')
    case 2
        params = '--nocheck';
end

if ischar(files)
    files = {files};
end

cwd = pwd;
status = nan(1, length(files));
tic;
for f=1:length(files)
    fprintf('datalad: drop %s\n', files{f})
    [filedir,filename,fileext] = fileparts(files{f});
    %[status, cmdout] = system(['chmod -R u+w ' filedir]);
    %if status>0
    %    cmdout
    %    error(['could not change directory permissions: ', filedir])
    %end
    %st = system('git config --add annex.pidlock true');
    if isfile(fullfile(datasetdir, files{f}))
        [status(f), dropresult] = system(['LD_LIBRARY_PATH= datalad -C ' datasetdir ' drop ' files{f} ' ' params]);
        if status(f)>0
            dropresult
            disp(['could not drop file: ', files{f}])
        end
    end
end
timetaken = toc;


