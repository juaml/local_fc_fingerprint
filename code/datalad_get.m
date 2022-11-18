function [datasetdir, timetaken, status, getresult, getcmd] = datalad_get(files, datasetdir, datalad_args)
    %clonetmp, tmpdir, ntries, fail, verbose, loglevel)

switch nargin
    case {0,1}
        error('need files and datasetdir')
    case 2
        datalad_args = struct();
end

% default arguments
clonetmp = true; tmpdir = ''; ntries=1; fail=true; verbose=true; loglevel='info'; getsh='';
% user provided
f = fields(datalad_args);
for i=1:length(f)
  switch f{i}
    case 'clonetmp'
        clonetmp = datalad_args.clonetmp;
    case 'tmpdir'
        tmpdir = datalad_args.tmpdir;
    case 'ntries'
        ntries=datalad_args.ntries;
    case 'fail'
        fail =datalad_args.fail;
    case 'verbose'
        verbose = datalad_args.verbose;
    case 'loglevel'
        loglevel = datalad_args.loglevel;
    case 'getsh'
        getsh = datalad_args.getsh;
        getnsh = datalad_args.getnsh;
  end
end

assert(ntries>0)

cwd = pwd;
if verbose
    fprintf('[datalad_get] pwd: %s\n', cwd)
end

PATH = getenv('PATH');

if ~isempty(getsh)
    if verbose
        fprintf('[datalad_get] getsh: %s\n', getsh)
    end
    assert(isfile(getsh), [getsh ': must be a file'])
    [pth, fn, fe] = fileparts(getsh);
    assert(strcmpi(fe, '.sh'))
    if isempty(pth)
        getsh = ['./' getsh];
    end
end

if ischar(files)
    files = {files};
end

if isempty(tmpdir)
    tmpdir = fullfile(filesep, 'tmp');
end

if isempty(verbose)
    verbose=true;
end

if isempty(loglevel)
    loglevel='info';
end

if verbose
    fprintf('[datalad_get] LD_LIBRARY_PATH=%s\n', getenv('LD_LIBRARY_PATH'))
    fprintf('[datalad_get] MATLAB_SHELL=%s\n', getenv('MATLAB_SHELL'))
    [~, datalad_bin] = system('which datalad');
    fprintf('[datalad_get] datalad=%s\n', datalad_bin)
    [~, ga_bin] = system('which git-annex');
    %[~,ga_ldd] = system(['LD_LIBRARY_PATH=; ldd ' ga_bin]);
    %fprintf('[datalad_get] ldd %s\n%s\n', ga_bin, ga_ldd)
end

datasetdirorig = datasetdir;
if clonetmp
    [~, dsname] = fileparts(tempname);
    datasetdir2 = fullfile(tmpdir, tempname);
    cmdout = '';
    clonecmd = ['LD_LIBRARY_PATH=; datalad --log-level "' loglevel '" clone ' datasetdir ' ' datasetdir2];
    if verbose
        clonecmd
        status = system(clonecmd);
    else
        [status, cmdout] = system(clonecmd);
    end
    if status > 0
        error(['datalad clone failed: ' datasetdir ' -> ' datasetdir2])
    end
    fprintf('[datalad_get] cloned %s -> %s\n', datasetdir, datasetdir2)
    datasetdir = datasetdir2;
end

status = nan(1, length(files));
tic;
for f=1:length(files)
    fprintf('[datalad_get] get %s\n', files{f});
    [filedir, filename, fileext] = fileparts(files{f});

    % set environmental variables as matlab system does not fetch them from rc files
    % for AWS
    %setenv('DATALAD_hcp_s3_secret_id', getenv('DATALAD_hcp_s3_secret_id'))
    %setenv('DATALAD_hcp_s3_key_id', getenv('DATALAD_hcp_s3_key_id'))
    % for direct access to connectomedb
    %setenv('DATALAD_hcp_db_user', getenv('DATALAD_hcp_db_user'))
    %setenv('DATALAD_hcp_db_password', getenv('DATALAD_hcp_db_password'))
    if isempty(getsh)
        getncmd = ['LD_LIBRARY_PATH=; datalad -C ' datasetdir ' -l "' loglevel '" get -n ' filedir];
        getcmd = ['LD_LIBRARY_PATH=; datalad -C ' datasetdir ' -l "' loglevel '" get ' files{f}];
    else
        assert(isfile(getsh), 'datalad_getnsh must be a file')
        assert(isfile(getsh), 'datalad_getsh must be a file')
        getncmd = [getnsh ' ' datasetdir ' ' loglevel ' ' filedir];
        getcmd = [getsh ' ' datasetdir ' ' loglevel ' ' files{f}];
    end
    if verbose
        fprintf('getncmd=%s\n', getncmd)
        fprintf('getcmd=%s\n', getcmd)
    end

    % get dir contents
    getresult = '';
    if verbose
        status = system(getncmd);
    else
        [status, getresult] = system(getncmd);
    end
    if status > 0
        if fail
            error(['datalad could not get dir contents: ' filedir])
        else
            disp(['datalad could not get dir contents: ' filedir])
            continue
        end
    end

    % check if file exists
    % exist and isfile dont work here
    filecheckcmd = ['ls ' fullfile(datasetdir, files{f})];
    [a,b] = system(filecheckcmd);
    if a > 0
        fprintf('[datalad_get] (skip) %d file not there: %s\n', a, files{f})
        continue
    end

    for itry=1:ntries % datalad sometimes needs multiple tries
        getresult = '';
        if verbose
            status(f) = system(getcmd);
        else
            [status(f), getresult] = system(getcmd);
        end
        if status(f) > 0 % fail
            if itry==ntries
                getresult
                fprintf('datalad command: %s\n', getcmd)
                fprintf('original file: %s\n', fullfile(datasetdirorig, files{f}))
                if fail
                    error(['datalad could not get file: ' files{f}])
                else
                    disp(['datalad could not get file: ' files{f}])
                    continue
                end
            else
                fprintf('getresult=%s\n', getresult)
                fprintf('try %d of %d failed :(\n', itry, ntries)
            end
        else
            if verbose
                fprintf('try %d of %d successful :)\n', itry, ntries)
            end
            %if clonetmp
                %savecmd = ['LD_LIBRARY_PATH=; datalad -C ' datasetdir ' save'];
                %if verbose
                %    fprintf('saving dataset: %s\n', savecmd)
                %    status(f) = system(savecmd);
                %else
                %    [status(f), saveresult] = system(savecmd);
                %end
            %end
            break
        end
    end
end
timetaken = toc;


