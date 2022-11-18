function [y, TR, yunfilt] = time_series_preprocess(y,itime,conf,Qreg,badTP,Filt,data, regextra, detrends, ...
                                                   nvolrmstart, nvolkeepstart, volblocktime, volblockindex)
% itime is timepoints to use, if empty them use all
% conf is the confound structure
% Qreg is which regressors to use from conf.reg
% badTP
% Filt is filter frequencies, either empty or two values > 0
% data, data type to use
% regextra extra confounds
% detrends two element attay with [TS CONF] detrend polynomial degree
switch nargin
    case {0,1,2,3,4,5}
         error('give me arguments')
    case 6
         data = ''; regextra = []; detrends=false; nvolrmstart = 0; nvolkeepstart=0; volblocktime = []; volblockindex = [];
    case 7
         regextra = []; detrends=false; nvolrmstart = 0; nvolkeepstart=0; volblocktime = []; volblockindex = [];
    case 8
         detrends=false; nvolrmstart = 0; nvolkeepstart=0; volblocktime = []; volblockindex = [];
    case 9
         nvolrmstart = 0; nvolkeepstart=0; volblocktime = []; volblockindex = [];
    case 10
         nvolkeepstart=0; volblocktime = []; volblockindex = [];
    case 11
        volblocktime = []; volblockindex = [];
    case 12
        volblockindex = [];
end

if ischar(conf) % assuming mat file
    conf = load(conf);
end

TR = NaN;
if ~isempty(conf) && isfield(conf, 'TR')
    TR = conf.TR;
end

if isnan(TR)
    fprintf('(default TR) ')
    if ~isempty(strfind(upper(data),'HCP'))
        disp data
        TR = 0.7200;
    elseif ~isempty(strfind(upper(data),'IBC'))
        TR = 2.0;
    else
        error(sprintf('TR not known for data %s and TR not in confounds file!', data))
    end
end
if ischar(TR)  % this is loaded from confounds
    TR = str2num(TR);
end
fprintf('(TR %g) ', TR)

assert(length(itime)==size(y,1))

if length(detrends)==1
    detrends = repmat(detrends, 1, 2);
end

regnow = [];
if ~isempty(Qreg)
    reg = conf.reg(itime,Qreg);
    %ii = ~any(isnan(reg));
    %regnow = reg(:,ii);
    reg(isnan(reg)) =0;
    regnow = reg;
end

if badTP && isfield(conf, 'badTP')
    badTPs = zeros(size(conf.reg,1),1);
    badTPs(conf.badTP) = 1;
    regnow = [regnow, badTPs];
end

if ~isempty(regextra)
    regnow = [regnow, regextra];
end

%% keep only volumes asked for
itimeuse = true(1, length(itime));

% remove volumes from the start
if nvolrmstart>0
    itimeuse(1:nvolrmstart) = false;
end

% time based, e.g. block design
if ~isempty(volblocktime)
    assert(size(volblocktime,2)==2)
    assert(~isnan(TR))
    itimeTR = linspace(0, TR*length(itime), length(itime));
    if max(itimeTR) < max(volblocktime(:))
        error(sprintf('block timing problem: TR=%g, length=%g, asked=%g', TR, max(itimeTR), max(volblocktime(:))))
    end
    itimeusevol = false(1, length(itime));
    for i=1:size(volblocktime,1)
        itimeusevol = itimeusevol | (itimeTR >= volblocktime(i,1) & itimeTR <= volblocktime(i,2));
    end
    itimeuse = itimeuse & itimeusevol;
end

% index based, e.g. block design
if ~isempty(volblockindex)
    assert(size(volblockindex,2)==2)
    iindex = 1:length(itime);
    itimeusevol = false(1, length(itime));
    for i=1:size(volblockindex,1)
        itimeusevol = itimeusevol | (iindex >= volblockindex(i,1) & iindex <= volblockindex(i,2));
    end
    itimeuse = itimeuse & itimeusevol;
end

% keep only volumes at the start
if nvolkeepstart>0
    if nvolkeepstart>size(y,1)
        % something is wrong just return empty stuff
        y = [];
        TR = 0;
        yunfilt = [];
        return
    end

    itimeuse((nvolkeepstart+1):end) = false;
end

fprintf('(%d/%d vol) ', sum(itimeuse), length(itimeuse))
assert(sum(itimeuse)>0)

y = y(itimeuse,:);
if ~isempty(regnow)
    regnow = regnow(itimeuse,:);
end

poolsize = 0;
ppool = gcp('nocreate'); % If no pool, do not create new one.
if ~isempty(ppool)
    poolsize = ppool.NumWorkers;
end

if detrends(1)>0 && detrends(2)>0
    fprintf('(dtTS...')
    if poolsize>1
        parfor(ts=1:size(y,2), poolsize)
            dts = detrend(y(:,ts), detrends(1));
            y(:,ts) = dts;
        end
    else
        for ts=1:size(y,2)
            dts = detrend(y(:,ts), detrends(1));
            y(:,ts) = dts;
        end
    end
    detrends(1) = 0;
    fprintf('done) ')
end

% regress out confounds
if ~isempty(regnow)
    if detrends(2)>0
        for rc=1:size(regnow,2)
            try
                regnow(:,rc) = detrend(regnow(:,rc), detrends(2));
            catch MEME
            end
        end
        detrends(2) = 0;
        fprintf('(dtConf) ')
    end

    regnow = [regnow, ones(size(regnow,1),1)];
    bb = regnow\y;
    y = y-regnow*bb;
end

if detrends(1)>0
    fprintf('(dtTS...')
    if poolsize>1
        parfor(ts=1:size(y,2), poolsize)
            dts = detrend(y(:,ts), detrends(1));
            y(:,ts) = dts;
        end
    else
        for ts=1:size(y,2)
            dts = detrend(y(:,ts), detrends(1));
            y(:,ts) = dts;
        end
    end
    fprintf('done) ')
end

yunfilt = y;

if length(Filt)==2
    assert(~isnan(TR))
    assert(Filt(1)<Filt(2))
    y = idealfilter_pass_col(TR,Filt,y);
end

