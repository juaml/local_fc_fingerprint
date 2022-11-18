function [hacf, rms] = high_amplitude_cofluctuations(ts, frackeep)
% https://www.pnas.org/content/early/2020/10/21/2005531117
% https://github.com/brain-networks/edge-ts
switch nargin
	case 0
		error('give me ts')
	case 1
		frackeep = [0.05, 0.1, 0.25, 0.5]; % default cutoff from the original code
end

[ntime,nnodes] = size(ts);
ts = zscore(ts);
[uu,vv] = find(triu(ones(nnodes),1));
idx = (vv - 1)*nnodes + uu;
% generate edge time series
ets = ts(:,uu).*ts(:,vv);
% calculate co-fluctuation amplitude at each frame
rms = sum(ets.^2,2).^0.5;
[~,idxsort] = sort(rms,'descend');

hacf = struct();
for i=1:length(frackeep)
    nkeep = round(ntime*frackeep(i));
    ikeep = idxsort(1:nkeep);
    frackeepName = matlab.lang.makeValidName(num2str(frackeep(i)));
    hacf.([frackeepName 'index']) = ikeep;
    if nkeep > 3
        cr = corr(ts(ikeep,:));
        hacf.(frackeepName) = cr(tril(true(size(cr)),-1))';
    end
end

