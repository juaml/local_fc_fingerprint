function y = idealfilter_pass_col(TR, Filt, y)
    % columns of data are assumed to be timeseries
    % TR: repetition time
    % Filt: filter frequency, two values
    % y : data matrix

if length(Filt) ~= 2
    return
end

if isnan(TR) or TR <= 0
    error('TR must be > 0')
end

m = mean(y);

y = timeseries(y, 0:TR:(TR*(size(y,1)-1)));
y = idealfilter(y, Filt, 'pass');
y = y.Data;

y = bsxfun(@plus,y,m);

