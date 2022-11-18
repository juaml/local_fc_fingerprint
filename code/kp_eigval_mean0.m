function e = kp_eigval_mean0(data, usegpu)

if nargin<2
    usegpu = false;
end

if usegpu
    data = gpuArray(data);
end

[n,d] = size(data);
m = mean(data);
data = data - repmat(m, n, 1);
e = svd(data,0).^2/n;

if usegpu
    e = gather(e);
end

