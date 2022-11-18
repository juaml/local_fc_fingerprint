function Y = ev(y)


[m n]   = size(y);
if m > n
    [v s v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(spm_atranspa(y'));
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);
Y1 = Y-mean(Y);
