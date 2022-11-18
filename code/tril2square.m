function sq = tril2square(x, diagval)
% assuming diagonal is not included, i.e. tril(x,-1) was used

switch nargin
    case 0
        error('give me x')
    case 1
        diagval = 0;
end

n = length(x);
N = (sqrt(8*n+1)+1)/2;

sq = zeros(N,N);

i = find(tril(ones(N,N), -1));

sq(i) = x;
sq = sq + sq';

if diagval > 0
    sq = sq + diagval * eye(size(sq));
end





