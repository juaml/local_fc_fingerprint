function look = lookup_rm(look, torm)

if isempty(torm)
    return
end

vars = fieldnames(look);
nsub = size(look.Cov,1)
assert(max(torm)<=nsub)

for i=1:length(vars)
   x  = look.(vars{i});
   if isempty(x)
      continue
   end

   if size(x,1)==nsub
      x(torm,:) = [];
   elseif numel(x)==nsub
      x(torm) = [];
   end
   look.(vars{i}) = x;
end

