function look = lookup_reorder(look, subs)

looksub = look.sub;
if ~any(ismember(looksub, subs))
    looksub = matlab.lang.makeValidName(looksub);
end
assert(all(ismember(subs, looksub)))

[~,ii] = ismember(subs,looksub);

vars = fieldnames(look);
nsub = size(look.Cov,1);
%assert(nsub==length(subs))
for i=1:length(vars)
   x  = look.(vars{i});
   if isempty(x)
      continue
   end

   if size(x,1)==nsub
      x = x(ii,:);
   elseif numel(x)==nsub
      x = x(ii);
   end
   look.(vars{i}) = x;
end




