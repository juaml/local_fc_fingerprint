function [look, torm] = lookup_filter(look, sub)

looksub = look.sub;

if ~iscell(sub)
    sub = looksub(sub);
end

if ~any(ismember(looksub, sub))
    looksub = matlab.lang.makeValidName(looksub);
end

[~,ii] =  ismember(sub, looksub);
torm = sort(setdiff(1:length(looksub), ii));

look = lookup_rm(look, torm); 
look = lookup_reorder(look, sub);

