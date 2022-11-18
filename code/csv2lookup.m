function look = csv2lookup(csv)
% convert csv file to a lookup structure

csv = readtable(csv, 'ReadVariableNames', true, 'FileType', 'text', 'Delimiter', '\t');
vars = csv.Properties.VariableNames;

assert(ismember('sub', vars))
assert(ismember('SubDir', vars))

look = [];
if iscell(csv.sub)
    look.sub = csv.sub;
else
    look.sub = cellstr(num2str(csv.sub));
end

if iscell(csv.SubDir)
    look.SubDir = csv.SubDir;
else
    look.SubDir = cellstr(num2str(csv.SubDir));
end

if ismember('hasRS', vars)
    look.hasRS = csv.hasRS;
else
    look.hasRS =  ones(1, length(look.sub));
end

if ismember('isPat', vars)
    look.isPat = csv.isPat;
else
    look.isPat =  zeros(1, length(look.sub));
end

covs = setdiff(vars, {'sub', 'SubDir', 'hasRS', 'isPat'});
if isempty(covs)
   look.Cov = ones(length(look.sub),1);
   look.Covariates = {'ones'};
else
  look.Cov = table2array(csv(:,covs));
  look.Covariates = vars;
end

