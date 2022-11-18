% 1nn sequencial accuracy
function [results, x1, x2] = kp_identification(x1, x2, ks, conn, distmetric, patt, conf, sub1, sub2, subset, feat)
% supports multisession data by averaging them
switch nargin
    case {0,1}
       error('need two matrices/files!')
    case 2
       ks = 1; conn = 'connm'; distmetric = 'correlation';  patt = ''; conf = []; sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 3
       conn = 'connm'; distmetric = 'correlation';  patt = ''; conf = []; sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 4
       distmetric = 'correlation';  patt = ''; conf = []; sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 5
         patt = ''; conf = []; sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 6
        conf = []; sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 7
        sub1 = {}; sub2 = {}; subset=0; feat = [];
    case 8
       sub2 = {}; subset=0; feat = [];
    case 9
        subset=0; feat = [];
    case 10
        feat = [];
end

results = [];

if length(sub1) && isempty(sub2)
    error('both sub1 and sub2 must be provided!')
end
if length(sub2) && isempty(sub1)
    error('both sub1 and sub2 must be provided!')
end

if ischar(sub1)
    fprintf('reading sub1 from: %s\n', sub1)
    sub1 = importdata(sub1);
    sub1 = matlab.lang.makeValidName(sub1);
end

if ischar(sub2)
    fprintf('reading sub2 from: %s\n', sub2)
    sub2 = importdata(sub2);
    sub2 = matlab.lang.makeValidName(sub2);
end

usersub1 = []
usersub2 = []
if length(sub1)
    usersub1 = sub1;
    usersub2 = sub2;
end

if isempty(distmetric)
    distmetric = 'correlation';
end
results.Distance = distmetric;

if ~iscell(x1)
   x1 = {x1};
end

if ~iscell(x2)
   x2 = {x2};
end

% check if XX needs to be read in as well
readXX = '';
if ~isempty(conf) && isstr(conf)
    fprintf('conf provided as string %s: will use it as spatial confound (XX)\n', conf)
    readXX = conf;
    conf = [];
    conf.XX = readXX;
end

if ~isempty(conf) && isfield(conf, 'XX')
    if isstr(conf.XX)
        fprintf('wil use %s as spatial confounds per dataset\n', conf.XX)
        readXX = conf.XX;
        conf.XX = {};
    end
end

results.data1 = 'x1'; results.data2 = 'x2';
results.gotfiles = {};

for ii=1:length(x1)
  if ischar(x1{ii}) && exist(x1{ii})
    assert(~isempty(conn))
    results.gotfiles = [results.gotfiles, x1{ii}];
    results.data1 = x1;
    if isfile(x1{ii})
        x1{ii} = load(x1{ii});
    elseif isfolder(x1{ii})
        x1{ii} = read_data_conn_sub(x1{ii}, patt);
    else
        error(['not sure how to read conn from: ', x1{ii}])
    end
  end
  % approximate field name match
  if length(usersub1)
    x1{ii} = rename_fields(x1{ii}, usersub1);
  end
  if isempty(sub1)
      sub1 = fields(x1{ii});
  else
      sub1 = intersect(sub1, fields(x1{ii}));
  end
end
fprintf('%d subjects in first dataset\n', length(sub1))

for ii=1:length(x2)
  if ischar(x2{ii}) && exist(x2{ii})
    assert(~isempty(conn))
    results.gotfiles = [results.gotfiles, x2{ii}];
    results.data1 = x2;
    if isfile(x2{ii})
        x2{ii} = load(x2{ii});
    elseif isfolder(x2{ii})
        x2{ii} = read_data_conn_sub(x2{ii}, patt);
    else
        error(['not sure how to read conn from: ', x2{ii}])
    end
  end
  % approximate field name match
  if length(usersub2)
    x2{ii} = rename_fields(x2{ii}, usersub2);
  end
  if isempty(sub2)
      sub2 = fields(x2{ii});
  else
      sub2 = intersect(sub2, fields(x2{ii}));
  end
end
fprintf('%d subjects in second dataset\n', length(sub2))

common = NaN;
assert(~isempty(sub1))
assert(~isempty(sub2))
common = intersect(sub1, sub2);
if isempty(common)
    fprintf('no common subjects found: will try to find a maximum match\n')
    [common, sub1max, sub2max] = max_common_substr(sub1, sub2, '_');
    if ~isempty(common)
            fprintf('found %d common subjects: example id %s\n', length(common), common{1})
            fprintf('renaming the original struct\n')
            for ii=1:length(x1)
                x1{ii} = cell2struct(struct2cell(x1{ii}), sub1max);
            end
            for ii=1:length(x2)
                x2{ii} = cell2struct(struct2cell(x2{ii}), sub2max);
            end
    else
        error('could not figure out common subjects')
    end
end

assert(~isempty(common))
sub1 = common;
sub2 = common;


if subset > 0
    fprintf('taking a subset of size: %d\n', subset)
    assert(subset < length(sub1));
    subset = randsample(1:length(sub1), subset);
    sub1 = sub1(subset);
    sub2 = sub1;
end

% actually get the measures
x1XX = {};
for ii=1:length(x1)
    fprintf('getting connectivity: %s %s\n', conn, readXX)
    if ~isempty(readXX)
        x1XX{ii} = get_glocal_sub(x1{ii}, readXX, sub1);
    end
    x1{ii} = get_glocal_sub(x1{ii}, conn, sub1);
end

x2XX = {};
for ii=1:length(x2)
    fprintf('getting connectivity: %s %s\n', conn, readXX)
    if ~isempty(readXX)
        x2XX{ii} = get_glocal_sub(x2{ii}, readXX, sub2);
    end
    x2{ii} = get_glocal_sub(x2{ii}, conn, sub2);
end

% align gradients if we need to average them
% this wont work if there are any nan / inf in the matrix
if contains(conn, 'grad')
  if length(x1) > 1
    fprintf('aligning gradients within first datasets\n')
    x1 = gradient_align(x1);
  end

  if length(x2) > 1
    fprintf('aligning gradients within second datasets\n')
    x2 = gradient_align(x2);
  end % x2
end % gradient alignment

if contains(readXX, 'grad')
  if length(x1XX) > 1
    fprintf('aligning gradients within first datasets\n')
    x1XX = gradient_align(x1XX);
  end

  if length(x2XX) > 1
    fprintf('aligning gradients within second datasets\n')
    x2XX = gradient_align(x2XX);
  end % x2
end % gradient alignment


% average connectivity within each dataset, e.g. across sessions
if length(x1)>1
  fprintf('averaging first datasets\n')
  x1 = plus(x1{:})./length(x1);
else
  x1 = x1{1};
end

if length(x2)>1
  fprintf('averaging second datasets\n')
  x2 = plus(x2{:})./length(x2);
else
  x2 = x2{1};
end
% also average XX confounds
if ~isempty(readXX)
    if length(x1XX)>1
        fprintf('averaging first datasets XX\n')
        x1XX = plus(x1XX{:})./length(x1XX);
    else
        x1XX = x1XX{1};
    end

    if length(x2XX)>1
        fprintf('averaging second datasets XX\n')
        x2XX = plus(x2XX{:})./length(x2XX);
    else
        x2XX = x2XX{1};
    end
end


% from here on x1 and x2 should be matrices

if ~isempty(sub1) && ~isempty(sub2)
  assert(length(sub1)==size(x1,1))
  assert(length(sub2)==size(x2,1))
  common = intersect(sub1, sub2);
  [~,i] = ismember(common, sub1);
  x1 = x1(i,:);
  if ~isempty(readXX)
      x1XX = x1XX(i,:);
  end
  [~,i] = ismember(common, sub2);
  x2 = x2(i,:);
  if ~isempty(readXX)
    x2XX = x2XX(i,:);
  end
end
results.sub1 = sub1;
results.sub2 = sub2;
results.sub = common;

fprintf('%d subjects in common\n', length(common))

% replace NaNs with 0: make a user option?
if ~contains(conn, 'grad')
  x1(isnan(x1)) = 0;
  x2(isnan(x2)) = 0;
end

if ~isempty(readXX) && ~contains(readXX, 'grad')
  x1XX(isnan(x1XX)) = 0;
  x2XX(isnan(x2XX)) = 0;
end

% remove any rows with nan
nans = [find(any(isnan(x1),2)); find(any(isnan(x2),2))];
nans = unique(nans);
if ~isempty(nans)
   fprintf('removing %d sub with NaN\n', length(nans))
   x1(nans,:) = [];
   x2(nans,:) = [];
end

assert(all(size(x1)==size(x2)))

if ~isempty(feat)
    fprintf('using %d features\n', length(feat))
    x1 = x1(:,feat);
    x2 = x2(:,feat);
end

% confound removal, both feature-wise and observation-wise
if ~isempty(readXX)
    fprintf('will remove spatial confounds per subject\n')
    assert(all(size(x1)==size(x1XX)))
    assert(all(size(x2)==size(x2XX)))
    conf.XX = {x1XX, x2XX};
end
confremoved = remove_confounds({x1, x2}, conf, results.sub);
x1 = confremoved{1};
x2 = confremoved{2};
clearvars confremoved;

% get Idiff
% https://www.nature.com/articles/s41598-018-25089-1
% also keep self correlations
cr = corr(x1', x2');
results.Iself.Pearson = diag(cr);
cr = cr + diag(repmat(nan,1,size(cr,1))); % set diagonal to nan
results.Imargin.Pearson = results.Iself.Pearson - max([max(cr); nanmax(cr, [], 2)'])';
results.Iothers.Pearson = cr;
results.Idiff.Pearson = (nanmean(results.Iself.Pearson) - nanmean(results.Iothers.Pearson(:))) * 100;

cr = corr(x1', x2', 'type', 'spearman');
results.Iself.Spearman = diag(cr);
cr = cr + diag(repmat(nan,1,size(cr,1)));
results.Imargin.Spearman = results.Iself.Spearman - max([max(cr); nanmax(cr, [], 2)'])';
results.Iothers.Spearman = cr;
results.Idiff.Spearman = (nanmean(results.Iself.Spearman) - nanmean(results.Iothers.Spearman(:))) * 100;


K = max(max(ks),2); % use minimum 2 so we can get distance-margin
results.K = K;
results.ks = ks;

if contains(conn, 'grad') || contains(distmetric, 'grad')
    fprintf('using gradient_distance metric\n')
    distmetric = str2func('gradient_distance');
    [results.x2query.idx, results.x2query.dist] = knnsearch(x1,x2,'K',K,'Distance',distmetric);
    [results.x1query.idx, results.x1query.dist] = knnsearch(x2,x1,'K',K,'Distance',distmetric);
elseif strcmp(distmetric, 'weighted_correlation')
    fprintf('using weighted_correlation_distance metric\n')
    distmetric = str2func('weighted_correlation_distance');
    [results.x2query.idx, results.x2query.dist] = knnsearch(x1,x2,'K',K,'Distance',distmetric, 'w', feat);
    [results.x1query.idx, results.x1query.dist] = knnsearch(x2,x1,'K',K,'Distance',distmetric, 'w', feat);
else
    fprintf('using %s metric\n', distmetric)
    [results.x2query.idx, results.x2query.dist] = knnsearch(x1,x2,'K',K,'Distance',distmetric);
    [results.x1query.idx, results.x1query.dist] = knnsearch(x2,x1,'K',K,'Distance',distmetric);
end

results.x1query.acc = nan(length(ks),1);
results.x2query.acc = nan(length(ks),1);

results.x1query.ident = nan(length(ks),1);
results.x2query.ident = nan(length(ks),1);

for k=1:length(ks)

 % x2 is query
 acc = 0;
 ident = 0; % identification accuracy using mean of empirical
 identMAX = 0; % identification accuracy using max of empirical
 identMIN = 0;  % identification accuracy using min of empirical
 for i=1:size(x2,1)
    acc = acc + ismember(i, results.x2query.idx(i,1:ks(k)));
    if ks(k)==1 % calculate ident and identMAX
        dd = results.x2query.dist(:,1);
        ddi = dd(i);
        dd(i) = nan; % set self to nan to get loo empirical
        dd(1:length(dd) ~= results.x2query.idx(:,1)') = nan; % set wrongly identified to nan effectively only using
        % correctly identified subjects to create the empirical distr
        if ddi < nanmean(dd) && results.x2query.idx(i,1)==i
            ident = ident+1;
        end
        if ddi < nanmax(dd) && results.x2query.idx(i,1)==i
            identMAX = identMAX+1;
        end
        if ddi < nanmin(dd) && results.x2query.idx(i,1)==i
            identMIN = identMIN+1;
        end
    end
 end % for x2
 results.x2query.acc(k) = acc/size(x2,1);
 results.x2query.ident(k) = ident/size(x2,1);
 results.x2query.identMAX(k) = identMAX/size(x2,1);
 results.x2query.identMIN(k) = identMIN/size(x2,1);

 % x1 is query
 acc = 0;
 ident = 0;
 identMAX = 0;
 identMIN = 0;
 for i=1:size(x1,1)
    acc = acc + ismember(i, results.x1query.idx(i,1:ks(k)));
    if ks(k)==1
        dd = results.x1query.dist(:,1);
        ddi = dd(i);
        dd(i) = nan;
        dd(1:length(dd) ~= results.x1query.idx(:,1)') = nan;
        if ddi < nanmean(dd) && results.x2query.idx(i,1)==i
            ident = ident+1;
        end
        if ddi < nanmax(dd) && results.x2query.idx(i,1)==i
            identMAX = identMAX+1;
        end
        if ddi < nanmin(dd) && results.x2query.idx(i,1)==i
            identMIN = identMIN+1;
        end
    end
 end % for x1
 results.x1query.acc(k) = acc/size(x1,1);
 results.x1query.ident(k) = ident/size(x1,1);
 results.x1query.identMAX(k) = identMAX/size(x1,1);
 results.x1query.identMIN(k) = identMIN/size(x1,1);

 results.avg.acc(k) = (results.x1query.acc(k) + results.x2query.acc(k))/2;

 fprintf('k=%d accuracy: x1query=%g x2query=%g avg=%g\n', k, results.x1query.acc(k), results.x2query.acc(k), results.avg.acc(k))
end % for k

