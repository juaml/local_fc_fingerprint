function run_conn(settings_file, subs)

switch nargin
    case {0}
        error('give me settings_file')
    case {1}
        subs = [];
end

NetFlags = settings_conn_default(settings_file);

% setup
if isfield(NetFlags, 'settings') && ~isempty(NetFlags.settings)
    sett = str2func(NetFlags.settings);
    sett = sett();
    fls = fields(sett);
    for i=1:length(fls)
        fprintf('setting %s\n', fls{i})
        NetFlags.(fls{i}) = sett.(fls{i});
    end
end

if isfield(NetFlags, 'parpool') && NetFlags.parpool
    p = gcp('nocreate');
    if isempty(p)
        fprintf('creating a new parpool: %d\n', NetFlags.parpool)
        parpool(NetFlags.parpool)
    else
        fprintf('using existing parpool: %d\n', p.NumWorkers)
        NetFlags.parpool = p.NumWorkers;
    end
end

if isempty(NetFlags.connectome_cov)
    NetFlags.connectome_cov = 'COV';
end

if isempty(NetFlags.lookup_file)
    NetFlags.lookup_file = fullfile(pwd,'Lookups',[NetFlags.Group '.mat']);
end
fprintf('using lookup file: %s\n', NetFlags.lookup_file)

[~,~,lookext] = fileparts(NetFlags.lookup_file);
switch lookext
    case '.mat'
        look = load(NetFlags.lookup_file);
    case {'.txt', '.csv', '.tab'}
        look = csv2lookup(NetFlags.lookup_file);
    otherwise
          error('unknown lookup file extension: ', lookext)
end

if ~isempty(subs)
    fprintf('calculating on %d subjects\n', length(subs))
    assert(max(subs)<=size(look.Cov,1))
    look = lookup_filter(look, subs);
    if length(subs)<=5
        disp('subjects:')
        look.sub
    end
end

datafield = NetFlags.subdatafield;
incl = find(look.(datafield)==1);
assert(~isempty(incl))
fprintf('using %d subjects\n', length(incl))
look = lookup_filter(look, incl);
NetFlags.sub = look.sub;
NetFlags.Cov = look.Cov;
NetFlags.isPat = look.isPat;
NetFlags.SubDir = look.SubDir;
NetFlags.Covariates = look.Covariates;
NetFlags.Covariates = regexprep(NetFlags.Covariates, '\s+', ''); % replace any spaces

if ~isempty(NetFlags.connectome_nii)
    fprintf('connectivity using %d files\n', length(NetFlags.connectome_nii))
    kp_conn_glocal_sub_noparfor(NetFlags.connectome_nii,'',NetFlags);
end


