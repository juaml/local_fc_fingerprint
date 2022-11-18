function [results, x1, x2] = kp_identification_run(folder1, folder2, conn, conf, sub, subset, patt, feat)
% folder1: can be string or cell array with strings indicating folder names with mat files
% folder2: same as above
% conn: connectivity measure to use
% conf: spatial confound (subject-wise) to remove from conn
% sub: which subjects to use
switch nargin
    case {0,1,2}
        error('need 3 arguments!')
    case 3
        conf = []; sub=''; subset=0; patt = ''; feat = [];
    case 4
        sub=''; subset=0; patt = ''; feat = [];
    case 5
        subset=0; patt = ''; feat = [];
    case 6
        patt=''; feat = [];
    case 7
        feat = [];
end

tStart = tic;

if ischar(folder1)
    folder1 = {folder1};
end

if ischar(folder2)
    folder2 = {folder2};
end

d1 = {};
for f=1:length(folder1)
    dd = read_data_conn_sub(folder1{f}, patt);
    if f==1
        s1 = fieldnames(dd);
    else
        s1 = intersect(s1, fieldnames(dd));
    end
    d1 = [d1, dd];
end

d2 = {};
for f=1:length(folder2)
    dd = read_data_conn_sub(folder2{f}, patt);
    if f==1
        s2 = fieldnames(dd);
    else
        s2 = intersect(s2, fieldnames(dd));
    end
    d2 = [d2, dd];
end


%d1 = read_data_conn_sub(folder1, patt);
%d2 = read_data_conn_sub(folder2, patt);

%[~, s1] = get_glocal_sub(d1, conn);
%[~, s2] = get_glocal_sub(d2, conn);

if isempty(sub)
    s = intersect(s1, s2);
else
    s = sub;
end
[results, x1, x2] = kp_identification(d1, d2, 1, conn, 'spearman', patt, conf, s, s, subset, feat);

results.timetaken = toc(tStart);
fprintf('time taken: %g\n', results.timetaken)


