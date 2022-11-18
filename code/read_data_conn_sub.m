function conn = read_data_conn_sub(folder, patt)
switch nargin
    case 0
        error('give me a folder')
    case 1
        patt = '*conn*';
end

if ~isfolder(folder)
    error(['not a folder: ' folder])
end

if isempty(patt)
    patt = '*';
end

files = dir(fullfile(folder, [patt '.mat']));
conn = struct();
for f=1:length(files)
    fl = fullfile(files(f).folder, files(f).name);
    x = load(fl);
    sub = fields(x);
    assert(length(sub)==1)
    sub = sub{1};
    if ismember(sub, fields(conn))
        error(['multiple subject files please set the pattern properly: ', sub])
    end
    conn.(sub) = x.(sub);
end

