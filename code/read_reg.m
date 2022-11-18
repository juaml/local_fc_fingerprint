function regextra = read_reg(fextra, basedir)

switch nargin
    case 0
        fextra = []; basedir = '';
    case 1
        basedir = '';
end

regextra = [];
for f=1:length(fextra)
    ff = fextra{f};
    if ~isempty(basedir)
        ff = fullfile(basedir, ff);
    end
    regextra = [regextra, textread(ff)];
end
regextra = regextra(:, var(regextra)>eps);

