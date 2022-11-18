function [Qreg, Nreg] = confound_indices_TS(confTS, regnames)
% confTS: a string describing which confounds to use
%         parts should be joined using '+'
%         e.g. WM+CSF+MOV to add WM mean, CSF mean and all movement regressors
%         to add squared values use 2, e.g. WM+WM2 will add WM and its squared value
%         to add PCA do the same as for square, e.g. WM+WMPCA

if ~iscell(confTS)
   confTS = split(confTS, '+');
end

% expand 
expanded = {};
for i=1:length(confTS)
    if ~isempty(regexpi(confTS{i}, '^WM'))
        xx = regexprep(confTS{i}, '^WM', 'white_matter', 'ignorecase');
        expanded = [expanded; xx];
    elseif ~isempty(regexpi(confTS{i}, '^GS'))
        xx = regexprep(confTS{i}, '^GS', 'global_signal', 'ignorecase');
        expanded = [expanded; xx];
    end
end
confTS = [confTS; expanded]

Nreg = {};
Qreg = [];
for i=1:length(confTS)
    for j=1:length(regnames)
      m = strcmpi(regnames{j},confTS{i}); % keep this case insensitive
      if ~m % no exact match so try .[/d]
          patt = ['^' regexptranslate('escape', confTS{i}) '\.(\d*)'];
          m = regexpi(regnames{j}, patt);
      end
      if m
         Qreg = [Qreg, j];
         Nreg = [Nreg, regnames{j}];
      end
    end
end

