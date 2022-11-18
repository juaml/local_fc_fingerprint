function [NetFlags, des] = analysis_settings(analysis, NetFlags)

GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
PCA  = {'PCA','NoPCA','FIX'};
Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
FILT = {'BP','HP'};

if ~isempty(analysis)
  try
      in = []; for xi=1:numel(GM); in(xi) = numel(strfind(analysis,['_' GM{xi}])); end
      NetFlags.gotMask = find(in);

      in = []; for xi=1:numel(PCA); in(xi) = numel(strfind(analysis,['_' PCA{xi}])); end
      NetFlags.gotPCA = find(in);

      in = []; for xi=1:numel(Glob); in(xi) = numel(strfind(analysis,['_' Glob{xi}])); end
      NetFlags.gotGlob = find(in);

      in = []; for xi=1:numel(FILT); in(xi) = numel(strfind(analysis,['_' FILT{xi}])); end
      NetFlags.gotFilt = find(in);
  catch
      error(['could not parse RS-FC settings', analysis])
  end
end

des = [GM{NetFlags.gotMask} '_' PCA{NetFlags.gotPCA} '_' Glob{NetFlags.gotGlob} '_' FILT{NetFlags.gotFilt}];


