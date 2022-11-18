function Qreg = confound_indices(gotPCA, gotGlob, datadir)
% regressors
%1. meanGM
%2. meanWM
%3. meanCSF
%4. global signal
%5. meanGM^2
%6. meanWM^2
%7. meanCSF^2
%8. global signal^2
%9.-32. movement parameter

% older and less flexible function which expects a fixed file structure

% WARN: using all confounds for HCP irrespective of what was asked!
if ~isempty(strfind(upper(datadir), 'HCP'))
   Qreg = 1:22;
   return
end

isfixaroma = 0;
if ~isempty(strfind(upper(datadir), 'FIX'))
  isfixaroma = 1;
elseif ~isempty(strfind(upper(datadir), 'AROMA'))
  isfixaroma = 2;
end

Qreg = [];
if gotPCA==1
    Qreg = [Qreg 33:37];
end

if gotGlob==1
    Qreg = [Qreg 1 2 3 5 6 7];
elseif gotGlob==2
    Qreg = [Qreg 1:3];
elseif gotGlob==3
    Qreg = [Qreg 2 3 6 7];
elseif gotGlob==4
    Qreg = [Qreg 2:3];
elseif gotGlob==5
    Qreg = [Qreg 4 8];
elseif gotGlob==6
    Qreg = [Qreg 4];
end

if gotPCA==3 || isfixaroma
   if isfixaroma==1 || isfixaroma==2 % fix or aroma
      Qreg = [Qreg 2:3];
   end
else
   Qreg = [9:32 Qreg];
end

Qreg = unique(Qreg);

