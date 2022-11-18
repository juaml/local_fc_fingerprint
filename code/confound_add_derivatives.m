function [conf, Qreg, Nreg, nDerivs] = confound_add_derivatives(conf, Qreg, Nreg, ConfoundsTSderiv) 

if ~iscell(ConfoundsTSderiv)
   confds = split(ConfoundsTSderiv, '+');
end

assert(length(Qreg)==length(Nreg))

nconf = size(conf.reg,2);

for i=1:length(confds)
   for j=1:length(Nreg)
       if strcmpi(Nreg{j}, confds{i})
           conf.reg = [conf.reg, [0; diff(conf.reg(:,Qreg(j)))]];
           Nreg = [Nreg, [Nreg{j} 'deriv']];
           Qreg = [Qreg, size(conf.reg,2)];
       end
   end
end

nDerivs = size(conf.reg,2) - nconf;


