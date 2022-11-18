% to be used for retults with connectDim_sub
function [cn, who] = get_glocal_sub(x, what, who, whr, perc, rmnvox)
  
  switch nargin
    case {0,1}
      error('give me arguments')
    case 2
      who = fieldnames(x); whr = 'connectDim'; perc = 10; rmnvox=false;
    case 3
      whr = 'connectDim'; perc = 10; rmnvox = false;
    case 4
      perc = 10; rmnvox = false;
    case 5
      rmnvox = false;
    end
    
    if isempty(who)
      who = fieldnames(x);
    end
    who = matlab.lang.makeValidName(who);
    
    if ~all(ismember(who, fields(x)))
      fprintf('some of the subjects that you asked for are not available\n')
      ii = find(~ismember(who, fields(x)));
      fprintf('for example: %s is not available\n', who{ii(1)})
      error('subject(s) not found')
    end
    
    what2 = what;
    left = false;
    right = false;
    add = false;
    wbrand = false;
    roi = false;
    varexpl = 0;
    dists = '';
    if ~isempty(regexp(what,'^WB'))
      what2 = 'WB';
      
      dists =  regexp(what, 'WB(\d+)mm', 'tokens');
      if ~isempty(dists)
        dists = [dists{1}{1} 'mm'];
        what = strrep(what, dists, '');
      end
      % check for '_' which means voxels within this distance were removed
      if isempty(dists)
        dists =  regexp(what, 'WB_(\d+)mm', 'tokens');
        if ~isempty(dists)
          dists = ['_' dists{1}{1} 'mm'];
          what = strrep(what, dists, '');
        end
      end
      
      if isempty(dists)
        dists = '';
      else
        dists
      end
      
      if ~isempty(regexp(what,'Left'))
        left = true;
        what = strrep(what,'Left','');
      elseif ~isempty(regexp(what,'Right'))	
        right = true;
        what = strrep(what,'Right','');
      end
      
      if ~isempty(regexp(what,'WBlaplacegeom'))
        what2 = 'WBlaplacegeom';
      elseif ~isempty(regexp(what,'WBlaplacevar'))
        what2 = 'WBlaplacevar';
      elseif ~isempty(regexp(what,'WBlaplace'))
        what2 = 'WBlaplace';
      elseif ~isempty(regexp(what,'WBgeom'))
        what2 = 'WBgeom';
      elseif ~isempty(regexp(what,'AUC'))
        what2 = 'WBAUC';
      elseif ~isempty(regexp(what,'Kaiser'))
        what2 = what;
      else
        varexpl = regexp(what,'\d+$','match');
        assert(length(varexpl)==1)
        varexpl = str2num(varexpl{1})/100;
      end
      
      if ~isempty(regexp(what,'Sum'))
        what = strrep(what,'Sum','');
        add = true;
      end
      
      if ~isempty(regexp(what,'rand'))
        what = strrep(what,'rand','');
        wbrand = true;
      end
      
      if ~isempty(regexp(what,'ROI'))
        what = strrep(what,'ROI','');
        roi = true;
      end
    end
   
    reverseStr = '';
    cn = [];
    for i=1:length(who)
      %fprintf('.')
      xx = [];  % this is actual measurement
      notfound = false;
      %if ismember(lower(what),{'hacf','connm','conne','connmean', 'grad','connmgrad', 'gradna', 'gradpca', 'connev','connmlaplace', 'connelaplace','topo', 'topology','degree','between','page','eigen','close','connm90','conne90','connm95','conne95','connm99','conne99'})
      switch lower(what)
        case {'edge', 'edgemean'}
          xx = x.(who{i}).edge.mean;
        case {'hacf0.05', 'hacf0.1', 'hacf0.25', 'hacf0.5'}
          frackeep = matlab.lang.makeValidName(strrep(lower(what), 'hacf', ''));
          xx = x.(who{i}).roipair.hacfmean.(frackeep);
        case {'connm', 'connmean'}
          xx = x.(who{i}).roipair.mean;
        case {'connmpartial', 'connmeanpartial'}
          % take correlation and convert it to partial
          xx = x.(who{i}).roipair.mean;
          xx(isnan(xx)) = eps;
          xx(isinf(xx)) = eps;
          xx = cor2pcor(tril2square(xx, 1.0)); % 1.0 is the diagonal
          xx = xx(tril(true(size(xx)),-1))';
        case {'grad','connmgrad','gradpca','gradna'}
          xx = x.(who{i}).roipair.mean;
          xx = tril2square(xx, 1.0); % get the symmetric full connectome
          xx(isnan(xx)) = 0;
          switch lower(what)
            case {'connmgrad', 'grad'}
              holder = evalc('gm = GradientMaps(''n_components'', 1)');
            case {'gradpca'}
              holder = evalc('gm = GradientMaps(''n_components'', 1, ''approach'', ''pca'')');
            case 'gradna'
              holder = evalc('gm = GradientMaps(''n_components'', 1, ''kernel'', ''na'')');
            otherwise
              error(['unknown gradient: ' what])
          end
          try
            if contains(lower(what), 'sparse0')
              holder = evalc('gm = gm.fit(xx, ''sparsity'', 0)');
            else
              holder = evalc('gm = gm.fit(xx)'); % evalc to avoid printed output
            end
            assert(length(gm.gradients{1}(:,1)) == size(xx,1))
            xx = gm.gradients{1}(:,1)';
          catch ME
            lowvar = var(xx) < eps;
            if sum(lowvar)
              xx = xx(~lowvar, ~lowvar);
              holder = evalc('gm = gm.fit(xx)');
              gr = nan(size(xx,1), 1);
              gr(~lowvar) = gm.gradients{1}(:,1);
              xx = gr';
            end
            %getReport(ME)
          end
        case {'conne', 'connev'}
          xx = x.(who{i}).roipair.ev;
        case {'connmlaplace', 'connm90', 'connm95', 'connm99'}
          xx = x.(who{i}).roipair.mean;
          xx = tril2square(xx, 1.0);
          xx(isnan(xx)) = 0;
          xx(isinf(xx)) = 0;
          eg = kp_eigval_mean0(xx);
          if strcmp(lower(what),'connmlaplace')
            xx = laplace_pca([], eg, size(xx,1), size(xx,2));
          elseif strcmp(lower(what),'connm90')
            xx = min(find(cumsum(eg)/sum(eg)>=0.90));
          elseif strcmp(lower(what),'connm95')
            xx = min(find(cumsum(eg)/sum(eg)>=0.95));
          elseif strcmp(lower(what),'connm99')
            xx = min(find(cumsum(eg)/sum(eg)>=0.99));
          end
        case {'connelaplace', 'conne90', 'conne95', 'conne99'}
          xx = x.(who{i}).roipair.ev;
          xx = tril2square(xx, 1.0);
          xx(isnan(xx)) = 0;
          xx(isinf(xx)) = 0;
          eg = kp_eigval_mean0(xx);
          if strcmp(lower(what),'connelaplace')
            xx = laplace_pca([], eg, size(xx,1), size(xx,2));
          elseif strcmp(lower(what),'conne90')
            xx = min(find(cumsum(eg)/sum(eg)>=0.90));
          elseif strcmp(lower(what),'conne95')
            xx = min(find(cumsum(eg)/sum(eg)>=0.95));
          elseif strcmp(lower(what),'conne99')
            xx = min(find(cumsum(eg)/sum(eg)>=0.99));
          end
        case {'topo', 'topology'}
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).topology;
        case 'degree'
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).degree;
        case 'between'
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).betweenness;
        case 'page'
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).pagerank;
        case 'eigen'
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).eigenvector;
        case 'close'
          xx = x.(who{i}).roipair.(['mean_centrality' num2str(perc) 'percent']).closeness;
        otherwise
          notfound = true;
      end % switch
      %end %if ismember
      
      if notfound
        gl = x.(who{i}).glocal;
        xx = nan(1, length(gl));
        for j=1:length(gl)
          try
            switch what2
              case 'kaiserwithin'
                eg = gl{j}.eigvalWITHIN;
                xx(j) = max(find(eg>=1));
              case 'laplacewithin'
                eg = gl{j}.eigvalWITHIN;
                xx(j) = laplace_pca([], eg, numel(eg), numel(eg)); 
              case 'eig1varwithin'
                xx(j) = gl{j}.eigvalWITHIN(1)/sum(gl{j}.eigvalWITHIN);
              case 'eig1'
                xx(j) = gl{j}.eigvalWB(1);
              case 'fcswb'
                xx(j) = gl{j}.FCdensity.fcswb;
              case 'fcswbabs'
                xx(j) = gl{j}.FCdensity.fcswbabs;
              case 'fcd06mean'
                xx(j) = gl{j}.FCdensity.fcd06mean;
              case 'fcd04mean'
                xx(j) = gl{j}.FCdensity.fcd04mean;
              case 'local14mm025mean'
                xx(j) = gl{j}.FCdensity.local14mm025mean;
              case 'distant14mm025mean'
                xx(j) = gl{j}.FCdensity.distant14mm025mean;
              case {'ReHo','reho'}
                xx(j) = gl{j}.reho;
              case {'ReHots', 'rehots'}
                xx(j) = gl{j}.rehots;
              case {'tsnr','tSNR'}
                xx(j) = gl{j}.tsnr;
              case {'sampen','SampEn'}
                xx(j) = gl{j}.sampen;
              case {'ReHo10','reho10'}
                xx(j) = gl{j}.reho10;
              case {'ALFF', 'alff'}
                xx(j) = gl{j}.alff;
              case {'FALFF', 'falff'}
                xx(j) = gl{j}.falff;
              case {'nvox', 'nvoxWB'}
                switch what2
                  case 'nvox'
                    xx(j) = gl{j}.nvox;
                  case 'nvoxWB'
                    xx(j) = gl{j}.nvoxWB;
                  end
                case {'WBAUC','WBlaplace','WB','WBlaplacegeom','WBlaplacevar','WBgeom','WBKaiser','WBKaiservar','WBKaisergeom'}
                  if wbrand
                    eigval = gl{j}.eigvalWB_rand;
                    nvoxWB = gl{j}.nvoxWB;
                  elseif left
                    eigval = gl{j}.eigvalWBLeft;
                    nvoxWB = gl{j}.nvoxWBLeft;
                  elseif right
                    eigval = gl{j}.eigvalWBRight;
                    nvoxWB = gl{j}.nvoxWBRight;
                  elseif roi
                    eigval = gl{j}.eigvalROIAV;
                    nvoxWB = length(gl); % this is #ROI
                  else
                    eigvalstr = ['eigvalWB' dists];
                    nvoxstr   = ['nvoxWB' dists];
                    eigval = gl{j}.(eigvalstr);
                    nvoxWB = gl{j}.(nvoxstr);
                  end
                  
                  if isempty(eigval) || all(isnan(eigval))
                    error('invalid eigval')
                  end
                  
                  nvox   = gl{j}.nvox;
                  evar   = cumsum(eigval)/sum(eigval);
                  switch what2
                    case 'WBKaiser'
                      xx(j) = max(find(eigval>=1));
                    case 'WBKaiservar'
                      ll = max(find(eigval>=1));
                      xx(j) = evar(ll);
                    case 'WBKaisergeom'
                      ll = max(find(eigval>=1));
                      xx(j) = geomean(eigval(1:ll));
                    case 'WBAUC'
                      xx(j) = trapz(evar);
                    case 'WBlaplace'
                      xx(j) = laplace_pca([], eigval, nvoxWB, nvox);
                    case 'WBlaplacegeom'
                      ll = laplace_pca([], eigval, nvoxWB, nvox);
                      xx(j) = geomean(eigval(1:ll));
                    case 'WBlaplacevar'
                      ll = laplace_pca([], eigval, nvoxWB, nvox);
                      xx(j) = evar(ll);
                    case 'WBgeom'
                      xx(j) = geomean(eigval);
                    case 'WB'
                      if add
                        jj = min(find(evar>=varexpl));
                        xx(j) = sum(eigval(1:jj));
                      else
                        xx(j) = min(find(evar>=varexpl));
                      end
                    end % switch what2
                  otherwise
                    xx = x.(who{i}).(whr).(what);
                  end
                  
                catch ME
                  fprintf('%s: %d %s\n', who{i}, j, ME.identifier)
                  xx(j) = NaN;
                end % try
              end % for j
              
              if size(xx,2)==1
                xx = xx';
              end
              
            end
            
      % add the measurement to the collection
      cn = [cn; xx];

	  % Display the progress
      reverseStr = print_progress(i, length(who), reverseStr);
    end % i=1:length(who)
    assert(size(cn,1) == length(who))

% remove number of voxel signal from each measurement  
if rmnvox
  nvox = get_glocal_sub(x, what, who, whr, perc, false);
  for i=1:size(cn,1)
    m = fitlm(nvox(i,:)', cn(i,:)');
    cn(i,:) = lm.Residuals.Raw;
  end
end
          
