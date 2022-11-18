function x = remove_confounds(x, conf, sub, checklowvar)

switch nargin
    case {0,1}
        error('I need at least two arguments!')
    case 2
        sub = []; checklowvar = true;
    case 3
        checklowvar = true;
end


if ~iscell(x)
	x = {x};
end

if ~isempty(sub) && isfield(conf, 'sub') && ~isempty(conf.sub)
	assert(iscell(conf.sub))
    confsub = {};
    for ii=1:length(conf.sub)
        assert(length(conf.sub{ii})==1)
        confsub{ii} = conf.sub{ii}{1};
        if isnumeric(confsub{ii})
            confsub{ii} = num2str(confsub{ii});
        end
    end
    conf.sub = confsub;
    fprintf('example subject id in conf.sub: %s\n', conf.sub{1})
    conf.sub = matlab.lang.makeValidName(conf.sub);
	fprintf('reordering confounds according to sub\n')
	[ii,jj] = ismember(sub, conf.sub);
	if any(ii==0)
        fprintf('some of the subjects not found in confounds\n')
        % see if we can find a match
        fprintf('checking if maximum substring can be matched\n')
        [common, sub, conf.sub] = max_common_substr(sub, conf.sub, '_', true);
        if length(common)==length(sub)
		    fprintf('found perfect match\n')
            [ii,jj] = ismember(sub, conf.sub);
        else
		    notfound = find(ii==0);
		    notfound = notfound(1:min(5, length(notfound)));
		    sub(notfound)
        end
	end
	if isfield(conf, 'X') && ~isempty(conf.X)
		conf.X = conf.X(jj,:);
		conf.X = [conf.X, ones(size(conf.X,1),1)];
	end
	if isfield(conf, 'XX') && ~isempty(conf.XX)
        if iscell(conf.XX)
            assert(length(conf.XX) == length(x))
            for iXX=1:length(conf.XX)
                conf.XX{iXX} = conf.XX{iXX}(jj,:);
            end
        else
		    conf.XX = conf.XX(jj,:);
        end
	end
end

if ~isfield(conf, 'data') || isempty(conf.data)
	conf.data = 1:length(x);
end

if isfield(conf, 'X') && ~isempty(conf.X)
	% now remove the confounds from each connection/feature
	for i=1:length(x)
		if ismember(i, conf.data)
			for j=1:size(x{i},2)
                doit = true;
                % low variance can lead to data leakage so avoid it
                if checklowvar && var(x{i}(:,j))<eps
                    doit = false;
                end

                if doit % we have already ones added above
				    [~,~,x{i}(:,j)] = regress(x{i}(:,j), conf.X);
                end
			end
		end
	end
end

if isfield(conf, 'XX') && ~isempty(conf.XX)
	% now remove the confounds from each subject/observation
	for i=1:length(x)
		if ismember(i, conf.data)
			for j=1:size(x{i},1)
                doit = true;
                if checklowvar && var(x{i}(j,:)')<eps
                    doit = false;
                end

                if doit
                    if iscell(conf.XX)
                        XX = conf.XX{i}(j,:);
                    else
                        XX = conf.XX(j,:);
                    end
				    XX = [ones(length(XX),1), XX'];
				    [~,~,XX] = regress(x{i}(j,:)', XX);
				    x{i}(j,:) = XX';
                end
			end
		end
	end
end

