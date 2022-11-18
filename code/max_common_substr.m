function [common, sub1max, sub2max] = max_common_substr(sub1, sub2, splitstr, sub1fixed)

    switch nargin
        case {0,1}
            error('I need 2 arguments')
        case 2
            splitstr = '_'; sub1fixed = false;
        case 3
            sub1fixed = false;
    end
    common = [];
    sub1max = [];
    sub2max = [];

    if ~sub1fixed
        subref1 = strsplit(sub1{1}, splitstr);
        subref2 = strsplit(sub2{1}, splitstr);
        if length(subref1)==1 || length(subref2)==1 || length(subref1)~=length(subref2)
            fprintf('could not split using %s: cant figure out common strings', splitstr)
            return
        end
        sub1split = regexp(sub1, splitstr, 'split');
    else
        subref2 = strsplit(sub2{1}, splitstr);
        if length(subref2)==1
            fprintf('could not split using %s: cant figure out common strings', splitstr)
            return
        end
    end
    sub2split = regexp(sub2, splitstr, 'split');

    for irm=(length(subref2)-1):-1:1
        sub1max = cell(length(sub1),1);
        sub2max = cell(length(sub2),1);
        for ii=1:length(sub1)
            if sub1fixed
                sub1max{ii} = sub1{ii};
            else
                sub1max{ii} = [sub1split{ii}{1:irm}];
            end
        end
        for ii=1:length(sub2)
            sub2max{ii} = [sub2split{ii}{1:irm}];
        end
        if length(sub1max)~=length(unique(sub1max))
            break
        end
        if length(sub2max)~=length(unique(sub2max))
            break
        end
        common = intersect(sub1max, sub2max);
        if ~isempty(common)
            break
        end
    end

