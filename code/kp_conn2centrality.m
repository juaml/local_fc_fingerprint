function [cent] = kp_conn2centrality(conn, dens, measures, takeabs)

switch nargin
    case 0
        error('please give me some arguments!')
    case 1
        dens = 0.1; measures = []; takeabs=false;
    case 2
        measures = []; takeabs=false;
    case 3
        takeabs=false;
end

if exist('clustering_coef_wu')~=2
    error('BCT toolbox is not in the path')
end

conn(isnan(conn)) = 0;

if isempty(dens) || isnan(dens)
    dens = 0.1;
end
assert(dens > 0)

if isempty(measures)
    measures = {'degree','closeness','betweenness','pagerank','eigenvector'}; %,'topology'};
end

cent = struct();
cent.dens = dens;
cent.measures = measures;
cent.topology = [];
for ith_m=1:length(measures)
    cent.(measures{ith_m}) = [];
end

if takeabs
    conn = abs(conn);
end

    if ismatrix(conn)
        confull = conn;
    else
        connfull = tril2square(conn);
    end
    conn = threshold_proportional(connfull, dens);
    conn = conn - diag(diag(conn)); % make diagonal zero
    gp   = graph(conn, 'lower');
    imp  = gp.Edges.Weight;
    cst  = 1./imp;
    for ith_m=1:length(measures)
        switch measures{ith_m}
            case {'degree','pagerank','eigenvector'}
        	    cn = centrality(gp, measures{ith_m},'Importance',imp)';
	        case {'closeness','betweenness'}
		        cn = centrality(gp, measures{ith_m},'Cost',cst)';
	        case 'topology'
    		    % get topology measurements
    		    % get measurements on random networks for normalization
   		        ccrand = nan(1,20);
    		    gerand = nan(1,20);
    		    rcrand = nan(1,20);
    		    for ii=1:20
        		    %using following two lines gives same values as original matrix
        		    %I think this is what this paper says that they do but perhaps not
        		    %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837693/
        		    %aa = randmio_und_connected(connfull, 3);
        		    %aa = threshold_proportional(aa, dens);
        		    %directly rewiring the thresholded matrix seems to provide reasonable numbers
        		    aa = randmio_und_connected(conn, 3);
        		    ccrand(ii) = nanmean(clustering_coef_wu(aa));
        		    gerand(ii) = efficiency_wei(aa);
        		    rcrand(ii) = nanmean(rich_club_wu(aa));
    		    end

    		    topo = [];
    		    ccnorm = nanmean(clustering_coef_wu(conn))/nanmean(ccrand);
    		    genorm = efficiency_wei(conn)/nanmean(gerand);
    		    topo = [topo, ccnorm];
    		    topo = [topo, genorm];
    		    topo = [topo, nanmean(rich_club_wu(conn))/nanmean(rcrand)];
    		    topo = [topo, ccnorm*genorm]; % this small world index
    		    [~,Q] = modularity_und(conn);
    		    topo = [topo, Q];
    		    topo = [topo, transitivity_wu(conn)];
    		    D = 1./conn;
    		    D(conn==0) = Inf;
    		    D(isnan(D)) = 0;
    		    D = D - diag(diag(D));
    		    D = distance_wei(D);
    		    [lambda, eff] = charpath(D,0,0); % the doc says argument must be distance matrix
    		    topo = [topo, lambda, eff];
    		    topo(isinf(topo)) = 0;
		        cn = topo;
	        otherwise
		        error('unknown measure')
	        end
        cent.(measures{ith_m}) = cn;
    end % for length(measures)

