% FC density as described in: Tomasi, Volkow, Functional connectivity density mapping, PNAS 2010
% as well as local and distant connectivity as described in Sepulcre et al., % PLoS Comp Biol 2010
% get FCS: here we average over the ROI where the original paper uses voxel-wise definition
% Coupling of functional connectivity and regional cerebral blood flow reveals a physiological basis
% for network hubs of the human brain.
% Liang X, Zou Q, He Y, Yang Y, PNAS 2013:
% Coupling of functional connectivity and regional cerebral blood flow reveals a physiological basis for network hubs of the human brain
% this is also similar to wGBC
% Cole et al., NeuroImage 2010: Identifying the brain's most globally connected regions


function fcden = kp_FCdensity_corr(z, wbXYZmm, roiXYZmm)

fcden = [];

fcden.fcd06mean = nanmean(nansum(z>0.6, 2));
fcden.fcd06max = nanmax(nansum(z>0.6, 2));

fcden.fcd04mean = nanmean(nansum(z>0.4, 2));
fcden.fcd04max = nanmax(nansum(z>0.4, 2));

fcden.fcs    = nanmean(nanmean(z,2)); % first average voxel-WB and then across voxels
fcden.fcsabs = nanmean(nanmean(abs(z),2));
fcden.fcspos = nanmean(z(z>0));
fcden.fcsneg = nanmean(z(z<0));

if ~isempty(wbXYZmm) && ~isempty(roiXYZmm)
    D = pdist2(roiXYZmm', wbXYZmm','euclidean','Smallest',1);
    assert(length(D)==size(z,2))
    DD = D;
    DD(D > 14 | D==0) = NaN; % local
    fcden.local14mm025mean = nanmean(nansum(z.*DD>0.25, 2));
    fcden.local14mm025max = nanmax(nansum(z.*DD>0.25, 2));
    DD = D;
    DD(D<=14) = NaN; % distance
    fcden.distant14mm025mean = nanmean(nansum(z.*DD>0.25, 2));
    fcden.distant14mm025max = nanmax(nansum(z.*DD>0.25, 2));
end
