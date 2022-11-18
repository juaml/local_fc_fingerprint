function settings = kp_get_settings(NetFlags)

settings = '';

% override Filt if provided directly with NetFlags
if isfield(NetFlags, 'filterFrequency') % this can be empty
    Filt = NetFlags.filterFrequency;
    assert(isempty(Filt) || length(Filt)==2)
else
    error('must have NetFlags.filterFrequency!!!')
end

if length(Filt)==2
    assert(Filt(1)<Filt(2))
    settings = [settings, 'filt', num2str(Filt)];
else
    settings = [settings, 'filtNO'];
end

regress_confounds = true;
if isfield(NetFlags, 'regress_confounds') && ~isempty(NetFlags.regress_confounds)
    regress_confounds = NetFlags.regress_confounds;
    if isfield(NetFlags, 'ConfoundsTS')
        settings = [settings, '_conf', NetFlags.ConfoundsTS];
    end

    if isfield(NetFlags, 'ConfoundsTSderiv')
        settings = [settings, '_confdrv', NetFlags.ConfoundsTSderiv];
    end
end

badTP = false;
if isfield(NetFlags, 'confound_badTP')
    badTP = NetFlags.confound_badTP;
end
if badTP
    settings = [settings, '_badTP'];
end

if isfield(NetFlags, 'global_mask_to_use')
    settings = [settings, '_mask', NetFlags.global_mask_to_use];
end

if isfield(NetFlags, 'subject_gm_mask') && NetFlags.subject_gm_mask>0
    settings = [settings, '_submask', num2str(NetFlags.subject_gm_mask)];
end

detrendTS = false;
if isfield(NetFlags, 'detrendTS') && ~isempty(NetFlags.detrendTS)
    detrendTS = NetFlags.detrendTS;
end
if detrendTS>0
    settings = [settings, '_dt' num2str(detrendTS)];
end

detrendConfTS = false;
if isfield(NetFlags, 'detrendConfTS') && ~isempty(NetFlags.detrendConfTS)
    detrendConfTS = NetFlags.detrendConfTS;
end
if detrendConfTS>0
    settings = [settings, '_dtconf' num2str(detrendConfTS)];
end

settings = matlab.lang.makeValidName(settings);

