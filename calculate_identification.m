%Definitions
parcellation=["Power_3mm_ds000115_fmriprep"];
meas=["ReHo", "connm", "falff","alff"];
zeroback={'ds000115_0back'};
twoback={'ds000115_2back'};
datasets=[zeroback,twoback];

currfolder=pwd;
cd('code')

% Identification
for k=1:length(meas)
    file1=fullfile(currfolder,'output',datasets{1},parcellation{1});
    file2=fullfile(currfolder,'output',datasets{2},parcellation{1});
    I = kp_identification_run(file1,file2,meas{k},[]);
    id=table(I.x1query.acc,I.x2query.acc,I.Idiff.Spearman,...
        'variablenames',{'identification_accuracy_file1_file2',...
        'identification_acccuracy_file2_file1',...
        'Idiff_spearman'});
    filename = fullfile(currfolder,'output',['identification_' meas{k} '.csv']);
    writetable(id,filename)
end
