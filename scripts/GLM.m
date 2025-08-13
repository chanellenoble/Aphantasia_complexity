%%OUTPUTS-
% All_glm_fc

%% GLM of timeseries to aphantasia task

cd '/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/'

%need to load 'All_onsets', 'All_dsmtx', & 'subject_id'
load('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/All_onsets.mat')

load('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/All_dsmtx.mat')

load('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/subject_id.mat')

load('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/run_name.mat')

%directory containing timeseries on server
cd '/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/timeseries'
for ID = 1:70
    
%file name example 'sub-01_run-1_timeseries.mat'
  if ID<10
    sub_name = sprintf('%s%d','sub-0',subject_id(ID,1));
  else
    sub_name = sprintf('%s%d','sub-',subject_id(ID,1));
  end
%   %run correspond to scanning session align to face onset times
  if run_name(ID,1)== 1
      run_newname= 2;
  else run_name(ID,1)== 2
      run_newname= 1;
  end

% run = sprintf('%s%d','run-',run_name(ID)); %faces
run = sprintf('%s%d','run-',run_newname); %places

% run = sprintf('%s%d','run-',run_newname); %for place onset

%filename of timeseries
file_name = sprintf('%s%s%s%s%s',sub_name,'_',run,'_','timeseries.mat');
if isfile(file_name) == 0
 sprintf('%s%s',file_name,' does not exist')
  % Skip to bottom of loop and continue with the loop
  continue;
end

load(file_name);


% standard functional connectivity
ts_corr = corr(ts); %calculate the Pearson's correlation between regions
nPairs = 502; %number of unique pairs from template
nTime = size(ts,1); %number of time points
he
%dsmtx truncated for length of ts

dsmtx_trunct = All_dsmtx(ID).dsmtx_faces;

x = [nTime+1:size(dsmtx_trunct,1)]';

dsmtx_trunct(x,:) = [];

for qq = 1:nPairs
beta_ts(:,qq) = glmfit(dsmtx_trunct,ts(:,qq)); 
%sprintf('%d',qq)
end

%save into data structure

%faces
% All_glm_fc(ID) = struct('sub',{sub_name},'fc_faces',{ts_corr},'glm_faces',{beta_ts});

%places
All_glm_fc(ID).glm_places = beta_ts; 
All_glm_fc(ID).fc_places = ts_corr;
All_glm_fc(ID).sub = sub_name; 

clear nTime clear ts clear dsmtx_trunct x ts_corr

end