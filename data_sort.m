% %Saves preprocessed fMRI dta into a data structure 
% 
%% OUTPUT: 
% 
%  1x70 data structure with 4 fields corresponding to subject ID, Run ID
%  (aligned to faces), and then two fields for the onset times for each block in task for face runs and for place runs. 



for ID=1:70
  % grows by 3 with each loop
  a = 1:3:210;
  c = a(ID)+2;

  %subject naming 
  if ID<10
    sub_name = sprintf('%s%d','sub-0',subject_id(ID,1));
  else
    sub_name = sprintf('%s%d','sub-',subject_id(ID,1));
  end
  %run correspond to scanning session align to face onset times
run = sprintf('%s%d','run-',run_name(ID,1));
All_onsets(ID) = struct('subject',{sub_name},'run',{run},'face_onset',{onsets(:,a(1,ID):c)},'place_onset',{onsets_places(:,a(1,ID):c)}); %changed from 'face_onset' to 'place_onset'...creates a variable 
clear c sub_name run
end
