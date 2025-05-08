%Cleans preprocessed fMRI data and saves into a data structure


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
