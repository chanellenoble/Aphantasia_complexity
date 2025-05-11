for ID=1:70
  if ID<10
    sub_name = sprintf('%s%d','sub-0',subject_id(ID,1));
  else
    sub_name = sprintf('%s%d','sub-',subject_id(ID,1));
  end

%run correspond to scanning session align to face onset times
run = sprintf('%s%d','run-',run_name(ID,1));
        
    
% ID - needs to correspond to subject
dsmtx_image = dsmtx_hrf_epoch(All_onsets(ID).face_onset(:,1),durations(:,1),3);

dsmtx_control = dsmtx_hrf_epoch(All_onsets(ID).face_onset(:,2),durations(:,2),3);

dsmtx_imagine = dsmtx_hrf_epoch(All_onsets(ID).face_onset(:,3),durations(:,3),3);

%Add zeroes to each dsmtx so they are all the same length
M = max([length(dsmtx_image) length(dsmtx_control) length(dsmtx_imagine)]);


if length(dsmtx_image) < M
    x = [length(dsmtx_image):M];
    dsmtx_image(x) = 0;
end

if length(dsmtx_control) < M
    x = [length(dsmtx_control):M];
    dsmtx_control(x) = 0;
end

if length(dsmtx_imagine) < M
    x = [length(dsmtx_imagine):M];
    dsmtx_imagine(x) = 0;
end

%dsmtx that contains different dsmtx for each condition, column 1= image
%dsmtx, col 2 = control dsmtx, col 3 = imagine dsmtx

dsmtx_all(:,1) = dsmtx_image;

dsmtx_all(:,2) = dsmtx_control;

dsmtx_all(:,3) = dsmtx_imagine;

%structural array that saves all dsmtx variables into one
All_dsmtx(ID) = struct('subject',{sub_name},'run',{run},'dsmtx_faces',{dsmtx_all});

clear dsmtx_all dsmtx_image dsmtx_control dsmtx_imagine

end


% added design matrix for place onset times
for ID=1:70   
    % ID - needs to correspond to subject
dsmtx_image = dsmtx_hrf_epoch(All_onsets(ID).place_onset(:,1),durations(:,1),3);

dsmtx_control = dsmtx_hrf_epoch(All_onsets(ID).place_onset(:,2),durations(:,2),3);

dsmtx_imagine = dsmtx_hrf_epoch(All_onsets(ID).place_onset(:,3),durations(:,3),3);

%Add zeroes to each dsmtx so they are all the same length
M = max([length(dsmtx_image) length(dsmtx_control) length(dsmtx_imagine)]);


if length(dsmtx_image) < M
    x = [length(dsmtx_image):M];
    dsmtx_image(x) = 0;
end

if length(dsmtx_control) < M
    x = [length(dsmtx_control):M];
    dsmtx_control(x) = 0;
end

if length(dsmtx_imagine) < M
    x = [length(dsmtx_imagine):M];
    dsmtx_imagine(x) = 0;
end

%dsmtx that contains different dsmtx for each condition, column 1= image
%dsmtx, col 2 = control dsmtx, col 3 = imagine dsmtx

dsmtx_all(:,1) = dsmtx_image;

dsmtx_all(:,2) = dsmtx_control;

dsmtx_all(:,3) = dsmtx_imagine;

%structural array that saves all dsmtx variables into one
%All_dsmtx(ID) = struct('subject',{sub_name},'run',{run},'dsmtx_faces',{dsmtx_all});
All_dsmtx(ID).dsmtx_places = dsmtx_all;

clear dsmtx_all dsmtx_image dsmtx_control dsmtx_imagine

end

% %plot of dsmtx
% figure
% plot(dsmtx_image)
% hold on
% plot(dsmtx_control)
% hold on
% plot(dsmtx_imagine)
% legend('dsmtx image', 'dsmtx control', 'dsmtx imagine')