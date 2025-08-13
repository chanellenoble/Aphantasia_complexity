% Set wd

    cd I:\Aphantasia_fMRI\Scripts

% Load in dsmtx and FC  (Pearson's Correlation)

    load('All_glm_fc.mat');

    load('All_onsets.mat');
    load('run_name.mat');
    load('subject_id.mat')
    load('template.mat');

    
    addpath 'functions\'
    load('temp_mask.mat');

    cd I:\Aphantasia_fMRI\Scripts\Complexity

    load idx_ach.mat
    load vviq.mat
 

% Predefine some variables
    
    nSubs_orig = 70; % number of subjects in original study
    nSubs = 66; % number of subjects
    nNodes = 502; % number of ROIs
    TR = 3;
    nEdges = (nNodes*nNodes-nNodes)/2;



% Create variables for aphantasics, controls, and hyperphantasics

    aphan = {'01','02','04','05','06','07','08','09','10','12','13','14','16','18','24','35','50','60','66','68','70','71','73'};
    
    ctrl = {'17','19','21','22','29','30','31','33','36','38','39','41','43','44','45','46','54','56','57','59'};
    
    hyper = {'23','25','26','27','28','32','34','37','40','42','47','48','49','52','53','55','58','61','62','63','64','65','67'};
    
    ach = horzcat(aphan,ctrl,hyper);



% ---- Calculate FC for imagery and perception for each group ----

cd I:\Aphantasia_fMRI\timeseries



% Calculate stats faces run

for ID = 1:nSubs_orig

    file_name = f_name(ID,subject_id,run_name); % loads file name. 1 = faces; 2 = places.

    if isfile(file_name) == 0;   % If file does not exist, skip file and move to next subject
       sprintf('%s%s',file_name,' does not exist');
       continue;
    end

    load(file_name); % load time series for faces run


    % Perception

        x = 1;    
        for ii = 1:4:length(All_onsets(1).face_onset);
         onset_temp(x,1) = All_onsets(ID).face_onset(ii,1);
         onset_temp(x,2) = All_onsets(ID).face_onset(ii,2);
         x = x + 1;
        end

        onset_temp = round(onset_temp/TR);

        data_temp = [];
        
        for ii = 1:length(onset_temp);
            if onset_temp(ii,1) > size(ts,1);
               continue
            else
                if onset_temp(ii,2) > size(ts,1);
                    x = size(data_temp,1) + 1;
                    y = x + (size(ts,1) - onset_temp(ii,1));
                    data_temp(x:y,:) = ts(onset_temp(ii,1):end,:);
                else
                    if isempty(data_temp);
                        x  = 1;
                        y = length(onset_temp(ii,1):onset_temp(ii,2));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    else
                        x = size(data_temp,1) + 1;
                        y = (size(data_temp,1) + 1) + (onset_temp(ii,2) - onset_temp(ii,1));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    end
                end
            end
        end



        for ii = 1:size(data_temp,2)
            var_perc_faces(ii,ID) = std(data_temp(:,ii)); % Calculate the standard deviation for each ROI in temp_lobe_mask
        end


        for ii = 1:size(data_temp,2)
            [acf] = autocorr(data_temp(:,ii),1); % Calculate the variance for each ROI in temp_lobe_mask
            ac_perc_faces(ii,ID) = acf(2);
        end

        clear data_temp onset_temp;



    % Imagery

        x = 1;    
        for ii = 1:4:length(All_onsets(1).face_onset);
            if ii < 37;
             onset_temp(x,1) = All_onsets(ID).face_onset(ii,3);
             onset_temp(x,2) = All_onsets(ID).face_onset(ii+4,1);
             x = x + 1;
            else
             onset_temp(x,1) = All_onsets(ID).face_onset(ii,3);
             onset_temp(x,2) = All_onsets(ID).face_onset(ii+3,3);
            end
        end


        onset_temp = round(onset_temp/TR);

        data_temp = [];

         for ii = 1:length(onset_temp);
            if onset_temp(ii,1) > size(ts,1);
               continue
            else
                if onset_temp(ii,2) > size(ts,1);
                    x = size(data_temp,1) + 1;
                    y = x + (size(ts,1) - onset_temp(ii,1));
                    data_temp(x:y,:) = ts(onset_temp(ii,1):end,find(temp_lobe_mask == 1));
                else
                    if isempty(data_temp);
                        x  = 1;
                        y = length(onset_temp(ii,1):onset_temp(ii,2));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    else
                        x = size(data_temp,1) + 1;
                        y = (size(data_temp,1) + 1) + (onset_temp(ii,2) - onset_temp(ii,1));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    end
                end
            end
        end



        for ii = 1:size(data_temp,2)
            var_imagery_faces(ii,ID) = std(data_temp(:,ii)); % Calculate the standard deviation for each ROI in temp_lobe_mask
        end


        for ii = 1:size(data_temp,2)
            [acf] = autocorr(data_temp(:,ii),1); % Calculate the variance for each ROI in temp_lobe_mask
            ac_imagery_faces(ii,ID) = acf(2);
        end

        clear data_temp onset_temp;

    end






% Calculate stats places run

for ID = 1:nSubs_orig

    [~,file_name] = f_name(ID,subject_id,run_name); % loads file name. 1 = faces; 2 = places.

    if isfile(file_name) == 0;   % If file does not exist, skip file and move to next subject
       sprintf('%s%s',file_name,' does not exist');
       continue;
    end
    load(file_name); % load time series for faces run


    % Perception

        x = 1;    
        for ii = 1:4:length(All_onsets(1).place_onset);
         onset_temp(x,1) = All_onsets(ID).place_onset(ii,1);
         onset_temp(x,2) = All_onsets(ID).place_onset(ii,2);
         x = x + 1;
        end

        onset_temp = round(onset_temp/TR);

        data_temp = [];
        
        for ii = 1:length(onset_temp);
            if onset_temp(ii,1) > size(ts,1);
               continue
            else
                if onset_temp(ii,2) > size(ts,1);
                    x = size(data_temp,1) + 1;
                    y = x + (size(ts,1) - onset_temp(ii,1));
                    data_temp(x:y,:) = ts(onset_temp(ii,1):end,find(temp_lobe_mask == 1));
                else
                    if isempty(data_temp);
                        x  = 1;
                        y = length(onset_temp(ii,1):onset_temp(ii,2));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    else
                        x = size(data_temp,1) + 1;
                        y = (size(data_temp,1) + 1) + (onset_temp(ii,2) - onset_temp(ii,1));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    end
                end
            end
        end



        for ii = 1:size(data_temp,2)
            var_perc_places(ii,ID) = std(data_temp(:,ii)); % Calculate the standard deviation for each ROI in temp_lobe_mask
        end


        for ii = 1:size(data_temp,2)
            [acf] = autocorr(data_temp(:,ii),1); % Calculate the variance for each ROI in temp_lobe_mask
            ac_perc_places(ii,ID) = acf(2);
        end

        clear data_temp onset_temp;

    


    % Imagery

        x = 1;    
        for ii = 1:4:length(All_onsets(1).place_onset);
            if ii < 37;
             onset_temp(x,1) = All_onsets(ID).place_onset(ii,3);
             onset_temp(x,2) = All_onsets(ID).place_onset(ii+4,1);
             x = x + 1;
            else
             onset_temp(x,1) = All_onsets(ID).place_onset(ii,3);
             onset_temp(x,2) = All_onsets(ID).place_onset(ii+3,3);
            end
        end


        onset_temp = round(onset_temp/TR);

        data_temp = [];

         for ii = 1:length(onset_temp);
            if onset_temp(ii,1) > size(ts,1);
               continue
            else
                if onset_temp(ii,2) > size(ts,1);
                    x = size(data_temp,1) + 1;
                    y = x + (size(ts,1) - onset_temp(ii,1));
                    data_temp(x:y,:) = ts(onset_temp(ii,1):end,find(temp_lobe_mask == 1));
                else
                    if isempty(data_temp);
                        x  = 1;
                        y = length(onset_temp(ii,1):onset_temp(ii,2));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    else
                        x = size(data_temp,1) + 1;
                        y = (size(data_temp,1) + 1) + (onset_temp(ii,2) - onset_temp(ii,1));
                        data_temp(x:y,:) = ts(onset_temp(ii,1):onset_temp(ii,2),find(temp_lobe_mask == 1));
                    end
                end
            end
        end



        for ii = 1:size(data_temp,2)
            var_imagery_places(ii,ID) = std(data_temp(:,ii)); % Calculate the standard deviation for each ROI in temp_lobe_mask
        end


        for ii = 1:size(data_temp,2)
            [acf] = autocorr(data_temp(:,ii),1); % Calculate the variance for each ROI in temp_lobe_mask
            ac_imagery_places(ii,ID) = acf(2);
        end

        clear data_temp onset_temp;

end

var_perc = horzcat(var_perc_places(:,idx_ach),var_perc_faces(:,idx_ach));
ac_perc = horzcat(ac_perc_places(:,idx_ach),ac_perc_faces(:,idx_ach));

var_imagery = horzcat(var_imagery_places(:,idx_ach),var_imagery_faces(:,idx_ach));
ac_imagery = horzcat(ac_imagery_places(:,idx_ach),ac_imagery_faces(:,idx_ach));


% Build a logical index of valid entries
 % Both must be non-zero

mean_var_perc_sub = mean(var_perc(:,:),1);
mean_var_imagery_sub = mean(var_imagery(:,:),1);


[h1,p1,~,stats1] = ttest(mean_var_imagery_sub,mean_var_perc_sub);

for ii = 1:size(var_perc,1)
    [rho_1_1(ii),p_1_1(ii)] = corr(vviq,var_perc(ii,:)','type','Spearman');
end

sig_1_1 = p_1_1 < 0.05;

for ii = 1:size(var_imagery,1)
    [rho_1_2(ii),p_1_2(ii)] = corr(vviq,var_imagery(ii,:)','type','Spearman');
end

sig_1_2 = p_1_2 < 0.05;






