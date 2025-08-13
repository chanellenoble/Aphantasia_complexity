% Loads BOLD time series for each subject, calculates MTD score between
% ROIs and then fits to GLM using HRF confovled DSMTX.
%
%% NB -- This file requires the following to work:
% 'flatten_mat.m'
% 'perm_code.m'
% 'coupling_pos.m' by Mac Shine (available at https://github.com/macshine/coupling)
%
% Reference: Shine JM, Koyejo O, Bell PT, Gorgolewski KJ, Gilat M, Poldrack RA. Estimation of dynamic functional connectivity using Multiplication of Temporal Derivatives. Neuroimage. 2015 Nov 15;122:399-407. doi: 10.1016/j.neuroimage.2015.07.064. Epub 2015 Jul 29. PMID: 26231247.


%Set current working directory - Make sure to include folder containing
%extra scripts needed
    cd 'Z:\PRJ-shine_hpc\Aphantasia_fMRI'
    addpath 'C:\Users\Chanelle\Documents\working\Scripts'


%Load dsmtx, list of subject IDs, and list of run names for subjet's first
%scan aligned to face sessions (i.e. if run_name = 1 --> patient's first
%scan was their face session)

    load('All_dsmtx.mat')
    load('idx_ach.mat');
    load('subject_id.mat')
    load('run_name.mat')


% predefine some variables
    nSubs = 70; %Number of subjects
    nNodes = 502; %Number of ROIs
    nEdges = (nNodes*nNodes-nNodes)/2;

% make a template that has a '1' every time there is a unique edge, and a
% '0' every time there is a redundant edge
    template = tril(ones(nNodes))-eye(nNodes);
    template = find(template);


%% ---- MTD of time series ----

%directory containing time series on server
    cd 'C:\Users\Chanelle\Documents\working\timeseries'


% Faces
    

    for ID = 1:nSubs
        
    %file name example 'sub-01_run-1_timeseries.mat'
      if ID<10
        sub_name = sprintf('%s%d','sub-0',subject_id(ID,1));
      else
        sub_name = sprintf('%s%d','sub-',subject_id(ID,1));
      end

    %run correspond to scanning session align to face onset times
      if run_name(ID,1)== 1
          run_newname= 2;
      else run_name(ID,1)== 2
          run_newname= 1;
      end
    
     run = sprintf('%s%d','run-',run_name(ID));
        
    %filename of time series
        file_name = sprintf('%s%s%s%s%s',sub_name,'_',run,'_','timeseries.mat');
        if isfile(file_name) == 0   % If file does not exist, skip file and move to next subject
          sprintf('%s%s',file_name,' does not exist');
          continue;
        end
        
        load(file_name);
        data = ts;
        nTR = size(data,1);
    
    %dsmtx truncated for length of ts
        dsmtx_trunct = All_dsmtx(ID).dsmtx_faces;
        x = [nTR+1:size(dsmtx_trunct,1)]';
        dsmtx_trunct(x,:) = [];
    
    % run MTD - note: coupling_pos is a newer version of the code that only
        % considers coupling when two nodes are both increasing together, so as
        % to avoid the weird cases when there is positive coupling between two
        % nodes that are dropping their activity...
        mtd = coupling_pos(data,10);
    % flatten MTD - should be nEdges x length of the scan
        mtd_flat = zeros(nEdges,nTR);
        
        for tt = 1:nTR
            mtd_temp = mtd(:,:,tt);
            mtd_flat(:,tt) = mtd_temp(template);
        end
        
    % load in the unique design matrix for the run/subject
        dsmtx1 = dsmtx_trunct;
        for jj = 1:nEdges
            mtd_beta(jj,:,ID) = glmfit(dsmtx1,mtd_flat(jj,:)'); % GLM predicting the effect of task on MTD score
        end

    
        nBeta = size(dsmtx_trunct,2)+1;
        for bb = 1:nBeta
            mtd_beta_mtx_faces(:,:,bb,ID) = squareform(mtd_beta(:,bb,ID)); % Turn flattened MTD_beta vector back into a matrix
        end
        
    % run a ticker to track how many subjects have finished
        sprintf('%s%s',file_name,' completed')
   
    
        clear ts dsmtx_trunct x data


% Places

     run = sprintf('%s%d','run-',run_newname); %places
    
    %filename of timeseries
        file_name = sprintf('%s%s%s%s%s',sub_name,'_',run,'_','timeseries.mat');

        if isfile(file_name) == 0
          sprintf('%s%s',file_name,' does not exist');
          continue;
        end

        load(file_name);
        data = ts;
        nTR = size(data,1);

    %dsmtx truncated for length of ts
    
        dsmtx_trunct = All_dsmtx(ID).dsmtx_places;
        
        x = [nTR+1:size(dsmtx_trunct,1)]';
        
        dsmtx_trunct(x,:) = [];
    
     % run MTD - note: coupling_pos is a newer version of the code that only
        % considers coupling when two nodes are both increasing together, so as
        % to avoid the weird cases when there is positive coupling between two
        % nodes that are dropping their activity...
        mtd = coupling_pos(data,10);
    % flatten MTD - should be nEdges x length of the scan
        mtd_flat = zeros(nEdges,nTR);
        
        for tt = 1:nTR
            mtd_temp = mtd(:,:,tt);
            mtd_flat(:,tt) = mtd_temp(template);
        end
        
    % load in the unique design matrix for the run/subject
        dsmtx1 = dsmtx_trunct;
        for jj = 1:nEdges
            mtd_beta(jj,:,ID) = glmfit(dsmtx1,mtd_flat(jj,:)');
        end
    
        nBeta = size(dsmtx_trunct,2)+1;
        for bb = 1:nBeta
            mtd_beta_mtx_places(:,:,bb,ID) = squareform(mtd_beta(:,bb,ID));
        end
        
        % run a ticker to track how many subjects have finished
        sprintf('%s%s',file_name,' completed')
    
  
    
    clear dsmtx_trunct x  data
    end


%% ---- Separating Subject Groups ----


    cd 'Z:\PRJ-shine_hpc\Aphantasia_fMRI\Sripts'

% Create variables for aphantasics, controls, and hyperphantasics

    aphan = {'01','02','04','05','06','07','08','09','10','12','13','14','16','18','24','35','50','60','66','68','70','71','73'};
    
    ctrl = {'17','19','21','22','29','30','31','33','36','38','39','41','43','44','45','46','54','56','57','59'};
    
    hyper = {'23','25','26','27','28','32','34','37','40','42','47','48','49','52','53','55','58','61','62','63','64','65','67'};
    
    ach = horzcat(aphan,ctrl,hyper);


    % Aphan

        mtd_beta_mtx_faces_aphan = mtd_beta_mtx_faces(:,:,:,idx_ach(1:size(aphan,2)));
    
        mtd_beta_mtx_places_aphan = mtd_beta_mtx_places(:,:,:,idx_ach(1:size(aphan,2)));
    
        mtd_beta_mtx_aphan = cat(4, mtd_beta_mtx_faces_aphan, mtd_beta_mtx_places_aphan);
    
    
        % Remove missing matrices
    
            % Create a logical vector to mark valid subjects
            num_subs = size(mtd_beta_mtx_aphan, 4);
            valid_subs = false(num_subs, 1);
            for ii = 1:num_subs
                mat = mtd_beta_mtx_aphan(:,:,:,ii);
                if any(mat(:) ~= 0)  % if there's any non-zero entry
                    valid_subs(ii) = true;
                end
            end
    
           mtd_beta_mtx_aphan = mtd_beta_mtx_aphan(:,:,:,valid_subs);
        
           clear num_subs valid_subs

    % Ctrl

        mtd_beta_mtx_faces_ctrl = mtd_beta_mtx_faces(:,:,:,idx_ach(size(aphan,2)+1:size(aphan,2)+size(ctrl,2)));
    
        mtd_beta_mtx_places_ctrl = mtd_beta_mtx_places(:,:,:,idx_ach(size(aphan,2)+1:size(aphan,2)+size(ctrl,2)));
    
        mtd_beta_mtx_ctrl = cat(4, mtd_beta_mtx_faces_ctrl, mtd_beta_mtx_places_ctrl);
    
    
        % Remove missing matrices
    
            % Create a logical vector to mark valid subjects
            num_subs = size(mtd_beta_mtx_ctrl, 4);
            valid_subs = false(num_subs, 1);
            for ii = 1:num_subs
                mat = mtd_beta_mtx_ctrl(:,:,:,ii);
                if any(mat(:) ~= 0)  % if there's any non-zero entry
                    valid_subs(ii) = true;
                end
            end
    
            
            mtd_beta_mtx_ctrl = mtd_beta_mtx_ctrl(:,:,:,valid_subs);
        
           clear num_subs valid_subs


     % Hyper

    mtd_beta_mtx_faces_hyper = mtd_beta_mtx_faces(:,:,:,idx_ach(size(aphan,2)+size(ctrl,2)+1:end));

    mtd_beta_mtx_places_hyper = mtd_beta_mtx_places(:,:,:,idx_ach(size(aphan,2)+size(ctrl,2)+1:end));

    mtd_beta_mtx_hyper = cat(4, mtd_beta_mtx_faces_hyper, mtd_beta_mtx_places_hyper);


    % Remove missing matrices

        % Create a logical vector to mark valid subjects
        num_subs = size(mtd_beta_mtx_hyper, 4);
        valid_subs = false(num_subs, 1);
        for ii = 1:num_subs
            mat = mtd_beta_mtx_hyper(:,:,:,ii);
            if any(mat(:) ~= 0)  % if there's any non-zero entry
                valid_subs(ii) = true;
            end
        end

       mtd_beta_mtx_hyper = mtd_beta_mtx_hyper(:,:,:,valid_subs);
    
       clear num_subs valid_subs




%% ---- Permutation Testing ----



    % Aphan minus Ctrl - Imagery
    
        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(mtd_beta_mtx_aphan,[],4,size(mtd_beta_mtx_aphan,4)); %Flattens regions
        test2 = reshape(mtd_beta_mtx_ctrl,[],4,size(mtd_beta_mtx_ctrl,4));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);

        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
            data1 = squeeze(test_lower(ii,4,:));
            data2 = squeeze(test2_lower(ii,4,:));
            [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_imagine_aphan_ctrl(1,:) = sig;
        perm_test_imagine_aphan_ctrl(2,:) = pval;
    
%         save('perm_test_imagine_aphan_ctrl.mat', 'perm_test_imagine_aphan_ctrl')
% 
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_imagine_aphan_ctrl.mat', 'perm_test_imagine_aphan_ctrl');



     % Aphan minus Ctrl - Perception
    
%         %creates matrix of lower triangle of full ROIxROI matrix
%         test = reshape(mtd_beta_mtx_aphan,[],4,size(mtd_beta_mtx_aphan,4)); %Flattens regions
%         test2 = reshape(mtd_beta_mtx_ctrl,[],4,size(mtd_beta_mtx_ctrl,4));
%         test_lower = test(template,:,:); %lower triangle
%         test2_lower = test2(template,:,:);
% 
%         sig = zeros(1, size(test_lower, 1));
%         pval = zeros(1, size(test_lower, 1));
%         
%         for ii = 1:size(test_lower,1)
%             data1 = squeeze(test_lower(ii,2,:));
%             data2 = squeeze(test2_lower(ii,2,:));
%             [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
%         end
%         
%         perm_test_perc_aphan_ctrl(1,:) = sig;
%         perm_test_perc_aphan_ctrl(2,:) = pval;
%         
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_perc_aphan_ctrl.mat', 'perm_test_perc_aphan_ctrl');



    % Hyper vs Ctrl - Imagery
    
        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(mtd_beta_mtx_hyper,[],4,size(mtd_beta_mtx_hyper,4)); %Flattens regions
        test2 = reshape(mtd_beta_mtx_ctrl,[],4,size(mtd_beta_mtx_ctrl,4));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);

        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
            data1 = squeeze(test_lower(ii,4,:));
            data2 = squeeze(test2_lower(ii,4,:));
            [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_imagine_hyper_ctrl(1,:) = sig;
        perm_test_imagine_hyper_ctrl(2,:) = pval;
    
        
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_imagine_hyper_ctrl.mat', 'perm_test_imagine_hyper_ctrl');



     % Hyper minus Ctrl - Perception
    
%         %creates matrix of lower triangle of full ROIxROI matrix
%         test = reshape(mtd_beta_mtx_hyper,[],4,size(mtd_beta_mtx_hyper,4)); %Flattens regions
%         test2 = reshape(mtd_beta_mtx_ctrl,[],4,size(mtd_beta_mtx_ctrl,4));
%         test_lower = test(template,:,:); %lower triangle
%         test2_lower = test2(template,:,:);
%         
%         sig = zeros(1, size(test_lower, 1));
%         pval = zeros(1, size(test_lower, 1));
%         
%         for ii = 1:size(test_lower,1)
%             data1 = squeeze(test_lower(ii,2,:));
%             data2 = squeeze(test2_lower(ii,2,:));
%             [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
%         end
%         
%         perm_test_perc_hyper_ctrl(1,:) = sig;
%         perm_test_perc_hyper_ctrl(2,:) = pval;
%     
        
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_perc_hyper_ctrl.mat', 'perm_test_perc_hyper_ctrl');



% Hyper vs Aphan - Imagery
    
        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(mtd_beta_mtx_hyper,[],4,size(mtd_beta_mtx_hyper,4)); %Flattens regions
        test2 = reshape(mtd_beta_mtx_aphan,[],4,size(mtd_beta_mtx_aphan,4));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);
        
        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
            data1 = squeeze(test_lower(ii,4,:));
            data2 = squeeze(test2_lower(ii,4,:));
            [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_imagine_hyper_aphan(1,:) = sig;
        perm_test_imagine_hyper_aphan(2,:) = pval;
    
       
%                 save ('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/Scripts/perm_test_imagine_hyper_aphan.mat', 'perm_test_imagine_hyper_aphan');
% 
% 
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_imagine_hyper_aphan.mat', 'perm_test_imagine_hyper_aphan');



% Hyper minus Aphan - Perception
    
%         %creates matrix of lower triangle of full ROIxROI matrix
%         test = reshape(mtd_beta_mtx_hyper,[],4,size(mtd_beta_mtx_hyper,4)); %Flattens regions
%         test2 = reshape(mtd_beta_mtx_aphan,[],4,size(mtd_beta_mtx_aphan,4));
%         test_lower = test(template,:,:); %lower triangle
%         test2_lower = test2(template,:,:);
%         
%         sig = zeros(1, size(test_lower, 1));
%         pval = zeros(1, size(test_lower, 1));
%         
%         for ii = 1:size(test_lower,1)
%             data1 = squeeze(test_lower(ii,2,:));
%             data2 = squeeze(test2_lower(ii,2,:));
%             [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
%         end
%         
%         perm_test_perc_hyper_aphan(1,:) = sig;
%         perm_test_perc_hyper_aphan(2,:) = pval;
    
%         save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_perc_hyper_aphan.mat', 'perm_test_perc_hyper_aphan');





%% ---- Calculating Mean Across Schaeffer 7-networks and Temporal Lobes ----


    % Load network labels for each cortical ROI 
        load('..\schaef_id14.mat'); 
        % 1 = visual
        % 2 = motor
        % 3 = dorsal attention
        % 4 = ventral attention
        % 5 = frontotemporal
        % 6 = frontoparietal
        % 7 = default


    load('..\temp_lobe_mask.mat')
        
    new_id = schaef_id14+2; %Create variable new_id with leaving two rows for temporal lobe
    
    new_id(temp_lobe_mask==1) = 1; %Apply temp_lobe_mask
    
    temp = temp_lobe_mask(1:200);%Left temporal lobe
    
    temp(201:400) = temp_lobe_mask(201:400).*2; %Right temporal lobe
    
    new_id(temp_lobe_mask>0) = temp(temp_lobe_mask>0);
    
    unique(new_id);
    
    [order_new,sort_new] = sort(new_id); 

    order_new_notemp = order_new(find(temp_lobe_mask == 0));



   % Aphan minus Ctrl - Imagery

       mean_mtd_beta_aphan_imagery = mean(mtd_beta_mtx_aphan(:,:,4,:),4); % Calculate mean MTD betas across subjects
       mean_mtd_beta_ctrl_imagery = mean(mtd_beta_mtx_ctrl(:,:,4,:),4);

       delta_mtd_beta_aphan_ctrl = mean_mtd_beta_aphan_imagery - mean_mtd_beta_ctrl_imagery; % aphan minus ctrl

       clear sig

       sig = squareform(perm_test_both_imagine(1,:)); % Create mask ROI w significant differences between aphan and ctrl
       delta_mtd_beta_aphan_ctrl = delta_mtd_beta_aphan_ctrl.*sig;

       delta_mtd_beta_aphan_ctrl_temp2rsn = delta_mtd_beta_aphan_ctrl(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_aphan_ctrl_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_aphan_ctrl_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_aphan_ctrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_aphan_ctrl_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_aphan_ctrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_aphan_ctrl_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end



    % Hyper minus Ctrl - Imagery

       mean_mtd_beta_hyper_imagery = mean(mtd_beta_mtx_hyper(:,:,4,:),4); % Calculate mean MTD betas across subjects
       mean_mtd_beta_ctrl_imagery = mean(mtd_beta_mtx_ctrl(:,:,4,:),4);

       delta_mtd_beta_hyper_ctrl =  mean_mtd_beta_hyper_imagery - mean_mtd_beta_ctrl_imagery; % hyper minus ctrl

       clear sig

       sig = squareform(perm_test_imagine_hyper_ctrl(1,:)); % Create mask ROI w significant differences between hyper and ctrl
       delta_mtd_beta_hyper_ctrl = delta_mtd_beta_hyper_ctrl.*sig;

       delta_mtd_beta_hyper_ctrl_temp2rsn = delta_mtd_beta_hyper_ctrl(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_hyper_ctrl_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_hyper_ctrl_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_hyper_ctrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_ctrl_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_hyper_ctrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_ctrl_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end




% Hyper minus Aphan - Imagery

       mean_mtd_beta_hyper_imagery = mean(mtd_beta_mtx_hyper(:,:,4,:),4); % Calculate mean MTD betas across subjects
       mean_mtd_beta_aphan_imagery = mean(mtd_beta_mtx_aphan(:,:,4,:),4);

       delta_mtd_beta_hyper_aphan =  mean_mtd_beta_hyper_imagery - mean_mtd_beta_aphan_imagery; % hyper minus aphan

       clear sig

       sig = squareform(perm_test_imagine_hyper_aphan(1,:)); % Create mask ROI w significant differences between hyper and aphan
       delta_mtd_beta_hyper_aphan = delta_mtd_beta_hyper_aphan.*sig;

       delta_mtd_beta_hyper_aphan_temp2rsn = delta_mtd_beta_hyper_aphan(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_hyper_aphan_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_hyper_aphan_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_hyper_aphan_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_aphan_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_hyper_aphan_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_aphan_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end


































%% ---- Perm test Deltas ----            

% Delta Aphan vs Delta Ctrl
    
    delta_aphan_mtd_beta = mtd_beta_mtx_aphan(:,:,4,:) - mtd_beta_mtx_aphan(:,:,2,:);
    delta_ctrl_mtd_beta = mtd_beta_mtx_ctrl(:,:,4,:) - mtd_beta_mtx_ctrl(:,:,2,:);

    delta_aphan_mtd_beta = squeeze(delta_aphan_mtd_beta);
    delta_ctrl_mtd_beta = squeeze(delta_ctrl_mtd_beta);


        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(delta_aphan_mtd_beta,[],size(delta_aphan_mtd_beta,3)); %Flattens regions
        test2 = reshape(delta_ctrl_mtd_beta,[],size(delta_ctrl_mtd_beta,3));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);

        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
            data1 = squeeze(test_lower(ii,:))';
            data2 = squeeze(test2_lower(ii,:))';
            [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_delta_aphan_ctrl_mtd_beta(1,:) = sig;
        perm_test_delta_aphan_ctrl_mtd_beta(2,:) = pval;
    

        save ('Z:\PRJ-shine_hpc\Aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_delta_aphan_ctrl_mtd_beta.mat', 'perm_test_delta_aphan_ctrl_mtd_beta');


% Delta hyper vs Delta Ctrl
    
    delta_hyper_mtd_beta = mtd_beta_mtx_hyper(:,:,4,:) - mtd_beta_mtx_hyper(:,:,2,:);
    delta_ctrl_mtd_beta = mtd_beta_mtx_ctrl(:,:,4,:) - mtd_beta_mtx_ctrl(:,:,2,:);

    delta_hyper_mtd_beta = squeeze(delta_hyper_mtd_beta);
    delta_ctrl_mtd_beta = squeeze(delta_ctrl_mtd_beta);


        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(delta_hyper_mtd_beta,[],size(delta_hyper_mtd_beta,3)); %Flattens regions
        test2 = reshape(delta_ctrl_mtd_beta,[],size(delta_ctrl_mtd_beta,3));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);

        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
            data1 = squeeze(test_lower(ii,:))';
            data2 = squeeze(test2_lower(ii,:))';
            [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_delta_hyper_ctrl_mtd_beta(1,:) = sig;
        perm_test_delta_hyper_ctrl_mtd_beta(2,:) = pval;
    

        save ('Z:\PRJ-shine_hpc\hypertasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_delta_hyper_ctrl_mtd_beta.mat', 'perm_test_delta_hyper_ctrl_mtd_beta');


    % Delta hyper vs Delta aphan
    
    delta_hyper_mtd_beta = mtd_beta_mtx_hyper(:,:,4,:) - mtd_beta_mtx_hyper(:,:,2,:);
    delta_aphan_mtd_beta = mtd_beta_mtx_aphan(:,:,4,:) - mtd_beta_mtx_aphan(:,:,2,:);

    delta_hyper_mtd_beta = squeeze(delta_hyper_mtd_beta);
    delta_aphan_mtd_beta = squeeze(delta_aphan_mtd_beta);


        %creates matrix of lower triangle of full ROIxROI matrix
        test = reshape(delta_hyper_mtd_beta,[],size(delta_hyper_mtd_beta,3)); %Flattens regions
        test2 = reshape(delta_aphan_mtd_beta,[],size(delta_aphan_mtd_beta,3));
        test_lower = test(template,:,:); %lower triangle
        test2_lower = test2(template,:,:);

        sig = zeros(1, size(test_lower, 1));
        pval = zeros(1, size(test_lower, 1));
        
        for ii = 1:size(test_lower,1)
           data1 = squeeze(test_lower(ii,:))';
           data2 = squeeze(test2_lower(ii,:))';
           [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
        end
        
        perm_test_delta_hyper_aphan_mtd_beta(1,:) = sig;
        perm_test_delta_hyper_aphan_mtd_beta(2,:) = pval;
    

        save ('Z:\PRJ-shine_hpc\aphantasia_fMRI\mtd_output\MTD_Perm_Test_25-07-16\perm_test_delta_hyper_aphan_mtd_beta.mat', 'perm_test_delta_hyper_aphan_mtd_beta');










         % Delta Aphan minus Delta Ctrl

       delta_aphan_delta_ctrl = mean(delta_aphan_mtd_beta(:,:,:),3) - mean(delta_ctrl_mtd_beta(:,:,:),3); % aphan minus ctrl

       clear sig

       sig = squareform(perm_test_delta_aphan_ctrl_mtd_beta(1,:)); % Create mask ROI w significant differences between aphan and ctrl
       delta_aphan_delta_ctrl = delta_aphan_delta_ctrl.*sig;

       delta_mtd_beta_aphan_ctrl_temp2rsn = delta_aphan_delta_ctrl(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_deltaaphan_deltactrl_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_deltaaphan_deltactrl_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_deltaaphan_deltactrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_aphan_ctrl_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_deltaaphan_deltactrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_aphan_ctrl_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end










         % Delta hyper minus Delta Ctrl

       delta_hyper_delta_ctrl = mean(delta_hyper_mtd_beta(:,:,:),3) - mean(delta_ctrl_mtd_beta(:,:,:),3); % hyper minus ctrl

       clear sig

       sig = squareform(perm_test_delta_hyper_ctrl_mtd_beta(1,:)); % Create mask ROI w significant differences between hyper and ctrl
       delta_hyper_delta_ctrl = delta_hyper_delta_ctrl.*sig;

       delta_mtd_beta_hyper_ctrl_temp2rsn = delta_hyper_delta_ctrl(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_deltahyper_deltactrl_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_deltahyper_deltactrl_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_deltahyper_deltactrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_ctrl_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_deltahyper_deltactrl_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_ctrl_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end





         % Delta hyper minus Delta aphan

       delta_hyper_delta_aphan = mean(delta_hyper_mtd_beta(:,:,:),3) - mean(delta_aphan_mtd_beta(:,:,:),3); % hyper minus aphan

       clear sig

       sig = squareform(perm_test_delta_hyper_aphan_mtd_beta(1,:)); % Create mask ROI w significant differences between hyper and aphan
       delta_hyper_delta_aphan = delta_hyper_delta_aphan.*sig;

       delta_mtd_beta_hyper_aphan_temp2rsn = delta_hyper_delta_aphan(order_new<3, order_new>2);

       %calculate left and right side temporal lobe avg. MTD for each
       %Resting-state network
            lt_delta_mtd_beta_deltahyper_deltaaphan_temp2rsn = zeros(14,1);
            rt_delta_mtd_beta_deltahyper_deltaaphan_temp2rsn = zeros(14,1);
            
            
            
            for xx = 1:14
            net_id = xx+2;
            lt_delta_mtd_beta_deltahyper_deltaaphan_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_aphan_temp2rsn(1:33,find(order_new_notemp==net_id))));
            
            rt_delta_mtd_beta_deltahyper_deltaaphan_temp2rsn(xx,1) = mean(mean(delta_mtd_beta_hyper_aphan_temp2rsn(34:66,find(order_new_notemp==net_id))));
            end


        
