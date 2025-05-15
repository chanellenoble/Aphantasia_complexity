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


% ------------------------------------------------
%% MTD of time series
% ------------------------------------------------

%directory containing time series on server
    cd 'C:\Users\Chanelle\Documents\working\timeseries'


%%-------------- FACES ------------------------

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
    

%%-------------- PLACES ------------------------ 
    
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
    
    
% ----------------------------------------------------------------------------------------


%%Separating subject groups

%Arrays with the subject numbers for each subject group
    
    aphan = {'01','02','04','05','06','07','08','09','10','12','13','14','16','18','24','35','50','60','66','68','70','71','73'};
    
    ctrl = {'17','19','21','22','29','30','31','33','36','38','39','41','43','44','45','46','54','56','57','59'};
    
    ac = horzcat(aphan, ctrl);

    load('idx_aphan.mat');
    load('idx_ac.mat'); %loads indexes of aphants and ctrls in the full set of 70 subjects according to the arrays. Aphants listed first.
    
    mtd_beta_mtx_faces_aphan = mtd_beta_mtx_faces(:,:,:,idx_aphan);
    mtd_beta_mtx_faces_ctrl = mtd_beta_mtx_faces(:,:,:,idx_ac(:,24:43);
    mtd_beta_mtx_places_aphan = mtd_beta_mtx_places(:,:,:,idx_aphan);
    mtd_beta_mtx_places_ctrl = mtd_beta_mtx_places(:,:,:,idx_ac(:,24:43);
    
    
    save('mtd_beta_mtx_places_aphan.mat','mtd_beta_mtx_places_aphan')
    
    save('mtd_beta_mtx_faces_aphan.mat','mtd_beta_mtx_faces_aphan')
    
    save('mtd_beta_mtx_faces_ctrl.mat','mtd_beta_mtx_faces_ctrl')


%% -------- Calculate means for each condition -------------------

    mean_mtd_beta_aphan_faces_imagine = mean(mtd_beta_mtx_faces_aphan(:,:,4,:),4);
    
    mean_mtd_beta_aphan_places_imagine = mean(mtd_beta_mtx_places_aphan(:,:,4,:),4);
    
    mean_mtd_beta_ctrl_faces_imagine = mean(mtd_beta_mtx_faces_ctrl(:,:,4,:),4);
    
    mean_mtd_beta_ctrl_places_imagine = mean(mtd_beta_mtx_places_ctrl(:,:,4,:),4);
    
    mean_mtd_beta_aphan_faces_image = mean(mtd_beta_mtx_faces_aphan(:,:,2,:),4);
    
    mean_mtd_beta_aphan_places_image = mean(mtd_beta_mtx_places_aphan(:,:,2,:),4);
    
    mean_mtd_beta_ctrl_faces_image = mean(mtd_beta_mtx_faces_ctrl(:,:,2,:),4);
    
    mean_mtd_beta_ctrl_places_image = mean(mtd_beta_mtx_places_ctrl(:,:,2,:),4);
    
    mean_mtd_beta_aphan_both_imagine = (mean_mtd_beta_aphan_places_imagine + mean_mtd_beta_aphan_faces_imagine)/2;
    
    mean_mtd_beta_aphan_both_image = (mean_mtd_beta_aphan_places_image + mean_mtd_beta_aphan_faces_image)/2;

   

%% ------- Subtracting means ---------


    mtd_dif_faces = ((mean_mtd_beta_aphan_faces_imagine - mean_mtd_beta_aphan_faces_image) - (mean_mtd_beta_ctrl_faces_imagine - mean_mtd_beta_ctrl_faces_image));
    
    mtd_delta_faces = mean_mtd_beta_aphan_faces_imagine - mean_mtd_beta_ctrl_faces_imagine;
    
    mtd_delta_faces_image = mean_mtd_beta_aphan_faces_image - mean_mtd_beta_ctrl_faces_image;
    
    mtd_dif_places = ((mean_mtd_beta_aphan_places_imagine - mean_mtd_beta_aphan_places_image) - (mean_mtd_beta_ctrl_places_imagine - mean_mtd_beta_ctrl_places_image));
    
    mtd_dif_both = ((mean_mtd_beta_aphan_both_imagine - mean_mtd_beta_aphan_both_image) - (mean_mtd_beta_ctrl_both_imagine - mean_mtd_beta_ctrl_both_image));
    
    mtd_delta_places = mean_mtd_beta_aphan_places_imagine - mean_mtd_beta_ctrl_places_imagine;
    
    mtd_delta_both_image = mean_mtd_beta_aphan_both(:,:,2) - mean_mtd_beta_ctrl_both_image; 



    
% --------------------------------- 
%% PERMUTATION TESTING 
% ---------------------------------

%% Places perception

%creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_places_aphan,[],4,23); %Flattens regions
    test2 = reshape(mtd_beta_mtx_places_ctrl,[],4,20);
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,2,:));
        data2 = squeeze(test2_lower(ii,2,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_places_image(1,:) = sig;
    perm_test_places_image(2,:) = pval;


    save ('perm_test_places_image.mat', 'perm_test_places_image');

%% Places imagine

%creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_places_aphan,[],4,23); %Flattens regions
    test2 = reshape(mtd_beta_mtx_places_ctrl,[],4,20);
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,4,:));
        data2 = squeeze(test2_lower(ii,4,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_places_imagine(1,:) = sig;
    perm_test_places_imagine(2,:) = pval;
    
    save ('perm_test_places_imagine.mat', 'perm_test_places_imagine');


%% Faces perception

%creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_faces_aphan,[],4,23); %Flattens regions
    test2 = reshape(mtd_beta_mtx_faces_ctrl,[],4,20);
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,2,:));
        data2 = squeeze(test2_lower(ii,2,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_faces_image(1,:) = sig;
    perm_test_faces_image(2,:) = pval;
    
    save ('perm_test_faces_image.mat', 'perm_test_faces_image');



%% Faces imagine

%creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_faces_aphan,[],4,23); %Flattens regions
    test2 = reshape(mtd_beta_mtx_faces_ctrl,[],4,20);
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,4,:));
        data2 = squeeze(test2_lower(ii,4,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_faces_imagine(1,:) = sig;
    perm_test_faces_imagine(2,:) = pval;
    
    save ('perm_test_faces_imagine.mat', 'perm_test_faces_imagine');

    
    
%% Both imagine

    mtd_beta_mtx_both_aphan = cat(4,mtd_beta_mtx_faces_aphan,mtd_beta_mtx_places_aphan);
    mtd_beta_mtx_both_ctrl = cat(4,mtd_beta_mtx_faces_ctrl,mtd_beta_mtx_places_ctrl);
    
    %creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_both_aphan,[],4,size(mtd_beta_mtx_both_aphan,4)); %Flattens regions
    test2 = reshape(mtd_beta_mtx_both_ctrl,[],4,size(mtd_beta_mtx_both_ctrl,4));
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,4,:));
        data2 = squeeze(test2_lower(ii,4,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_both_imagine(1,:) = sig;
    perm_test_both_imagine(2,:) = pval;
    
    save ('perm_test_both_imagine.mat', 'perm_test_both_imagine');
    

%% Both image


    %creates matrix of lower triangle of full ROIxROI matrix
    test = reshape(mtd_beta_mtx_both_aphan,[],4,size(mtd_beta_mtx_both_aphan,4)); %Flattens regions
    test2 = reshape(mtd_beta_mtx_both_ctrl,[],4,size(mtd_beta_mtx_both_ctrl,4));
    test_lower = test(template,:,:); %lower triangle
    test2_lower = test2(template,:,:);
    
    
    for ii = 1:size(test_lower,1)
        data1 = squeeze(test_lower(ii,2,:));
        data2 = squeeze(test2_lower(ii,2,:));
        [sig(1,ii),pval(1,ii)] = perm_code(data1,data2,5000);
    end
    
    perm_test_both_imagine(1,:) = sig;
    perm_test_both_imagine(2,:) = pval;
    
    save ('perm_test_both_image.mat', 'perm_test_both_image');


% ------------------------------------------------------------------
%% Calculating mean across Schaeffer 7-networks
% ------------------------------------------------------------------

% Load network labels for each cortical ROI 
    load('schaef_id14.mat'); 
    % 1 = visual
    % 2 = motor
    % 3 = dorsal attention
    % 4 = ventral attention
    % 5 = frontotemporal
    % 6 = frontoparietal
    % 7 = default

%% Aphan - faces imagine    
    
    load('perm_test_faces_imagine.mat')
    sig_perm_faces_imagine = squareform(perm_test_faces_imagine(1,:));     
    aphan_face_sig_imagine = zeros(502);
    
    aphan_face_sig_imagine(sig_perm_faces_imagine==1) = mtd_delta_faces(sig_perm_faces_imagine == 1); % Isolate ROIs with statistically significant differences

    
% mean per network
    mat_net_imagine = zeros(14);
    for ii = 1:14
        for jj = 1:14
            mat_net_imagine(ii,jj) = mean(mean(aphan_face_sig_imagine(schaef_id14 == ii,schaef_id14 == jj)));
        end
    end
    
    
%% Aphan - faces image      
    
    load('perm_test_faces_image.mat')
    sig_perm_faces_image = squareform(perm_test_faces_image(1,:));
    
    aphan_face_sig_image = zeros(502);
    
    aphan_face_sig_image(sig_perm_faces_image==1) = mtd_delta_faces(sig_perm_faces_image==1);
    
    
 % mean per network
    mat_net_image = zeros(14);
    for ii = 1:14
        for jj = 1:14
            mat_net_image(ii,jj) = mean(mean(aphan_face_sig_image(schaef_id14==ii,schaef_id14==jj)));
        end
    end
    
    
    
%%---- Temporal2rsn for Aphantasics places ----
    
    load('temp_lobe_mask.mat')
    
    new_id = schaef_id14+2; %Create variable new_id with leaving two rows for temporal lobe
    
    new_id(temp_lobe_mask==1) = 1; %Apply temp_lobe_mask
    
    temp = temp_lobe_mask(1:200);%Left temporal lobe
    
    temp(201:400) = temp_lobe_mask(201:400).*2; %Right temporal lobe
    
    new_id(temp_lobe_mask>0) = temp(temp_lobe_mask>0);
    
    unique(new_id);
    
    [order_new,sort_new] = sort(new_id); 
    
    
    aphan_places_delta = (mean_mtd_beta_aphan_places_imagine - mean_mtd_beta_aphan_places_image);
    aphan_places_delta_temp2rsn = aphan_places_delta(order_new<3, order_new>2);
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state network
    lt_aphan_places_delta_temp2rsn = zeros(14,1);
    rt_aphan_places_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_aphan_places_delta_temp2rsn(xx,1) = mean(mean(aphan_places_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_aphan_places_delta_temp2rsn(xx,1) = mean(mean(aphan_places_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
%%---- Temporal2rsn for Aphantasics faces ----
    
    aphan_faces_delta = (mean_mtd_beta_aphan_faces_imagine - mean_mtd_beta_aphan_faces_image);
    aphan_faces_delta_temp2rsn = aphan_faces_delta(order_new<3, order_new>2);
    
    

%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_aphan_faces_delta_temp2rsn = zeros(14,1);
    rt_aphan_faces_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_aphan_faces_delta_temp2rsn(xx,1) = mean(mean(aphan_faces_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_aphan_faces_delta_temp2rsn(xx,1) = mean(mean(aphan_faces_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
%%---- Temporal2rsn for Aphantasics both ----
    
    aphan_both_delta = (mean_mtd_beta_aphan_both_imagine - mean_mtd_beta_aphan_both_image);
    aphan_both_delta_temp2rsn = aphan_both_delta(order_new<3, order_new>2);
    
    
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_aphan_both_delta_temp2rsn = zeros(14,1);
    rt_aphan_both_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_aphan_both_delta_temp2rsn(xx,1) = mean(mean(aphan_both_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_aphan_both_delta_temp2rsn(xx,1) = mean(mean(aphan_both_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
%%---- Temporal2rsn for ctrl faces ----
    
    ctrl_faces_delta = (mean_mtd_beta_ctrl_faces_imagine - mean_mtd_beta_ctrl_faces_image);
    ctrl_faces_delta_temp2rsn = ctrl_faces_delta(order_new<3, order_new>2);
    
    
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_ctrl_faces_delta_temp2rsn = zeros(14,1);
    rt_ctrl_faces_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_ctrl_faces_delta_temp2rsn(xx,1) = mean(mean(ctrl_faces_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_ctrl_faces_delta_temp2rsn(xx,1) = mean(mean(ctrl_faces_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
    
%%---- Temporal2rsn for ctrl places ----
    
    ctrl_places_delta = (mean_mtd_beta_ctrl_places_imagine - mean_mtd_beta_ctrl_places_image);
    ctrl_places_delta_temp2rsn = ctrl_places_delta(order_new<3, order_new>2);
    
    
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_ctrl_places_delta_temp2rsn = zeros(14,1);
    rt_ctrl_places_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_ctrl_places_delta_temp2rsn(xx,1) = mean(mean(ctrl_places_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_ctrl_places_delta_temp2rsn(xx,1) = mean(mean(ctrl_places_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
%%---- Temporal2rsn for ctrl both ----
    
    ctrl_both_delta = (mean_mtd_beta_ctrl_both_imagine - mean_mtd_beta_ctrl_both_image);
    ctrl_both_delta_temp2rsn = ctrl_both_delta(order_new<3, order_new>2);
    
    
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_ctrl_both_delta_temp2rsn = zeros(14,1);
    rt_ctrl_both_delta_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_ctrl_both_delta_temp2rsn(xx,1) = mean(mean(ctrl_places_delta_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_ctrl_both_delta_temp2rsn(xx,1) = mean(mean(ctrl_places_delta_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
    
    
    
    
    
    
    
% -------------------------- -------------------------- -------------------------- --------------------------
%% Aphan vs Ctrl perception
%  -------------------------- -------------------------- -------------------------- --------------------------
    
    mtd_delta_image_all = (mean_mtd_beta_aphan_both_image - mean_mtd_beta_ctrl_both_image);
    sig = squareform(perm_test_both_image(1,:));
    mtd_delta_image_all = mtd_delta_image_all.*sig;
    
    mtd_delta_image_all_temp2rsn = mtd_delta_image_all(order_new<3, order_new>2);
    
    

%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_mtd_delta_image_all_temp2rsn = zeros(14,1);
    rt_mtd_delta_image_all_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_mtd_delta_image_all_temp2rsn(xx,1) = mean(mean(mtd_delta_image_all_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_mtd_delta_image_all_temp2rsn(xx,1) = mean(mean(mtd_delta_image_all_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
    
% -------------------------- -------------------------- -------------------------- --------------------------
%% Aphan vs Ctrl Imagery
% -------------------------- -------------------------- -------------------------- --------------------------
    
    mtd_delta_imagine_all = (mean_mtd_beta_aphan_both_imagine - mean_mtd_beta_ctrl_both_imagine);
    sig = squareform(perm_test_both_imagine(1,:));
    mtd_delta_imagine_all = mtd_delta_imagine_all.*sig;
    
    mtd_delta_imagine_all_temp2rsn = mtd_delta_imagine_all(order_new<3, order_new>2);
    
    
    
%calculate left and right side temporal lobe avg. MTD for each
%Resting-state networks
    lt_mtd_delta_imagine_all_temp2rsn = zeros(14,1);
    rt_mtd_delta_imagine_all_temp2rsn = zeros(14,1);
    
    
    
    for xx = 1:14
    net_id = xx+2;
    lt_mtd_delta_imagine_all_temp2rsn(xx,1) = mean(mean(mtd_delta_imagine_all_temp2rsn(1:33,find(order_new_notemp==net_id))));
    
    rt_mtd_delta_imagine_all_temp2rsn(xx,1) = mean(mean(mtd_delta_imagine_all_temp2rsn(34:66,find(order_new_notemp==net_id))));
    end
    
