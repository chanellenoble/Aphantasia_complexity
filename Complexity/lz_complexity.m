% Calculates the Lempel-Ziv Complexity of TS.
%
%% NB:
% Requires 'calc_lz_complexity.m' by Quang Thai in order to work!
%
% Citation: 'Quang Thai (2025). calc_lz_complexity (https://www.mathworks.com/matlabcentral/fileexchange/38211-calc_lz_complexity), MATLAB Central File Exchange. Retrieved May 15, 2025.'



%%LOADING IN VARIABLES
    cd 'C:\Users\Chanelle\Sydney Uni Work\Neuroscience\Shine-Labs\Aphantasia Dataset'
    
    
    load All_dsmtx.mat
    load subj_id.mat % 1:23 -> aphan; 24:43 -> ctrl;
    load idx_ac.mat
    load sub_id_aphan.mat
    load sub_id_ctrl.mat
    load run_name.mat
    load temp_lobe_mask.mat
    addpath 'Scripts/'

% predefine some variables (note: you'll need to update these!)
    nSubs = 43;

% -------------------------------------------------------------------------------------------------------------------


%Create run_name1 s.t. it is a subset of run_name only containing subjects from ac. These are the runs that are aligned to the faces scanning sessions
    
    for qq = 1:nSubs
        run_name1(qq) = run_name(idx_ac(qq));
    end
    
    idx_check = zeros(nSubs,2); %creates a check list of the indexes for which ts is being added. Column 1 -faces; column 2 - places
    
    ac = horzcat(aphan, ctrl);

% %adding "sub-" to the front of each entry of ac
    for qq = 1:nSubs
        ac{1,qq} = sprintf('%s%s','sub-',ac{1,qq});
    end
    ac = char(ac);

%Remove hyperphantasics from All_dsmtx

    for ID = 1:nSubs
        All_dsmtx_ac(ID) = All_dsmtx(idx_ac(ID));
    end

% -------------------------------------------------------------------------------------------------------------------
%%LZ COMPLEXITY ANALYSIS
% -------------------------------------------------------------------------------------------------------------------

%FACES
    beta_comp = [];
    
    temp_comp = [];
    
    cd('C:\Users\Chanelle\Sydney Uni Work\Neuroscience\Shine-Labs\Aphantasia Dataset\timeseries')


    for ID = 1:nSubs
        file_name = sprintf('%s%s%s%d%s%s',ac(ID,:),'_','run-',run_name1(ID),'_','timeseries.mat');
    
        %Checks that file is in timeseries folder
        if isfile(file_name) == 0;
          sprintf('%s%s',file_name,' does not exist');
          % Skip to bottom of loop and continue with the loop
          continue;
        end
    
        idx_check(ID,1) = 1; %Ticks off on list to track the files that it has added
    
        bold = load(file_name); %loads TS
    
        perc_block = All_dsmtx_ac(ID).dsmtx_faces(:,1)>0.5;
        imag_block = All_dsmtx_ac(ID).dsmtx_faces(:,3)>0.5;
        bold_bin = bold.ts'>0; % Creates a region x time matrix with logical values s.t. 1 iff BOLD > 1
        for tt = 1:size(bold_bin,2)
            temp_comp(tt,ID) = calc_lz_complexity(flatten_mat(bold_bin(temp_lobe_mask==1,tt)),'primitive',0);% Column corresponds to ID
        end
        beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to ID
    
     end


    beta_comp_faces = beta_comp;
    
    temp_comp_faces = temp_comp;
    
    clear temp_comp beta_comp

%PLACES    -------------------------------------------------------------------------

    beta_comp = [];
    
    temp_comp = [];
    
    for ID = 1:nSubs
            
          if run_name1(ID)== 1;
              run_newname= 2;
          else run_name1(ID)== 2;
              run_newname= 1;
          end
        
          file_name = sprintf('%s%s%s%d%s%s',ac(ID,:),'_','run-',run_newname,'_','timeseries.mat');
        
          %Checks that file is in timeseries folder
          if isfile(file_name) == 0;
             sprintf('%s%s',file_name,' does not exist');
             % Skip to bottom of loop and continue with the loop
             continue;
          end
          
        idx_check(ID,2) = 1; %Ticks off on list to track the files that it has added
        
        bold = load(file_name); %loads the file
        
        perc_block = All_dsmtx_ac(ID).dsmtx_places(:,1)>0.5;
        imag_block = All_dsmtx_ac(ID).dsmtx_places(:,3)>0.5;
        bold_bin = bold.ts'>0;
        for tt = 1:size(bold_bin,2)
        temp_comp(tt,ID) = calc_lz_complexity(flatten_mat(bold_bin(temp_lobe_mask==1,tt)),'primitive',0);% Column corresponds to sub_id
        end
        beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to sub_id
    
        end
    
    beta_comp_places = beta_comp;
      
    temp_comp_places = temp_comp;




% -------------------------------------------------------------------------------------------------------------------
%% CREATE DIFFERENCE SCORES FOR TL COMPLEXITY COMPARING IMAGINATION AND PERCEPTION
% -------------------------------------------------------------------------------------------------------------------
%% 1:23 --> aphantasics; 24:end --> controls

    beta_comp_faces = beta_comp_faces(:,any(beta_comp_faces));
    beta_comp_places = beta_comp_places(:,any(beta_comp_places));
    
  delta_lz_faces = beta_comp_faces(3,:) - beta_comp_faces(2,:);
    delta_lz_places = beta_comp_places(3,:) - beta_comp_places(2,:);


% -------------------------------------------------------------------------------------------------------------------
%% PERMUTATION TEST
% -------------------------------------------------------------------------------------------------------------------


    iter = 5000;
    
    aphan_lz = horzcat(delta_lz_faces(1,1:23),delta_lz_places(1,1:23));
    ctrl_lz = horzcat(delta_lz_faces(1,24:end),delta_lz_places(1,24:end));
    %% 1:23 --> aphantasics; 24:end --> controls
    
    size1 = size(aphan_lz,2);
    size2 = size(ctrl_lz,2);
    
    null_delta = zeros(iter,1);
    
    data_combo = horzcat(aphan_lz,ctrl_lz);
    grp_combo(1,1:size1) = 1;
    grp_combo(1,size1+1:size1+size2) = 2;
    orig_delta = nanmean(data_combo(grp_combo==1))-nanmean(data_combo(grp_combo==2));
    
    for x = 1:iter
        rand_vec = rand(size1+size2,1);
        [~,sort_rand] = sort(rand_vec);
        grp_rand = grp_combo(sort_rand); 
        null_delta(x,1) = nanmean(data_combo(grp_rand==1))-nanmean(data_combo(grp_rand==2));
    end
    
    
        pos_thr = prctile(null_delta,97.5);
        neg_thr = prctile(null_delta,2.5);
        
        pos_sig = double(orig_delta > pos_thr);
        neg_sig = double(orig_delta < neg_thr);
        both_sig = [pos_sig,neg_sig];
        sig = max(both_sig); 
    %     sig = max(pos_sig,neg_sig);
       
        pos_pval = sum(null_delta>orig_delta)/iter;
        neg_pval = sum(null_delta<orig_delta)/iter;
    
        both_pval = [pos_pval,neg_pval];
        pval = min(both_pval);

