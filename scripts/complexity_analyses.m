%% Calculate the LZ-Complexity and Kolmogorov complexity within predefined ventral temporal lobe mask


% Load in Variables

cd 'I:\Aphantasia_fMRI'
addpath 'I:\Aphantasia_fMRI\Scripts'


load All_dsmtx.mat
load temp_mask.mat

cd 'I:\Aphantasia_fMRI\Scripts\Complexity'


load('valid_idx.mat')
load('ach.mat')
load('run_name.mat')
load('idx_ach.mat')
load('vviq.mat')



% Predefine Variables
nSubs = 58;
nAphan = 21;
nCtrl = 17;
nHyper = 20;

ach = char(ach); %convert ach to character

%% ---- LZ Complexity -----    



    % Faces

        beta_comp = [];
        
        temp_comp = [];
        
        cd('I:\Aphantasia_fMRI\timeseries')
        
        for ID = 1:nSubs
        
            file_name = sprintf('%s%s%s%d%s%s',ach(ID,:),'_','run-',run_name(ID),'_','timeseries.mat');
        
            %Checks that file is in timeseries folder
            if isfile(file_name) == 0;
              sprintf('%s%s',file_name,' does not exist');
              sprintf('%s%s',file_name,' does not exist');
              idx_check(ID,1) = 0;
            else
                idx_check(ID,1) = 1; %Ticks off on list to track the files that it has added
    
                bold = load(file_name); %loads TS
        
                perc_block = All_dsmtx(idx_ach(ID)).dsmtx_faces(:,1)>0.5;
                imag_block = All_dsmtx(idx_ach(ID)).dsmtx_faces(:,3)>0.5;
                bold_bin = bold.ts'>0;
                for tt = 1:size(bold_bin,2)
                    temp_comp(tt,ID) = calc_lz_complexity(flatten_mat(bold_bin(temp_lobe_mask==1,tt)),'primitive',0);% Column corresponds to ID
                end
                beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to ID
    
                end
        
            
        
        end
        
        beta_comp_faces = beta_comp;
        
        temp_comp_faces = temp_comp;
        
        clear temp_comp beta_comp idx_check


        % Places
    
            beta_comp = [];
            
            temp_comp = [];
            
           for ii = 1:size(run_name,2)
               if run_name(1,ii) == 1;
                  run_name1(1,ii) = 2;
               else
                   run_name1(1,ii) = 1;
               end
           end
            
            for ID = 1:nSubs
            
                file_name = sprintf('%s%s%s%d%s%s',ach(ID,:),'_','run-',run_name1(ID),'_','timeseries.mat');
            
                %Checks that file is in timeseries folder
                if isfile(file_name) == 0;
                  sprintf('%s%s',file_name,' does not exist');
                  idx_check(ID,1) = 0;
                else
                     idx_check(ID,1) = 1; %Ticks off on list to track the files that it has added
        
                bold = load(file_name); %loads TS
        
                perc_block = All_dsmtx(idx_ach(ID)).dsmtx_places(:,1)>0.5;
                imag_block = All_dsmtx(idx_ach(ID)).dsmtx_places(:,3)>0.5;
                bold_bin = bold.ts'>0;
                for tt = 1:size(bold_bin,2)
                    temp_comp(tt,ID) = calc_lz_complexity(flatten_mat(bold_bin(temp_lobe_mask==1,tt)),'primitive',0);% Column corresponds to ID
                end
                beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to ID

                end
            
               
            
            end
            
            beta_comp_places = beta_comp;
            
            temp_comp_places = temp_comp;
            
            clear temp_comp beta_comp



%% LZ - Spearman's Rho for all VVIQ and Complexity - Imagery

beta_comp_imagery = horzcat(beta_comp_places(3,:),beta_comp_faces(3,:)); % perception 2; imagery 3

[rho_lz, pval_lz] = corr(vviq, beta_comp_imagery', 'Type', 'Spearman');

% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    [rho, ~] = corr(vviq_rand, beta_comp_imagery', 'Type', 'Spearman');
    null_rho_lz(x) = rho;
end


perm_pval_lz = mean(null_rho_lz >= rho_lz);

    
%% ---- Make LZ vs VVIQ Figure Imagery ----

aphan_idx(1:2.*nSubs) = false;
aphan_idx(1:nAphan) = true;
aphan_idx(nSubs+1:nSubs+nAphan) = true;

ctrl_idx(1:2.*nSubs) = false;
ctrl_idx(nAphan+1:nAphan+nCtrl) = true;
ctrl_idx(nSubs+nAphan+1:nSubs+nAphan+nCtrl) = true;

hyper_idx(1:2.*nSubs) = false;
hyper_idx(nAphan+nCtrl+1:nSubs) = true;
hyper_idx(nSubs+nAphan+nCtrl+1:end) = true;

figure
hold on
h1 = scatter(vviq(aphan_idx), beta_comp_imagery(aphan_idx), 50, 'r', 'filled');
h2 = scatter(vviq(ctrl_idx), beta_comp_imagery(ctrl_idx), 50, 'o', 'filled');
h3 = scatter(vviq(hyper_idx), beta_comp_imagery(hyper_idx), 50, 'y', 'filled');

text(min(vviq+0.05*range(vviq)), max(beta_comp_imagery)-0.1*range(beta_comp_imagery),sprintf('Spearman\\rho = %.2f\\newlinep = %.3f', rho_lz, perm_pval_lz),'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Add the legend
legend([h1 h2 h3], {'Aphantasic', 'Control', 'Hyperphantasic'},'Location', [0.513928573663745,0.775642858845847,0.24107142410108,0.126190472784496])


title('Lempel-Ziv Beta Coefficients vs VVIQ Score');
xlabel('VVIQ Score');
ylabel('Lempel-Ziv Complexity Beta Coefficient');
hold off




%% LZ - Spearman's Rho for all VVIQ and Complexity Perception

beta_comp_perc = horzcat(beta_comp_places(2,:),beta_comp_faces(2,:)); % perception 2; perc 3

[rho_lz_perc, pval_lz_perc] = corr(vviq, beta_comp_perc', 'Type', 'Spearman');

% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    [rho, ~] = corr(vviq_rand, beta_comp_perc', 'Type', 'Spearman');
    null_rho_lz_perc(x) = rho;
end


perm_pval_lz_perc = mean(null_rho_lz_perc >= rho_lz_perc);

    
%% ---- Make LZ vs VVIQ Figure ----

aphan_idx(1:2.*nSubs) = false;
aphan_idx(1:nAphan) = true;
aphan_idx(nSubs+1:nSubs+nAphan) = true;

ctrl_idx(1:2.*nSubs) = false;
ctrl_idx(nAphan+1:nAphan+nCtrl) = true;
ctrl_idx(nSubs+nAphan+1:nSubs+nAphan+nCtrl) = true;

hyper_idx(1:2.*nSubs) = false;
hyper_idx(nAphan+nCtrl+1:nSubs) = true;
hyper_idx(nSubs+nAphan+nCtrl+1:end) = true;

figure
hold on
h1 = scatter(vviq(aphan_idx), beta_comp_perc(aphan_idx), 50, 'r', 'filled');
h2 = scatter(vviq(ctrl_idx), beta_comp_perc(ctrl_idx), 50, 'o', 'filled');
h3 = scatter(vviq(hyper_idx), beta_comp_perc(hyper_idx), 50, 'y', 'filled');

text(min(vviq+0.05*range(vviq)), max(beta_comp_perc)-0.1*range(beta_comp_perc),sprintf('Spearman\\rho = %.2f\\newlinep = %.3f', rho_lz_perc, perm_pval_lz_perc),'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Add the legend
legend([h1 h2 h3], {'Aphantasic', 'Control', 'Hyperphantasic'},'Location', [0.513928573663745,0.775642858845847,0.24107142410108,0.126190472784496])




%% ---- K Complexity -----    



    % Faces

        beta_comp = [];
        
        temp_comp = [];
        
        cd('I:\Aphantasia_fMRI\timeseries')
        
        for ID = 1:nSubs
        
            file_name = sprintf('%s%s%s%d%s%s',ach(ID,:),'_','run-',run_name(ID),'_','timeseries.mat');
        
            %Checks that file is in timeseries folder
            if isfile(file_name) == 0;
              sprintf('%s%s',file_name,' does not exist');
              sprintf('%s%s',file_name,' does not exist');
              idx_check(ID,1) = 0;
            else
                idx_check(ID,1) = 1; %Ticks off on list to track the files that it has added
    
                bold = load(file_name); %loads TS
        
                perc_block = All_dsmtx(idx_ach(ID)).dsmtx_faces(:,1)>0.5;
                imag_block = All_dsmtx(idx_ach(ID)).dsmtx_faces(:,3)>0.5;
                bold_bin = bold.ts'>0;
                for tt = 1:size(bold_bin,2)
                    temp_comp(tt,ID) = kolmogorov(flatten_mat(bold_bin(temp_lobe_mask==1,tt)));% Column corresponds to ID
                end
                beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to ID
    
                end
        
            
        
        end
        
        beta_comp_faces = beta_comp;
        
        temp_comp_faces = temp_comp;
        
        clear temp_comp beta_comp idx_check


        % Places
    
            beta_comp = [];
            
            temp_comp = [];
            
           for ii = 1:size(run_name,2)
               if run_name(1,ii) == 1;
                  run_name1(1,ii) = 2;
               else
                   run_name1(1,ii) = 1;
               end
           end
            
            for ID = 1:nSubs
            
                file_name = sprintf('%s%s%s%d%s%s',ach(ID,:),'_','run-',run_name1(ID),'_','timeseries.mat');
            
                %Checks that file is in timeseries folder
                if isfile(file_name) == 0;
                  sprintf('%s%s',file_name,' does not exist');
                  idx_check(ID,1) = 0;
                else
                     idx_check(ID,1) = 1; %Ticks off on list to track the files that it has added
        
                bold = load(file_name); %loads TS
        
                perc_block = All_dsmtx(idx_ach(ID)).dsmtx_places(:,1)>0.5;
                imag_block = All_dsmtx(idx_ach(ID)).dsmtx_places(:,3)>0.5;
                bold_bin = bold.ts'>0;
                for tt = 1:size(bold_bin,2)
                    temp_comp(tt,ID) = kolmogorov(flatten_mat(bold_bin(temp_lobe_mask==1,tt)));% Column corresponds to ID
                end
                beta_comp(:,ID) = glmfit(horzcat(perc_block(1:size(bold_bin,2),1),imag_block(1:size(bold_bin,2),1)),temp_comp(:,ID),'normal'); % Column corresponds to ID

                end
            
               
            
            end
            
            beta_comp_places = beta_comp;
            
            temp_comp_places = temp_comp;
            
            clear temp_comp beta_comp


%% K - Spearman's Rho for all VVIQ and Complexity -- Imagery

beta_comp_imagery = horzcat(beta_comp_places(3,:),beta_comp_faces(3,:));

[rho_k, pval_k] = corr(vviq, beta_comp_imagery', 'Type', 'Spearman');

% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    [rho, ~] = corr(vviq_rand, beta_comp_imagery', 'Type', 'Spearman');
    null_rho_k(x) = rho;
end


perm_pval_k = mean(null_rho_k >= rho_k);

%% ---- Make k vs VVIQ Figure Imagery ----


figure
hold on
h1 = scatter(vviq(aphan_idx), beta_comp_imagery(aphan_idx), 50, 'r', 'filled');
h2 = scatter(vviq(ctrl_idx), beta_comp_imagery(ctrl_idx), 50, 'o', 'filled');
h3 = scatter(vviq(hyper_idx), beta_comp_imagery(hyper_idx), 50, 'y', 'filled');

text(min(vviq+0.05*range(vviq)), max(beta_comp_imagery)-0.1*range(beta_comp_imagery),sprintf('Spearman\\rho = %.2f\\newlinep = %.3f', rho_k, perm_pval_k),'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Add the legend
legend([h1 h2 h3], {'Aphantasic', 'Control', 'Hyperphantasic'},'Location', [0.513928573663745,0.775642858845847,0.24107142410108,0.126190472784496])


title('Kolmogorov vs VVIQ Score');
xlabel('VVIQ Score');
ylabel('Kolmogorov Complexity');
hold off


%% k - Spearman's Rho for all VVIQ and Complexity Perception

beta_comp_perc = horzcat(beta_comp_places(2,:),beta_comp_faces(2,:)); % perception 2; perc 3

[rho_k_perc, pval_k_perc] = corr(vviq, beta_comp_perc', 'Type', 'Spearman');

% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    [rho, ~] = corr(vviq_rand, beta_comp_perc', 'Type', 'Spearman');
    null_rho_k_perc(x) = rho;
end


perm_pval_k_perc = mean(null_rho_k_perc >= rho_k_perc);

    
%% ---- Make k vs VVIQ Figure - Perception ----

aphan_idx(1:2.*nSubs) = false;
aphan_idx(1:nAphan) = true;
aphan_idx(nSubs+1:nSubs+nAphan) = true;

ctrl_idx(1:2.*nSubs) = false;
ctrl_idx(nAphan+1:nAphan+nCtrl) = true;
ctrl_idx(nSubs+nAphan+1:nSubs+nAphan+nCtrl) = true;

hyper_idx(1:2.*nSubs) = false;
hyper_idx(nAphan+nCtrl+1:nSubs) = true;
hyper_idx(nSubs+nAphan+nCtrl+1:end) = true;

figure
hold on
h1 = scatter(vviq(aphan_idx), beta_comp_perc(aphan_idx), 50, 'r', 'filled');
h2 = scatter(vviq(ctrl_idx), beta_comp_perc(ctrl_idx), 50, 'o', 'filled');
h3 = scatter(vviq(hyper_idx), beta_comp_perc(hyper_idx), 50, 'y', 'filled');

text(min(vviq+0.05*range(vviq)), max(beta_comp_perc)-0.1*range(beta_comp_perc),sprintf('Spearman\\rho = %.2f\\newlinep = %.3f', rho_k_perc, perm_pval_k_perc),'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Add the legend
legend([h1 h2 h3], {'Aphantasic', 'Control', 'Hyperphantasic'},'Location', [0.513928573663745,0.775642858845847,0.24107142410108,0.126190472784496])



