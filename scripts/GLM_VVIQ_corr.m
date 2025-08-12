cd 'I:\Aphantasia_fMRI\Scripts\Complexity'

load vviq.mat
load idx_ach.mat

cd 'I:\Aphantasia_fMRI'

load 'All_glm_fc'



% Create matrices with GLM betas for percpetion and imagery aligned to
% idx_ach

all_glm_faces_perc = zeros(size(idx_ach,2),502);


for ii = 1:size(idx_ach, 2)
    x = idx_ach(ii);

    % Check if the field exists and has at least 2 rows
    if ~isfield(All_glm_fc(x), 'glm_faces') || isempty(All_glm_fc(x).glm_faces)
        continue
    elseif size(All_glm_fc(x).glm_faces, 1) < 2
        continue
    else
        all_glm_faces_perc(ii,:) = All_glm_fc(x).glm_faces(2,:);
    end
end


all_glm_places_perc = zeros(size(idx_ach,2),502);

for ii = 1:size(idx_ach, 2)
    x = idx_ach(ii);

    % Check if the field exists and has at least 2 rows
    if ~isfield(All_glm_fc(x), 'glm_places') || isempty(All_glm_fc(x).glm_places)
        continue
    elseif size(All_glm_fc(x).glm_places, 1) < 2
        continue
    else
        all_glm_places_perc(ii,:) = All_glm_fc(x).glm_places(2,:);
    end
end

all_glm_perc = vertcat(all_glm_places_perc,all_glm_faces_perc);


all_glm_faces_imagery = zeros(size(idx_ach,2),502);

for ii = 1:size(idx_ach, 2)
    x = idx_ach(ii);

    % Check if the field exists and has at least 2 rows
    if ~isfield(All_glm_fc(x), 'glm_faces') || isempty(All_glm_fc(x).glm_faces)
        continue
    elseif size(All_glm_fc(x).glm_faces, 1) < 2
        continue
    else
        all_glm_faces_imagery(ii,:) = All_glm_fc(x).glm_faces(4,:);
    end
end


all_glm_places_imagery = zeros(size(idx_ach,2),502);

for ii = 1:size(idx_ach, 2)
    x = idx_ach(ii);

    % Check if the field exists and has at least 2 rows
    if ~isfield(All_glm_fc(x), 'glm_places') || isempty(All_glm_fc(x).glm_places)
        continue
    elseif size(All_glm_fc(x).glm_places, 1) < 2
        continue
    else
        all_glm_places_imagery(ii,:) = All_glm_fc(x).glm_places(4,:);
    end
end

all_glm_imagery = vertcat(all_glm_places_imagery,all_glm_faces_imagery);




%% Calculate correlation between perception and VVIQ

% Spearman's Correlation
for ii = 1:size(all_glm_perc,2)
    [rho_glm_vviq_perc(ii), ~] = corr(vviq, all_glm_perc(:,ii), 'Type', 'Spearman');
end


% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    for ii = 1:size(all_glm_perc,2)
        [rho(ii),~] = corr(vviq_rand, all_glm_perc(:,ii), 'Type', 'Spearman');
    end    
null_rho(:,x) = rho(:);
end

for ii = 1:size(all_glm_perc,2);
   pval_glm_vviq_perc(ii) = sum(abs(null_rho(ii,:)) >= abs(rho_glm_vviq_perc(ii))) / 5000;;
end

%% Calculate correlation between imagery and VVIQ

% Spearman's Correlation
for ii = 1:size(all_glm_imagery,2)
    [rho_glm_vviq_imagery(ii), ~] = corr(vviq, all_glm_imagery(:,ii), 'Type', 'Spearman');
end

% Permutation test

for x = 1:5000;
    vviq_rand = vviq(randperm(size(vviq,1)));
    for ii = 1:size(all_glm_imagery,2)
        [rho(ii),~] = corr(vviq_rand, all_glm_imagery(:,ii), 'Type', 'Spearman');
    end    
null_rho(:,x) = rho(:);
end

for ii = 1:size(all_glm_imagery,2);
   pval_glm_vviq_imagery(ii) = sum(abs(null_rho(ii,:)) >= abs(rho_glm_vviq_imagery(ii))) / 5000;;
end



%% Deltas

delta_glm = nanmean(all_glm_imagery(:,:),1) - nanmean(all_glm_perc(:,:),1);

for ii = 1:size(all_glm_perc,2);
    [sig_delta_glm(ii),pval_delta_glm(ii)] = perm_code(all_glm_perc(:,ii),all_glm_imagery(:,ii),5000);
end
