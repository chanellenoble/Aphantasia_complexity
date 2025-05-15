% Plot vector of Schaeffer lodings onto FS surface plot with Oldam's ROI
% boundaries
% 
% NOTE: this function requires the 32k FS average .gii files and the matching
% Schaefer parcellation .gii files to work.
% 




    %  -------- Set Working Directory and Add paths --------   

    % desktop
%         cd('Z:\PRJ-shine_hpc\Aphantasia_fMRI');
%         addpath('Z:\PRJ-shine_hpc\Aphantasia_fMRI\Scripts');
%         addpath('Z:\PRJ-shine_hpc\Aphantasia_fMRI\Scripts\StuartJO-plotSurfaceROIBoundary');
%         addpath('Z:\PRJ-shine_hpc\Aphantasia_fMRI\colour maps');
%         addpath('Z:\PRJ-shine_hpc\Aphantasia_fMRI\Scripts\surf_schaef');
%         addpath('Z:\PRJ-shine_hpc\Aphantasia_fMRI\Scripts\surf_schaef\gifti-main');
% % 
% %      % Mac
        cd('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI');
        addpath('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/Scripts');
        addpath('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/Scripts/StuartJO-plotSurfaceROIBoundary');
        addpath('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/colour maps');
        addpath('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/Scripts/surf_schaef');
%         
% 
% 



    
    % -------- Load Files --------
        
        % Gifti doesn't work unless you have in in directory
        cd('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/Scripts/surf_schaef');

        % Schaefer Parcellation
        left = gifti('schaef_left.func.gii');
        
        % Average surface
        struc_L = gifti('Conte69.L.midthickness.32k_fs_LR.surf.gii');

        cd('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI/timeseries');

        load('sub-01_run-1_timeseries.mat');
        
        cd('/Volumes/PRJ-shine_hpc/Aphantasia_fMRI');

        load('temp_lobe_mask.mat');

    
    % Set some Variables
    
    surface.vertices = struc_L.vertices;
    
    surface.faces = struc_L.faces;
    
    vertex_id = left.cdata;

    bold_bin = ts'>0; % Creates a region x time matrix with logical values s.t. 1 iff BOLD > 1

    data = bold_bin(1:200,:);
    data(find(temp_lobe_mask(1:200) == 0)) = 0;



    %Set regions outside of temp lobe mask to 0 in vertex_id

    cmap = [1 1 1; 0.1289 0.9733 0.3822];

    vertex_id2 = vertex_id;

    vec = length(data);
    
    for ii = 1:vec;
        if temp_lobe_mask(ii) == 0;
            vertex_id2(find(vertex_id2 == ii)) = 0;
        else
            continue
        end
    end
    
    vertex_id = vertex_id2;
    
    new_ROI_index = find(vertex_id2 ~= 0);
    new_ROI_index = vertex_id2(new_ROI_index);
    new_ROI_index = unique(new_ROI_index);

    cmap = [1 1 1; 0.1289 0.9733 0.3822];

% cmap = [0.2,0.4,1.0; 1.0,0.9,0.1];

% Sequence of VTL plots to explain methods

for ii = 10:12; % Want to run 1:12, but this gets a bit too much so better to do in batches.

    data = bold_bin(1:200,ii);
    data(find(temp_lobe_mask(1:200) == 0)) = 0;

    data2 = data(new_ROI_index);
    clear data
    data = double(data2);;
    data(:,2) = new_ROI_index(:,1);

climits = [min(data(:,1)) max(data(:,1))];

% Ventral
figure('Position',[-1022 751 220 326])
p = plotSurfaceROIBoundary(surface,vertex_id,data,'midpoint',cmap,1,climits);
view([-180 -90])
camlight(80,-10);
camlight(-80,-10);
axis off
axis image

end
    

    
