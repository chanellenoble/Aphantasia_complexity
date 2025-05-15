data1 = [beta_comp_places(3,1:23), beta_comp_faces(3,1:22)];
data1(all(~data1, 1)) = [];
data1 = data1';

data2 = [beta_comp_places(3,24:end), beta_comp_faces(3,23:end)];
data2(:,all(~data2, 1)) = [];
data2 = data2';

% Box chart
data = [data1; data2];
group = [repelem(1.25,length(data1))'; repelem(1.75,length(data2))'];
cgroup = group;

figure;
b = boxchart(group,data,'GroupByColor',cgroup,'MarkerStyle','none');
b(1).BoxFaceColor = [ 192/225,0,0];
b(2).BoxFaceColor = [57/255,83/225,165/225];


hold on
RGB_color = [ 192/225,0,0; 57/255,83/225,165/225];
RGB_temp = [];
RGB_temp(:,1) = repelem(RGB_color(1,:),45);
RGB_temp = reshape(RGB_temp,45,3);

RGB_temp2 = [];
RGB_temp2(:,1) = repelem(RGB_color(2,:),34);
RGB_temp2 = reshape(RGB_temp2,34,3); 

RGB_color2 = [RGB_temp; RGB_temp2];

group2 = [repelem(1,length(data1))'; repelem(2,length(data2))'];
scatter(group2,data,20,RGB_color2,'filled','jitter','on','jitterAmount',0.1);

ylim([-3 3]);
set(gca,'FontSize',16,'FontName','Calibri','linew',1.5, ...
    'XTick',[1 2],'XTickLabel',[],'box','off');
set(findobj(gca,'type','line'),'linew',2);







