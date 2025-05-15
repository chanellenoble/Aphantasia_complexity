%Calculates the correlation between MTD scores and K-complexity

beta_comp_aphan_imagine = horzcat(kbeta_comp_places(3,1:23),kbeta_comp_faces(3,1:23))';
beta_mtd_imagine = cat(4, mtd_beta_mtx_places_aphan(:,:,4,:), mtd_beta_mtx_faces_aphan(:,:,4,:));
% beta_mtd_imagine = squeeze(beta_mtd_imagine);

for r = 1:502;
    for s = 1:502;
        bms = squeeze(beta_mtd_imagine(r,s,:));
        [rcm_k(r,s,:),p(r,s,:)] = corr(bms,beta_comp_aphan_imagine);
    end
end


 rcm2_k = fillmissing(rcm_k, 'constant', 0);
mean_rcm2_k = mean(rcm2_k,1);