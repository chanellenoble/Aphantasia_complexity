%onset defines the start time of the event in secs
%duration defines the duration for the event in secs
%tr is the number of repeats
%% NB:
% Requires  @(#)spm_hrf.m	2.8 Karl Friston 02/07/31 in order to run.
% SPM2 - Statistical Parametric Mapping. Wellcome Centre for Human Neuroimaging, UCL. Available at: https://www.fil.ion.ucl.ac.uk/spm/
%
% Citation: Friston, K.J., et al. (2003). Statistical parametric mapping: The analysis of functional brain images. Academic Press.
 

function dsmtx = dsmtx_hrf(onsets,tr)

    hrf = spm_hrf(tr);
    
    nOnsets = size(onsets,1);
    
    onsets_tr = onsets./(tr);
    
    time_temp = zeros(round(max(onsets_tr)),1);
    
    for n = 1:nOnsets
        time_temp(round(onsets_tr(n)),1) = 1;
    end

    dsmtx = conv(hrf,time_temp);

end
