# Aphantasia_fMRI
Code repository for the paper - 'Seeing through the Static: Reduced Imagery Vividness in Aphantasia is Associated with Impaired Temporal Lobe Signal Complexity'

1. After pre-processing using fMRIPrep's pipeline, additional denoising is performed to regress six headmotion parameters and their temporal derivatives, as well as applying a band-pass filter and z-scoring. 
2. use 'onset_times' to create a data structure - onset times from 'Timing.
3. Onset times are convolved by HRF using 'dsmtx.m' and 'dsmtx_hrf_epoch'
4. 
