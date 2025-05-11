# Aphantasia_fMRI
Code repository for the paper - 'Seeing through the Static: Reduced Imagery Vividness in Aphantasia is Associated with Impaired Temporal Lobe Signal Complexity'

1. After pre-processing using fMRIPrep's pipeline, additional denoising is performed to regress six headmotion parameters and their temporal derivatives, as well as applying a band-pass filter and z-scoring. 
2. Data is run through 'data_clean' to save into a data structure before further analysis.
