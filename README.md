# The multiplex dynamic functional connectivity (mdFCN) MATLAB Package
The multiplex dynamic functional connectivity network (mdFCN) analysis pipeline for objectively uncovering latent network structures from cortical activity. The material in this repository is provided to supplement the following paper:
Al?Sa'd, M., Vanhatalo, S. and Tokariev, A., ‚ÄúA Multiplex dynamic networks in the newborn brain disclose latent links with neurobehavioral phenotypes‚Äù, *Human Brain Mapping*, (2024), https://doi.org/10.1002/hbm.26610.
![plot](./ Results/ hbm26610-fig-0001-m.jpg)

The MATLAB scripts, functions, and data listed in this repository are used to produce results, and supporting figures illustrated in the paper.

## Demo Scripts:
The developed mdFCN package contains the following demo scripts within its directory:
### Demo_1_mdFC_pipeline.m
-   This demo script produces generates the study overview in Fig. 1. Note that the neonateís raw EEG and neurocognitive scores are not supplied, and the code uses random numbers as example.
### Demo_2_correlation_analysis.m
-   This demo script generates the correlation analysis results in Fig. 4. Note that the neonateís raw EEG and neurocognitive scores are not supplied, and the code uses random numbers as example. ### Demo_3_multiplexity_analysis.m
-   This demo script generates the multiplexity analysis results in Fig. 5.
### Demo_4_pipeline_verification.m
-   This demo script generates the mdFCN verification results in Fig. 6.

## Main Scripts:
The developed mdFCN package contains the following main scripts within its directory:
### Main_1_FC.m
-   This main script calculates the static and dynamic functional connectivity (FC) of each neonate's EEG. Note that the neonateís raw EEG and neurocognitive scores are not supplied.
-   The FCs of the HC group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/ERoEGI1JET1JriyvItu7sDcBh27eTqihky6cY1BSY-PHrQ?e=Wz7VSR
-   The FCs of the AED group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EfRENcju9JJLoxtC2hAUGZgBSQisD_WXuGcDV6TGwoqXog?e=mpN7JP
### Main_2_scores.m
-   This main script generates masks for the missing scores in each group and each sleep state. Note that the neonateís raw EEG and neurocognitive scores are not supplied.
### Main_3_extraction.m
-   This main script extracts dynamic FC latent networks by NMF and selects the best model order.
-   The extracted networks of the HC group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/Ebvj0IVn-T9LpqhKZITmhgMBcbUeCiCyrbBnfCdNT7hEew?e=WhuMmG
-   The extracted networks of the AED group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EVNFKvv3ReFBoXlRzluCa0sBBroVECsYYfgXp0Fa5B8ljg?e=5gROTH
### Main_4_decomposition.m
-   This main script decomposes the latent networks by CPD and selects the best model order.
-   The decomposed networks of the HC group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EdGmhXzNFu5Ki1B7HwD6pDIB6pWJmikLaL6fxgy-wkoYew?e=oghBjF
-   The decomposed networks of the AED group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EU-vQTmDt39NrjGiEUb9dFAB0HBZMf1Wv1PJnm8XwBKwJg?e=JpeXZx
### Main_5_reconstruction.m
-   This main script reconstructs the decomposed networks and generates null distributions for all other reconstruction configurations. Note that the neonateís raw EEG and neurocognitive scores are not supplied.
-   The reconstructed networks of the HC group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EZiSNevWKkNMmmKeDIY0ZP8BkLcsyqL-T0xk7jYoaoeTrQ?e=sDgIki
-   The reconstructed networks of the AED group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/Eagqz-gZb1RPjHlC0zkgsucBA8a2SCw7yuW1vI6eBumRCg?e=793hGn
### Main_6_system_verification.m
-   This main script verifies the mdFCN analysis pipeline. Note that the neonateís raw EEG and neurocognitive scores are not supplied.
-   The verified results for the HC group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EdK61BuuZu5IuE-isHp6CoIBMLHjTWWHWGtNQm54GFSj4Q?e=viBYhX
-   The verified results for the AED group can be downloaded from: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EbIlt5pdq_dInE1ucgMtF4sBu5g6toR07_a_3ZYAyv-LGA?e=PF6kYW

## Functions:
The developed mdFCN package is comprised of the following MATLAB functions that are in specific folders within the Functions directory:
### Preprocessing
-   *compute_parcel_signals.m*: It runs the preprocessing stage in Fig. 1.
-   *iir_filters.m*: It generates Bandpass filters via cascaded zero-phase low-pass and high-pass infinite impulse response filters.
-   *parcel_seg.m*: It segments the parcel signals.
### Functional Connectivity
-   *phase_lag_index.m*: It computes the phase-lag index (PLI) measures of connectivity.
-   *static_pli_measures.m*: It computes the static PLI measures of connectivity for different frequency bands.
-   *dynamic_pli_measures.m*: It computes the dynamic PLI measures of connectivity for different frequency bands.
### NMF
-   You need to request and download *The non-negative matrix factorization toolbox for biological data mining* from: https://sites.google.com/site/nmftool/home/source-code
-   *nmfnnls_mod.m*: It is a modified version of nmfnnls.m in \nmfv1_4.
### CPD
-   You need to request and download *Tensorlab | A Matlab package for tensor computations* from: https://www.tensorlab.net/
-   *cpd_network_time_avg.m*: It averages the mdFCN temporal variations.
-   *cpd_positive_als.m*: It is a modified version of cpd_als.m in \tensorlab_2016-03-28.
### Multiplexity
-   *multiplex_measures.m*: It calculates various measures of structural multiplexity.
### Correlation
-   *corr_density.m*: It calculates the density of significant positive/negative correlations.
-   *make_FDR.m*: It applies multiple comparisons correction using the BenjaminiñHochberg procedure.
-   *partialcorr_score_tensor.m*: It calculates correlation using Spearman's rank coefficient with the infants' conceptional age acting as a covariate. 
-   *tensor_correct_median.m*: It standardizes the connectivity patterns of every subject by removing their median and setting their median absolute deviation to one.
### Plotting
-   *plot3d_atlas.m*: It plots the brain atlas and the cortical regions in 3D.
-   *plot_mult_signals.m*: It plots multichannel signals at different frequencies.
-   *plot_signals.m*: It plots multichannel signals.
-   *violinPlot.m*: It is a modified code for generating violin plots for visualizing multiple distributions. The original code is developed by Jonas Dorn (jonas.dorn@gmail.com).
### Other
-   *fold_upper_matrix.m*: It matricizes the vectorized connectivity values into a square matrix.
-   *unfold_upper_matrix.m*: It vectorizes the upper triangle of a square matrix.
-   *fold_connectivity.m*: It folds the connectivity values back to their square form with a supplied fidelity matrix.
-   *order_selection.m*: It selects a model order using an entropy-based technique to avoid under and over-fitting problems.
## Head Model:
*The developed mdFCN package contains various atlas and head models used for the cortical source reconstruction task and for plotting*.

## Results:
The developed mdFCN package contains the following data that are in specific folders within Results\HC or Results\AED directories:
### FC: This folder holds the output of the main script *Main_1_FC.m*; the FC for each neonate in the HC and AED groups. Alternatively, the FCs can be downloaded from:
-   The FCs of the HC group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/ERoEGI1JET1JriyvItu7sDcBh27eTqihky6cY1BSY-PHrQ?e=Wz7VSR
-   The FCs of the AED group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EfRENcju9JJLoxtC2hAUGZgBSQisD_WXuGcDV6TGwoqXog?e=mpN7JP
### NMF: This folder contains the output of the main script *Main_3_extraction.m*; the extracted latent networks for short and long-term scores, for the five frequency bands, and for the different model orders. Alternatively, the latent networks can be downloaded from:
-   The extracted networks of the HC group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/Ebvj0IVn-T9LpqhKZITmhgMBcbUeCiCyrbBnfCdNT7hEew?e=WhuMmG
-   The extracted networks of the AED group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EVNFKvv3ReFBoXlRzluCa0sBBroVECsYYfgXp0Fa5B8ljg?e=5gROTH
### CPD: This folder contains the output of the main script *Main_4_decomposition.m*; the decomposed latent networks for short and long-term scores and for the different model orders. Alternatively, the decomposed networks can be downloaded from:
-   The decomposed networks of the HC group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EdGmhXzNFu5Ki1B7HwD6pDIB6pWJmikLaL6fxgy-wkoYew?e=oghBjF
-   The decomposed networks of the AED group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EU-vQTmDt39NrjGiEUb9dFAB0HBZMf1Wv1PJnm8XwBKwJg?e=JpeXZx
### Selection: This folder holds the output of the main script *Main_5_reconstruction.m*; the reconstructed networks for short and long-term scores along with their correlations to the scores and null distributions. Alternatively, this folder can be downloaded from:
-   The reconstructed networks of the HC group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EZiSNevWKkNMmmKeDIY0ZP8BkLcsyqL-T0xk7jYoaoeTrQ?e=sDgIki
-   The reconstructed networks of the AED group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/Eagqz-gZb1RPjHlC0zkgsucBA8a2SCw7yuW1vI6eBumRCg?e=793hGn
### Verification: This folder contains the output of the main script *Main_6_system_verification.m*; the mdFCN verification results. Alternatively, this folder can be downloaded from:
-   The verified results for the HC group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EdK61BuuZu5IuE-isHp6CoIBMLHjTWWHWGtNQm54GFSj4Q?e=viBYhX
-   The verified results for the AED group: https://helsinkifi-my.sharepoint.com/:u:/g/personal/almohamm_ad_helsinki_fi/EbIlt5pdq_dInE1ucgMtF4sBu5g6toR07_a_3ZYAyv-LGA?e=PF6kYW
