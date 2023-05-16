# SCCN
This is a repository that contains codes for the **S**patially **C**onstrained and **C**onnected **N**etworks (SCCN) model. The main focous is to detect covariate-related sub-regions **between** two large brain regions, where codes are stored in the "**SCCN_between**" folder. We also provide codes for detecting covariate-related sub-regions **within** a large region of interest in the "**SCCN_within**" folder.

Specifically, the "**SCCN_between**" folder includes:

1.1. Matlab codes ("*SCCN_alg.m*" and "*Reshuffle_W.m*") for implementing the SCCN model to extract densely altered sub-area pairs from a pair of regions using brainimaging data (e.g., fMRI)

1.2. A live-script (which is matlab version of Rmarkdown) for SCCN-between anlaysis on a synthetic data ("*Toy_example.m*")

Note that the optimization is non-convex and may depend on the intitialization of the algorithms.

***************************************************************************************************************************************************************************************************
For who are interested to detect sub-networks within a regions, see "**SCCN_within**" folder which includes:

2.1. Matlab codes ("**SCCN_within.m**"）for implementing the SCCN model to extract densely altered sub-area pairs from a pair of regions using brainimaging data (e.g., fMRI)
"SCCN_within.m" requires auxiliry functions "SICERS_skip.m", "generate_data.m","pick_case_idx.m","custom_statistic.m", all of which are included in the folder.

2.2. A live-script (which is matlab version of Rmarkdown) for SCCN-within anlaysis on a synthetic data ("*toy_example.m*")


*****************************************************************************************************************************************************************************************************
Here is an overview of SCCN (between-regions):

\
![Image](/SCCN_between/SCCN_pipeline.png)


