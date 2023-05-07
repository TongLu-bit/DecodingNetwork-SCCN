# SCCN
This is a repository that contains codes for the **S**patially **C**onstrained and **C**onnected **N**etworks (SCCN) model. The main focous is to detect covariate-related sub-regions **between** two large brain regions, where codes are stored in the "**SCCN_between**" folder.

Specifically, it includes:
1. Matlab codes ("*SCCN_alg.m*" and "*Reshuffle_W.m*") for implementing the SCCN model to extract densely altered sub-area pairs from a pair of regions using brainimaging data (e.g., fMRI)

2. A live-script (which is matlab version of Rmarkdown) for SCCN anlaysis on a synthetic data ("*Toy_example.m*")

Note that the optimization is non-convex and may depend on the intitialization of the algorithms.

Here is an overview of SCCN (between-regions):

\
![Image](/SCCN_between/SCCN_pipeline.png)


In the "**SCCN_within**" folder, we provide codes for detecting covariate-related sub-regions **within** a large region of interest. 
