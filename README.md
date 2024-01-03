# Shape-metrics-Gliomas
2D Fractal Dimension and Lacunarity - Calculation and Analysis in Gliomas

This data is a part of the paper uploaded to BioRxiv titled: "Fractality and Lacunarity Measures of Glioma Subcomponents are Discriminative of IDH Status: A Quantitative Radiogenomics in Gliomas", doi: https://doi.org/10.1101/2023.12.28.573519

& in Neuro-Oncology titled the same.

# Contents:
1. All codes used in the paper for visualization, calculation and analysis.
2. Datasheets containing all the demographic and calculated information.
3. Glioma mask images 

## Codes - File Description
1. Python files - FD_Calculator and Lac_Calculator - to be initially run as python files on a bash terminal to calculate fractal dimension (FD) of the three glioma subcomponents as mentioned in the manuscript of either a single subject (to test the integrity of the code) or all the subjects at once (to be run once by following the commands in the code).
2. Python Notebook 1 (Glioma Final Image): To visualize the MR images of gliomas and the generated masks (Mask numbers and color codes are mentioned within the notebook).
3. Python Notebook 2 (Violin Plots): Visualization and statistical analysis of the differences between the FD and Lacunarity of the 3 combinations in IDH mutant and wildtype molecular subtypes and in groups with combination of different IDHa nd MGMT status.
4. Python Notebook 3 (Statistics): Some miscellaneous statistics concerned with the subject data and calculated information.
5. Python Notebook 4 (ML Classification): Detailed steps of the Machine Learning classification algorithm (train and test sets) and generation of relevant images and plots.
6. Python Notebook 5 (ML Classification analysis): Statistics validating the accuracy and sensitivity of the machine learning algorithms.
7. R Notebook 1 (CPH Model Hazard Ratio Calculation): Code implementing the Cox Proprtional Hazards model to find out the hazard for each group divided by the log rank statistics (mentioned in manuscipt, Methods: Statistical Analyses).
8. R Notebook 2 (Survival Analysis (KM curve)): Survival analysis using the Kaplam Meier estimator.
9. TCGA_LGG_GBM_radiomicFeatures_clinicalDetails.xlsx --> Excel sheet containing clinical information of all the subjects
    
   fractal_lac_data.csv --> Calculated FD and lacunarity of the three subcomponents

   Dictionary.csv --> Dictionary explaining the abbreviations used in the fractal_lac_data.csv document.
10. Folder 
