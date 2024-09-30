# Shape-metrics-Gliomas
2D and 3D Fractal Dimension and Lacunarity - Calculation and Analysis in Gliomas

The data and codes are a part of the paper uploaded to BioRxiv titled: "Fractality and Lacunarity Measures of Glioma Subcomponents are Discriminative of IDH Status: A Quantitative Radiogenomics in Gliomas", doi: https://doi.org/10.1101/2023.12.28.573519 and "Fractal Dimension and Lacunarity Measures of Glioma Subcomponents are Discriminative of Grade of Gliomas and IDH Status" accepted in NMR in Biomedicine.


# Contents:
1. All codes used in the paper for visualization, calculation and analysis.
2. Datasheets containing all the demographic and calculated information.
3. Glioma mask images 

## Codes - File Description
### Codes related to calculation of Fractal Dimension and Lacunarity
1. Python files - FD_Calculator and Lac_Calculator (Calculates 2D fractal Dimension) - to be initially run as python files on a bash terminal to calculate fractal dimension (FD) of the three glioma subcomponents as mentioned in the manuscript of either a single subject (to test the integrity of the code) or all the subjects at once (to be run once by following the commands in the code).
2. 3D Fractal Dimension Calculator (3D Frac) - Calculated Fractal Dimension across all three planes of the MR Image.
3. 3D Lacunarity Calculator (3D Lac) - Calculates Lacunarity across all three planes of the MR Image. `3DLac_cor.ipynb` and `3DLac_sag.ipynb` are supporting notebooks required to run `3DLac.ipynb`.

> The output of 3DFrac.ipynb and 3DLac.ipynb is fractal_lac_data.csv

### Codes used for visualization, training and analysis
1. Python Notebook 1 (Glioma Final Image): To visualize the MR images of gliomas and the generated masks (Mask numbers and color codes are mentioned within the notebook).
2. Python Notebook 2 (Violin Plots): Visualization and statistical analysis of the differences between the FD and Lacunarity of the 3 combinations in IDH mutant and wildtype molecular subtypes and in groups with combination of different IDHa nd MGMT status.
3. Python Notebook 3 (Statistics): Some miscellaneous statistics concerned with the subject data and calculated information.
4. Python Notebook 4 (ML Classification): Detailed steps of the Machine Learning classification algorithm (train and test sets) and generation of relevant images and plots.
5. Python Notebook 5 (ML Classification analysis): Statistics validating the accuracy and sensitivity of the machine learning algorithms.
6. R Notebook 1 (CPH Model Hazard Ratio Calculation): Code implementing the Cox Proprtional Hazards model to find out the hazard for each group divided by the log rank statistics (mentioned in manuscipt, Methods: Statistical Analyses).
7. R Notebook 2 (Survival Analysis (KM curve)): Survival analysis using the Kaplan Meier estimator.

8.
    i. TCGA_LGG_GBM_radiomicFeatures_clinicalDetails.xlsx --> Datasheet of Clinical Information of all subjects.
    
   ii. fractal_lac_data.csv --> Calculated FD and lacunarity of the three subcomponents (including 2D and 3D calculations).

  iii. Dictionary.csv --> Dictionary explaining the abbreviations used in the fractal_lac_data.csv document.

11. Folder `All_Tumor_masks_nii_files` contains the glioma masks (glioma region in the brain) of all subjects used in this study in NIfTI file format. The glioma masks were used as the input for calculating fractal dimension and lacunarity of the glioma subcomponents.
