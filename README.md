# CYFOR
Decoding circadian gene regulation via statistical modeling and machine learning

# Overview
Welcome to CYFOR, a machine-learning framework to decode circadian gene regulation through a combination of statistical modeling, ChIP-seq peak list annotation, random forest machine-learning models, and analyses on existing experimental data!

# Code Organization
The code for this project is organized into the following files:
- CYFOR_Oscillation - Modeling gene expression activity across time and identifying circadian genes
- Overlap_Analysis - Annotating ChIP-Seq peak lists and creating training data frame by determining which genes are regulated by which transcription factors
- RF_Model - Training a random forest machine learning model to predict the effects of regulatory activity on circadian oscillation
- MYC_Disruption - Analyzing existing public experimental data to determine whether MYC overexpression disrupts circadian regulation

# Data
- Public time-course RNA-sequencing data from 64 different tissues in baboons was downloaded (paper: https://www.science.org/doi/10.1126/science.aao0318) and used to fit oscillations and train random forests machine learning models.
- All available human ChIP-Seq peak lists were downloaded from the Encyclopedia of DNA Elements (ENCODE). A randomly selected list for each ChIP-Seq list was annotated used to create the training dataset. MYC and ARNTL lists were also used to perform overlap to see if their binding sites were largely the same.
- Time-course MYC overexpression data was downloaded (paper: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010904#sec011) and used to determine whether experimental data supported that MYC overexpression disrupts circadian regulation.
