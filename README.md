# Apoptosis-Genome-Resequencing-Data
6/5/2018
Erin Roberts
Dr. Gomez-Chiarri Lab, University of Rhode Island, Kingston, RI, USA. 

This repository contains code I used to process whole genome resequencing data to characterize structural variation in the GIMAP gene family of the eastern oyster, *Crassostrea virginica*.

The Major steps in the pipeline are as follows:

1. Pre-processing of Illumina WGS reads
2. Alignment to reference genome via BWA-MEM
3. Filtering of alignment
4. Detection of Structural Variants using LUMPY
5. Genotyping of Structural Variants using SVTyper
6. Identification of consensus genotypes using SURVIVOR
7. Further analysis using VCFTools

For questions regarding this work please contact Erin Roberts, University of Rhode Island, Kingston RI at erin_roberts@my.uri.edu


