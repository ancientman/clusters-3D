This pipeline analyzes clusters of fluorophore labeled protein in the nuclei of _Drosophila_ embryos. The input files are either in CZI format of in mutistack TIFF format. For CZI format files the pipeline uses bioformats, so the appropriate files need to be downloaded first. For positional data, the pipeline also expects an input image stack with the front and the rear end of the embryos taken at the midsaggital plane of the embryo, and calculates the absolute position of the embryo extremities using those. 

The main function is coreFunctions/starterFun3.m for TIFF stacks, and coreFunctions/starterFunCZI2.m for CZI files. The functions take input folder which is currently /Analysis_RM. 

The pipeline also takes a text file for arguments which is then parsed. All arguments must be changed in this text file. The variable names in this file and throughout the pipeline are descriptive.
Several functions although redundant are retained, as some pertain to certain imaging conditions not used in the final figures of the paper, but might be useful for a broader use. Function variables are almost always explicitly passed, to make them easily used as a standalone function. 
Lastly, use of Matlab r2023a or later is recommended to run the files.
