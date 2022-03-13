### NOTE

Files in this repository are merged into a new repository-[Oslops](https://github.com/gangwug/Oslops). Please visity the new repository for detailed information and future updates.

## Files

RonCYCLOPSv3: the source julia code file from CYCLOPS (julia version 0.3.12)

runCYCLOPS_Eigen.jl: the code file for getting the Eigen genes

runOscope.R: the code file for selecting Eigen clusters

runCYCLOPS_Order.jl: the code file for ordering

clockList.csv: the input seed gene list for ordering

CYCLOPSout: the output folder for CYCLOPS ordering

example.ToRunCYCLOPS.csv: an example of expression data

skinBenchmarkMatrixARNTLorder.csv: the benchmark correlation matrix used for down-sampling (minor different from the 17 clock and clock-associated genes, which may be helpful for avoiding the bias of using the same 17 seed gene list during ordering)

randomSampling.R: the code file for maximizing the oscillation signal by down-sampling

## Advice about running this pipeline

Step 1: prepare the input data file (e.g., example.ToRunCYCLOPS.csv). Make sure the first column is gene name/id, other columns are samples, and make sure that there are no duplicate gene name/id in the first column. If the input expression data file does not have a functional clock, tested by [nCV](https://github.com/gangwug/nCV) using 'skinBenchmarkMatrixARNTLorder.csv' as the reference matrix, randomSampling.R may be helpful for selecting subset of samples with a functional clock. 

Step 2: run 'runCYCLOPS_Eigen.jl' to generate the eigen genes.

Step 3: run 'runOscope.R' to select clusters of cycling eigen genes with different phase.

Step 4: run 'runCYCLOPS_Order.jl' to order the samples using selected eigen gene clusters.

Step 5: Check and evaluate the output results in 'CYCLOPSout' folder. 

## License
This package is free and open source software, licensed under GPL(>= 2).

