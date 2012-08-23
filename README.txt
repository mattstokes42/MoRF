MoRF_Runner.exe is compiled from four program files, each of which is described briefly below:
MoRF_Runner.cs - the overall experiment architecture, runs several replications of several models and generates summary statistics
ReliefUtils.cs - contains the actual MoRF algorithm and associated component functions
FileUtils.cs - contains tools for reading data files, writing output, etc.
CommandLineUtils.cs - contains tools for parsing command line flags

The program MoRF_Runner.exe will run the Modular Relief Framework on up to 100 data files (replications) for up to 70 genetic models. MoRF scores are calculated for each data file (replication) individually. Once the attribute ranking for each data file of a model has been found, the overall power for the model is computed as the percentage of replications which have both informative SNPs ranked above a threshold. To run the experiment, unzip the data file packages (described below) in the same directory as the executable, and run the .exe file.

INPUT: The data files required for this experiment were used by Moore et al. in the evaluation of SURF*, and may be found at http://discovery.dartmouth.edu/epistatic_data/. The data files are available in four sample sizes. For each sample size, there are 7,000 data files. There are 70 different genetic models that vary in both minor allele frequency and disease heritability, grouped into sets of 5. Each genetic model has 100 data files (replications). Each of these data files contains two informative SNPs (0 and 1), as well as 998 uninformative SNPs. The data packages should be unzipped in the same directory as the executable, so that the actual data files reside two folder levels below the executable. An example file location would be {executable directory}/200/00/00.200.0.txt.

OUTPUT: There are two types of output files. The first type is produced for each input data file (replication), and contains the raw MoRF score of each SNP. The second type of file is for summary statistics, and only one is produced for each genetic model (set of replications).

MoRF score output file - 
{Neighbor weighting function}.{p1}.{samples}.{model}.{replication}.txt
Example filename: SWRF.4.200.00.0.txt (SWRF neighbor weighting function with scale factor 4, run on 200 samples of model 00, replication 0)
Each line contains a SNP number (0-999) and final MoRF score, separated by a tab.

Power statistics summary file - 
{Neighbor weighting function}.{p1}.{samples}.{model}.{startRep}to{endRep}.power.txt
Example filename: ReliefF.10.400.05.0to99.power.txt (ReliefF neighbor weighting function with 10 neighbors, run on 400 samples of model 05, for replications 0 to 99)
Contains the power of an algorithm (fraction of datasets with SNPs 0 and 1 ranked above a threshold), taken over all the replications run for a particular model. Power is listed by percentile. That is, the first number is the power of the algorithm at the 1% threshold (i.e. the fraction of replications run that have SNPs 0 and 1 in the top-ranked 10 SNPs). The second number is the power at the 2% threshold, and so forth. This string of numbers can easily be plotted in Excel or Matlab to visualize the power of each algorithm. In Excel, load the file in tab-delimited format, and plot using standard tools. In Matlab, use the command plot(1:100, dlmread(‘power_filename.txt’)).


FLAGS: There are several command line flags, all of which are optional. Flags are marked by a hyphen and should be followed by a space and the parameter value.
-s : Samples - the number of samples (individuals) in the dataset, valid values {200, 400, 800, 1600}
-sm : Start Model - the first model to run, valid values 0-69
-em : End Model - the last model to run, valid values 0-69
-sr : Start Replication - the first replication to run, valid values 0-99
-er : End Replication - the last replication to run, valid values 0-99
-f : Neighbor Weighting function - which neighbor weighting function should be used by MoRF. See figure 2 of the paper. Valid values are {SWRF, SWRFstar, SURF, SURFstar, ReliefF, ReliefFstar, ramp}
-p1 : Relief parameter 1 - some neighbor weighting functions require a parameter, entered via this flag. For ReliefF, the parameter represents the number of nearest hits/misses to use (default value 10). For SWRF, the parameter represents the width scaling of the sigmoid function (default value 4) - larger numbers make the sigmoid narrower, while smaller numbers make it wider.
-p2 : Unused parameter, allows expansion of code to include weighting functions that require more than one user parameter.
-p3 : Unused parameter, see p2

Code written by Matthew E. Stokes
University of Pittsburgh Biomedical Informatics Program
University of Pittsburgh Intelligent Systems Program
2012