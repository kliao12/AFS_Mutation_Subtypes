Repository containing code to perform analysis in paper: The Effect of Mutation Subtypes on the Allele Frequency Spectrum and Population Genetics Inference. Note AFS for each of the 96 mutation subtypes can be downloaded from: https://zenodo.org/record/6643451#.YtWeOuzMI-Q

Files include:
1) Main_Analysis.R: Main script. Imports AFS for each subtype and computes various summary statistics (proportion of singletons to doubletons, Tajima's D, etc).
2) selection_MST_abundant.R: Computes Tajima's D in 100kb windows 
3) Dadi_bottleneck.py and Dadi_expGrowth.py: Code to run DaDi 
4) D2_estimator.R and run_fastsimcoal.R : Computes new D-2 estimator and simulates neutral SFS for each subtype 
