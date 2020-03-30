# Src Hsp90 Deep Mutational Scan


Analyses done so far (VN)
- Replicate QC=0.89
- Comparing distribution of WT and synonymous variants between treatments shows that the distributions are similar in both treatments
- Plot of distributions of all point mutants and synonymous variants showed that the range of values is different in the two conditions
- Different radidicol sensitivity metrics
  - Ratio of the two scores (this doesn't work bc Enrich2 uses a log2 scoring so the distribution becomes very narrow at 0)
  - Difference (log2-log2, not sure if this is the best method)
  - Product (again, doesn't work bc of log2 scoring)
  - Residual
- Recoloring of PDB structures based on residuals and ratios (residual coloring seems more interpretable)
- Heatmaps of residuals, ratios (with and without trimming outliers)
- Dendrogram of heatmap of residuals
- Separation of residuals by SD (works bc residuals are normally distributed) to look at >2SD tails of residual distribution for highly sensitized mutants
- Correlagrams of SASA and by amino acid type did not reveal anything striking
