#### The Goal:
The following example is a comparison between 7 of the most used single-sample Gene Set Analysis methods, in order to find the methods with the best sensitivity, specificity, and precision.
* [GSVA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/) - "Gene set variation analysis". 	 
* [SSGSEA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2783335/) - "Single sample gene set enrichment analysis".
* [PLAGE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1261155/) - "Pathway level analysis of gene expression".
* [ZSCORE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2563693/) - a classification method based on pathway activities inferred for each patient. 
* [individPath](https://academic.oup.com/bib/article/17/1/78/1742633) - "Individualized Pathway Coordination Analysis". 
* [GRAPE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5485588/) - "Gene ranking analysis of pathway expression".
* [Pathifier](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3631698/) - an algorithm designed to model the heterogeneity in samples based on their genomic data.


#### The Procedure:
The following is our benchmarking procedure: 
We chose a group of transcriptomic datasets related to a given respiratory disease (see Tab 1 and 2) and a group of "target pathways", which are the specific pathways corresponding to 
that disease (see Tab 2). We also selected a pathway database which contains the target pathways and many others (see Tab 2). Then we apply all seven methods to our datasets. The results 
are compared to the group of target pathways, which is considered a gold standard. After that comparison, we count true and false positives, as well as true and false negatives; such values 
allow us to measure sensitivity, specificity, and precision per dataset and per method. Out last step is to plot the distribution of all datasets sensitivities per method, and compare between methods.
The procedure (and R code) to reproduce this example can be found here: 
* [GSVA Tutorial](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/GSVA_GSE10245.ipynb) 
* [Pathifier Tutorial](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/pathifier_GSE10245.ipynb)
* [GRAPE Tutorial](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/GRAPE_GSE10245.ipynb)
