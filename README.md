# Automated generation of gene set enrichment analyses 

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/augere.gsea.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/augere.gsea.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/augere.gsea/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/augere.gsea.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/augere.gsea.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/augere.gsea/)|

Implements pipeline functions to generate parametrized Rmarkdown reports for gene set enrichment analyses.
We support analyses based on pre-computed tables of differential expression results,
or by using **limma** to test for differential gene set activity from the original expression data.
Each report contains all of the R commands required to reproduce the analysis.
