# Simple single-cell analyses

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/augere.solo.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/augere.solo.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/augere.solo/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/augere.solo.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/augere.solo.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/augere.solo/)|

Implements pipeline functions to generate parametrized Rmarkdown reports for simple single-cell (RNA sequencing) analyses.
This uses **scrapper** for routine steps such as quality control, normalization, feature selection, clustering and marker detection.
We also implement a pipeline for automatic cell type annotation against a labelled reference with **SingleR**.
Each report contains all of the R commands required to reproduce its analysis.
