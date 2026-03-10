# Simple single-cell analyses

Implements pipeline functions to generate parametrized Rmarkdown reports for simple single-cell (RNA sequencing) analyses.
This uses **scrapper** for routine steps such as quality control, normalization, feature selection, clustering and marker detection.
We also implement a pipeline for automatic cell type annotation against a labelled reference with **SingleR**.
Each report contains all of the R commands required to reproduce its analysis.
