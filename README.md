# Smith_SpC_2018

A large-scale spectral counting data analysis using PAW pipeline, edgeR, R, and Jupyter notebooks.

### Phil Wilmarth, OHSU
#### February 2019

The data is from a [recent study](https://www.sciencedirect.com/science/article/pii/S0002939418301193) where human retinal and choroidal endothelial cells were compared. The study was 5 donor eyes where retinal and choroidal cells were collected and cultured in a paired design. The cell lysates from each of the 10 cell cultures were profiled using large-scale separations with a fast-scanning linear ion trap. There were about half a million MS2 scans per sample for a dataset size of a little over 5 million. The data are available at the PRIDE archive ([PXD005972](https://www.ebi.ac.uk/pride/archive/projects/PXD005972)).

There is a direct link to the rendered [notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/Smith_2018_edgeR.html).

---

A few relevant files from the archive are present in the repository:

- analysis_overview.pptx - summary of the data analysis steps
- quant_protein_summary_8.txt - a grouped protein summary file
- edgeR_input.txt - data extracted from results file for edgeR analysis
- edgeR_results.txt - statistical testing results
- HCEC_HREC_quant_protein_summary_sprot.xlsx - final summary sheet
  - proteomics data from quant_protein_summary_8.txt
  - statistical results from edgeR
  - extra protein annotations
