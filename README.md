<div align='justify'>

# TOM-complex-analysis
Data analysis and visualization for an article:

## A large-scale map of transient and constituive interactor of human TOM comple enabled by crosslinking

Mayra A. Borrero-Landazabal $^{1,6}$, Vanessa Linke $^{1,6}$, Tereza Kadavá $^{2}$, Zeshi Li $^{3}$, Piotr Draczkowski $^{4,5}$, 
Kacper Kaszuba&nbsp;$^{1}$, Remigiusz Serwa $^{1}$, Albert J. R. Heck $^{2*}$, Agnieszka Chacińska $^{1*}$

$^1$ IMol Polish Academy of Sciences, Warsaw, Poland <br>
$^2$ Biomolecular Mass Spectrometry and Proteomics, Bijvoet Center for Biomolecular 
Research and Utrecht Institute for Pharmaceutical Sciences, University of Utrecht, 
Padualaan 8, Utrecht 3584 CH, the Netherlands <br>
$^3$ Chemical Biology & Drug Discovery, Utrecht Institute for Pharmaceutical Sciences,
 University of Utrecht, Padualaan 8, Utrecht 3584 CH, the Netherlands <br>
$^4$ Department of Biochemistry and Biophysics, National Bioinformatics Infrastructure 
Sweden, Science for Life Laboratory, Stockholm University, Solna, Sweden <br>
$^5$ Department of Synthesis and Chemical Technology of Pharmaceutical Substances, 
Medical University of Lublin, Lublin, Poland <br>
$^6$ These authors contributed equally: Mayra A. Borrero-Landazabal, Vanessa Linke <br>
*Corresponding authors: Albert J. R. Heck (a.j.r.heck@uu.nl); Agnieszka Chacinska (a.chacinska@imol.institute)

## Abstract

The core architecture of the translocase of the outer mitochondrial membrane (TOM) that
serves as the central entry gate for nuclear-encoded proteins is evolutionarily conserved from
yeast to human. Here, we define the interactome of the human TOM complex that deals with a
vastly more complex environment. Encouraged by substantial overlap with known yeast
homologs, we expanded the interactome to include transient and labile interactors stabilized by 
a&nbsp;membrane-permeable crosslinker. Mapping of 24 unique inter-protein crosslinks to TOM core
subunits provided structural insight into interaction interfaces. Crosslinking further enhanced
recovery of peripheral components such as the receptor TOMM70 and several TOM-associated
quality control factors, including ATAD1, the human homolog of Msp1, recently found to prevent
complex clogging. We also mapped novel human TOM interactors, including FKBP8, an
immunophilin linked to mitophagy and apoptosis. Our findings provide the basis to uncover
human-specific regulatory mechanisms behind the complex human mitochondrial network.

# Content 📁

```
│   .gitattributes
│   README.md
│
└───code
        analysis.Rmd
        plotting.R
        prepare_the_data.R
        supplementary_file.R
```

# Description

This repository contains the code used to analyze the mass spectrometry data for TOM-PD. Each file is responsible for the following calculations:
1) The file, **`analysis.Rmd`**, contains the following elements of analysis:
    - The initial step involves the cleansing of missing data.
    - Subsequent to this step is the imputation of relative proteins.
    - The third step is the calculation of the fold change.
    - The final step is the conducting of a t-test.
2) The **`prepare_the_data.R`** is utilized for the preparation of the results from *1)* to generate the plots that are included in the article.
3) The **`plotting.R`** is utilized for the purpose of generating a plot.
4) The **`supplementary_file.R`** combines all of the utilized datasets into a single Excel file.

# Code authors 🖥️

1. Vanessa Linke - [vanilink](https://github.com/vanilink)
2. Kacper Kaszuba - [KacperKaszuba0608](https://github.com/KacperKaszuba0608)

</div>