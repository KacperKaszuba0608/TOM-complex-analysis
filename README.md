<div align='justify'>

# TOM-complex-analysis
Data analysis and visualization for an article:

* ## A dynamic gateway: uncovering the expanded human TOM complex interactome and its regulatory complexity 
* ## The Dynamic Network of the Human TOM Complex.

**Authors:**
Mayra A. Borrero-Landazabal $^{1}$†, Vanessa Linke $^{1}$†‡, Tereza Kadavá $^{2,3}$, Zeshi Li $^{4}$, Piotr Draczkowski $^{5,6}$, <br>
Kacper Kaszuba $^{1}$, Lea Bertgen $^{1}$, Remigiusz A. Serwa $^{1}$, Albert J. R. Heck $^{2,3}$\*, Agnieszka Chacinska $^{1}$\*


**Affiliations:** 
$^1$ IMol Polish Academy of Sciences, Warsaw, Poland. <br>
$^2$ Biomolecular Mass Spectrometry and Proteomics, Bijvoet Center for Biomolecular Research and Utrecht Institute for Pharmaceutical Sciences, University of Utrecht, Padualaan 8, Utrecht 3584 CH, the Netherlands. <br>
$^3$ Netherlands Proteomics Center, Padualaan 8, Utrecht 3584 CH, the Netherlands. <br>
$^4$ Chemical Biology & Drug Discovery, Utrecht Institute for Pharmaceutical Sciences, Utrecht University, Universiteitsweg 99, Utrecht 3584 CG, the Netherlands. <br>
$^5$ Department of Biochemistry and Biophysics, National Bioinformatics Infrastructure Sweden, Science for Life Laboratory, Stockholm University, Solna, Sweden. <br>
$^6$ Department of Synthesis and Chemical Technology of Pharmaceutical Substances, Medical University of Lublin, Lublin, Poland. <br>
\* Corresponding authors: a.j.r.heck@uu.nl (Albert J. R. Heck ); a.chacinska@imol.institute (Agnieszka Chacinska ) <br>
† These authors contributed equally. <br>
‡ Present address: Mass Spectrometry Facility, IN-MOL-CELL Infrastructure, International Institute of Molecular and Cell Biology in Warsaw, 02‑109, Warsaw, Poland

## Abstract

The translocase of the outer mitochondrial membrane (TOM) is the conserved entry 
gate for nuclear-encoded proteins. While structurally similar from yeast to humans, 
the human TOM complex operates in a cellular environment of vastly greater complexity. 
Here, we present a high-confidence map of the human TOM interactome using a membrane-permeable
crosslinker to capture both stable and transient interactors. Excitingly, alongside 
extensive overlap with known yeast partners, we uncover a set of human-specific 
interactors including regulatory factors and newly identified TOM-associated proteins. 
Mapping unique inter-protein crosslinks reveals conformational flexibility of 
the receptor TOM20 and enhanced recovery of peripheral components such as TOM70 
and several associated quality control factors. Notably, we identify FKBP8 as a 
novel human-specific interactor that binds multiple TOM subunits and promotes 
organization of the complex. Our work redefines the human TOM complex as a dynamic, 
multifaceted hub coordinating biogenesis, quality control, and signaling. This 
expanded TOM landscape offers a rich resource for exploring mitochondrial regulation 
in health and disease.

# Content 📁

```
│   .gitattributes
│   README.md
│
└───code
        analysis.Rmd
        fkbp8_analysis.R
        fxn.R
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
2) The **`fkbp8_analysis.R`** is utilized for the analysis of the FKBP8 protein, which includes the following steps:
    - The initial step involves the cleansing of missing data.
    - Subsequent to this step is the imputation of relative proteins.
    - The third step is the calculation of the fold change.
    - The final step is the conducting of a t-test.
3) The **`fxn.R`** contains the functions that are utilized in the **`fkbp8_analysis.R`** file.
4) The **`plotting.R`** is utilized for the purpose of generating a plot.
5) The **`prepare_the_data.R`** is utilized for the preparation of the results from *1)* to generate the plots that are included in the article.
6) The **`supplementary_file.R`** combines all of the utilized datasets into a single Excel file.

# Code authors 🖥️

1. Vanessa Linke - [vanilink](https://github.com/vanilink)
2. Kacper Kaszuba - [KacperKaszuba0608](https://github.com/KacperKaszuba0608)

</div>