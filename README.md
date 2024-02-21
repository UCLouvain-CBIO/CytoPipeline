## Automation and visualization of flow cytometry data analysis pipelines

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/UCLouvain-CBIO/CytoPipeline/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/UCLouvain-CBIO/CytoPipeline/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/CytoPipeline/coverage.svg?branch=main)](https://codecov.io/github/UCLouvain-CBIO/CytoPipeline?branch=main)
[![license](https://img.shields.io/badge/license-GPL3.0-blue)](https://opensource.org/licenses/GPL-3.0)

### What is CytoPipeline?

`CytoPipeline` is a package that provides support for automation and 
visualization of flow cytometry data analysis pipelines. In the current 
state, the package focuses on the preprocessing and quality control part.   

The framework is based on two main S4 classes, i.e. `CytoPipeline` and 
`CytoProcessingStep`. The `CytoProcessingStep` defines the link between
pipeline step names and corresponding R functions that are either provided in
the CytoPipeline package itself, or exported from a third party package,
or coded by the user her/himself. The processing steps need to be specified 
centrally and explicitly using either a json input file or through step by step 
creation of a `CytoPipeline` object with dedicated methods.  

After having run the pipeline, obtained results at all steps can be retrieved 
and visualized thanks to file caching (the running facility uses a 
`BiocFileCache` implementation). The package provides also specific 
visualization tools like pipeline workflow summary display, and 1D/2D comparison 
plots of obtained flowFrames at various steps of the pipeline.

### License

The `CytoPipeline` code is provided under [GPL license version 3.0 or 
higher](https://opensource.org/licenses/GPL-3.0). The documentation, 
including the manual pages and the vignettes, are distributed under a [CC BY-SA 
4.0 license](https://creativecommons.org/licenses/by-sa/4.0/).

### Citation

If you use `Cytopipeline` in your research, please use the following citation:

>Hauchamps P, Bayat B, Delandre S, Hamrouni M, Toussaint M, Temmerman S, 
> Lin D, Gatto L (2024).
> “CytoPipeline and CytoPipelineGUI: a Bioconductor R package suite 
> for building and visualizing automated pre-processing pipelines 
> for flow cytometry data.” _BMC Bioinformatics_, *25*(1), 80. 
> doi:10.1186/s12859-024-05691-z <https://doi.org/10.1186/s12859-024-05691-z>.

or run `citation("CytoPipeline")` to get the bibtex entry.
