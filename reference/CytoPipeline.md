# CytoPipeline package

`CytoPipeline` is a package that provides support for automation and
visualization of flow cytometry data analysis pipelines. In the current
state, the package focuses on the preprocessing and quality control
part.

The framework is based on two main S4 classes, i.e. `CytoPipeline` and
`CytoProcessingStep`. The `CytoProcessingStep` defines the link between
pipeline step names and corresponding R functions that are either
provided in the CytoPipeline package itself, or exported from a third
party package, or coded by the user her/himself. The processing steps
need to be specified centrally and explicitly using either a json input
file or through step by step creation of a `CytoPipeline` object with
dedicated methods.

After having run the pipeline, obtained results at all steps can be
retrieved and visualized thanks to file caching (the running facility
uses a `BiocFileCache` implementation). The package provides also
specific visualization tools like pipeline workflow summary display, and
1D/2D comparison plots of obtained flowFrames at various steps of the
pipeline.

For a step by step example using `CytoPipeline`, please have a look at
the vignette!

## See also

[CytoPipelineClass](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md),
[CytoProcessingStep](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)

## Author

**Maintainer**: Philippe Hauchamps <philippe.hauchamps@uclouvain.be>
([ORCID](https://orcid.org/0000-0003-2865-1852))

Authors:

- Laurent Gatto <laurent.gatto@uclouvain.be>
  ([ORCID](https://orcid.org/0000-0002-1520-2268))

Other contributors:

- Dan Lin <dan.8.lin@gsk.com> \[contributor\]
