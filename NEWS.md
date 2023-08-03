# CytoPipeline 1.1

## CytoPipeline 1.1.3
- in subSample(), renamed parameter 'nSamples' into 'nEvents', and added 
possibility for passing unused parameters, in order to support 
the use of the function as a processing step.

## CytoPipeline 1.1.2
- storage of phenoData into cache upon execution of CytoPipeline object
(and back into CytoPipeline object when re-built from cache)
- changed the default behaviour of estimateScaleTransforms() so that
the default method for scatter channels is now "none" instead of 
"linearQuantile"
- changed default behaviour of ggplotEvents() and ggplotFilterEvents(), when
logicle scale is used but no logicle parameters provided, these are now
estimated using flowCore::estimatLogicle(), instead of explicit default values

## CytoPipeline 1.1.1
- tiny modifications to support upgrade to Bioc 3.18

## CytoPipeline 0.99

## CytoPipeline 0.99.6
- corrected the OMIP021Samples fcs data in order to keep the original file
name
- bug correction: error message on execution with no sample file
- added `phenoData` slot in CytoPipeline object
- updated `readSamples()` to allow passing a `pData` parameters
- updated `compensateFromMatrix()` to allow passing a mapping 
based on a `pData` variable
- updated `readSamples()` to allow selecting a random number of samples and
removed `selectSamples()`
- vignette with demo and links to videos

## CytoPipeline 0.99.5
- reactivated unit tests for ggplot2 objects
- added man page for CytoPipeline package
- a few modifs in the vignette related to Bioc review process
- replaced withr::local_tempdir() by base::tempdir()
- removed extraneous whitespaces in CytoPipeline show() method
- removed `LazyData: true` in DESCRIPTION file 
- replaced `paste0(path, "/", filename)` by `file.path(path, filename)`
- updated License field in DESCRIPTION file

## CytoPipeline 0.99.4

- improved CytoPipeline constructors (`experimentName` and `sampleFiles` are
now parameters of all constructor version)
- centralized the production of standard outputs during pipeline execution, 
set all tuning parameters in execute() instead of slots 
in CytoPipeline object.  

## CytoPipeline 0.99.3

- some minor changes for BiocCheck()

## CytoPipeline 0.99.2

- removed dependencies to a number of packages, moved corresponding 
implementations of CytoProcessingSteps (wrappers) into `CytoPipelineUtils` 
package

## CytoPipeline 0.99.1

- Maintenance due to Bioc version change (3.17)
- removed use of openCyto::gate_tail() (disappeared w/o deprecation), replaced
by flowDensity::deGate()
- implemented export of pre-processed file (writeFlowFrame as a 
CytoProcessingStep implementation)
- extended readSampleFiles : mapping between channels and markers
- selectRandomSamples (new CytoProcessing step implementation)


## CytoPipeline 0.99.0

- Prior to Bioconductor submission
