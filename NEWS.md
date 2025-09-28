# CytoPipeline 1.9

## CytoPipeline 1.9.2
- sample files can now have duplicate base names 
(provided full paths are different)
- `pData<-` is now more liberal.   
1. It can accept new pData containing more rows than existing sample names 
(the corresponding subset of pData is taken).
2. It can accept pData with row names pointing to either sample file full paths 
or base file names
3. It can accept pData with no row names provided the number of rows correspond 
to the number of sample files. Row names are then set by default to sample 
file base names (if unique), or sample file full paths.

## CytoPipeline 1.9.1
- upgraded to GHA cache v4

# CytoPipeline 1.7
(no devel)

# CytoPipeline 1.5

## CytoPipeline 1.5.2
- updated unit test to account for flowAI version change
- now suggesting `CytoPipelineGUI` package

## CytoPipeline 1.5.1
- updated processing step argument matching using `phenoData`

# CytoPipeline 1.3

## CytoPipeline 1.3.6
- `execute()` now stores the nb of events retained at each pre-processing step,
to speed-up `collectNbOfRetainedEvents()`

## CytoPipeline 1.3.5
- added CITATION file

## CytoPipeline 1.3.4
- in execute(), when cache does already exist, make it clean before executing 
the pipeline steps (= preventing inconsistent cache upon crash/forced 
interruption)

## CytoPipeline 1.3.3
- added `collectNbOfRetainedEvents()` function

## CytoPipeline 1.3.2
- systematically override `pData` in cache upon execute() to allow running 
consistently running several times for increasing number of samples
- `sampleFiles<-` and `pData<-`: make sure that order of sample files 
follow the one of `pData` rownames if `pData` exists.
- added 'verbose' argument in `estimateScaleTransforms()`

## CytoPipeline 1.3.1
- refactored documentation files

# CytoPipeline 1.1

## CytoPipeline 1.1.5
- areSignalCols(), are FluoCols() can now accept a flowSet as input, 
on top of a flowFrame
- applyScaleTransforms() processing step has been improved (can take flowSet 
as input, checks channel concordance between transList and data object)

## CytoPipeline 1.1.4
- updated fcs files, that are at the source of OMIP021Samples dataset

## CytoPipeline 1.1.3
- in subSample(), renamed parameter 'nSamples' into 'nEvents', and added 
possibility for passing unused parameters, in order to support 
the use of the function as a processing step. Also amended the function
as to keep the original order of the events (keep chronology).
Finally, adds a 'keepOriginalCellIDs' parameter (default=TRUE).
- simplified the arguments of execute() related to the storage of the results
after last pre-processing step.

## CytoPipeline 1.1.2
- storage of phenoData into cache upon execution of CytoPipeline object
(and back into CytoPipeline object when re-built from cache)
- changed the default behaviour of estimateScaleTransforms() so that
the default method for scatter channels is now "none" instead of 
"linearQuantile"
- changed default behaviour of ggplotEvents() and ggplotFilterEvents(), when
logicle scale is used but no logicle parameters provided, these are now
estimated using flowCore::estimateLogicle(), instead of explicit default values

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
