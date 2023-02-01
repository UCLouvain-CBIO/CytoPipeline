# CytoPipeline 0.99

## CytoPipeline 0.99.5

(no changes yet)

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
