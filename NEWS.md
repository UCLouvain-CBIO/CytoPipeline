# CytoPipeline 0.99

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
