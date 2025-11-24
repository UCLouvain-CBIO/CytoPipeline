# Package index

## All functions

- [`show(`*`<CytoPipeline>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`CytoPipeline(`*`<missing>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`CytoPipeline(`*`<list>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`CytoPipeline(`*`<character>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`as.list(`*`<CytoPipeline>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`experimentName()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`` `experimentName<-`() ``](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`sampleFiles()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`` `sampleFiles<-`() ``](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`pData()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`` `pData<-`() ``](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`sampleDisplayNames()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  [`sampleNameFromDisplayName()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoPipelineClass.md)
  : CytoPipeline class
- [`CytoProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`show(`*`<CytoProcessingStep>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`executeProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`getCPSName()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`getCPSFUN()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`getCPSARGS()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`as.list(`*`<CytoProcessingStep>`*`)`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`as.json.CytoProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  [`from.json.CytoProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/CytoProcessingStep.md)
  : Cyto Processing step
- [`OMIP021Samples`](https://uclouvain-cbio.github.io/CytoPipeline/reference/OMIP021Samples.md)
  : OMIP021Samples dataset
- [`aggregateAndSample()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/aggregateAndSample.md)
  : Aggregate and sample multiple flow frames of a flow set together
- [`appendCellID()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/appendCellID.md)
  : append 'Original_ID' column to a flowframe
- [`applyScaleTransforms()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/applyScaleTransforms.md)
  : apply scale transforms
- [`areFluoCols()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/areFluoCols.md)
  : find flow frame columns that represent fluorochrome channel
- [`areSignalCols()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/areSignalCols.md)
  : find flow frame columns that represent true signal
- [`compensateFromMatrix()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/compensateFromMatrix.md)
  : compensation of fcs file(s) from matrix
- [`computeScatterChannelsLinearScale()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/computeScatterChannelsLinearScale.md)
  : compute linear transformation of scatter channels found in ff, based
  on 5% and 95% of referenceChannel, set as target. If there is a
  transformation defined in transList for referenceChannel, it is
  applied first, before computing quantiles. Then the computed linear
  transformations (or each scatter channel) are added into the
  transfo_list. -A channels are computed, and same linear transformation
  is then applied to corresponding -W and -H channels (if they exist in
  ff).
- [`estimateScaleTransforms()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/estimateScaleTransforms.md)
  : estimates scale tranformations
- [`execute()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/execute.md)
  : executing CytoPipeline object
- [`export2JSONFile()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/exportCytoPipeline.md)
  : exporting CytoPipeline objects
- [`findTimeChannel()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/findTimeChannel.md)
  : find time channel in flowSet/flowFrame
- [`getAcquiredCompensationMatrix()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/getAcquiredCompensationMatrix.md)
  : extract compensation matrix from a flowCore::flowFrame
- [`getChannelNamesFromMarkers()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/getChannelNamesFromMarkers.md)
  : get channel names from markers
- [`getFCSFileName()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/getFCSFileName.md)
  : get fcs file name
- [`getTransfoParams()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/getTransfoParams.md)
  : get tranformation parameters for a specific channel
- [`ggplotEvents()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/ggplotEvents.md)
  : plot events in 1D or 2D, using ggplot2
- [`ggplotFilterEvents()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/ggplotFilterEvents.md)
  : plot filtered events in 2D, using ggplot
- [`ggplotFlowRate()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/ggplotFlowRate.md)
  : plot flow rate as a function of time, using ggplot2
- [`addProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`removeProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`getNbProcessingSteps()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`getProcessingStep()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`getProcessingStepNames()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`cleanProcessingSteps()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  [`showProcessingSteps()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/handlingProcessingSteps.md)
  : handling processing steps in CytoPipeline objects
- [`getCytoPipelineExperimentNames()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`getCytoPipelineObjectFromCache()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`getCytoPipelineObjectInfos()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`getCytoPipelineFlowFrame()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`getCytoPipelineScaleTransform()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`plotCytoPipelineProcessingQueue()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  [`collectNbOfRetainedEvents()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/inspectCytoPipelineObjects.md)
  : inspect CytoPipeline results objects
- [`deleteCytoPipelineCache()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/interactingWithCytoPipelineCache.md)
  [`buildCytoPipelineFromCache()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/interactingWithCytoPipelineCache.md)
  [`checkCytoPipelineConsistencyWithCache()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/interactingWithCytoPipelineCache.md)
  : interaction between CytoPipeline object and disk cache
- [`qualityControlFlowAI()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/qualityControlFlowAI.md)
  : perform QC with flowAI
- [`qualityControlPeacoQC()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/qualityControlPeacoQC.md)
  : perform QC with PeacoQC
- [`readRDSObject()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/readRDSObject.md)
  : read RDS object
- [`readSampleFiles()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/readSampleFiles.md)
  : Read fcs sample files
- [`removeChannels()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/removeChannels.md)
  : remove channels from a flowFrame
- [`removeDeadCellsManualGate()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/removeDeadCellsManualGate.md)
  : remove dead cells from a flowFrame using manual gating
- [`removeDebrisManualGate()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/removeDebrisManualGate.md)
  : remove debris from a flowFrame using manual gating
- [`removeDoubletsCytoPipeline()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/removeDoubletsCytoPipeline.md)
  : remove doublets from a flowFrame, using CytoPipeline custom
  algorithm
- [`removeMarginsPeacoQC()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/removeMarginsPeacoQC.md)
  : remove margin events using PeacoQC
- [`resetCellIDs()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/resetCellIDs.md)
  : reset 'Original_ID' column in a flowframe
- [`runCompensation()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/runCompensation.md)
  : compensate with additional options
- [`singletsGate()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/singletsGate.md)
  : Clean doublet events from flow cytometry data
- [`subsample()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/subsample.md)
  : sub-sampling of a flowFrame
- [`updateMarkerName()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/updateMarkerName.md)
  : update marker name of a given flowFrame channel
- [`writeFlowFrame()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/writeFlowFrame.md)
  : write flowFrame to disk
