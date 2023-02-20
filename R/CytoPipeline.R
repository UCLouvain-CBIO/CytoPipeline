# CytoPipeline - Copyright (C) <2022>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoPipeline) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).


#' @title CytoPipeline package
#' 
#' @name CytoPipeline
#'
#' @rdname CytoPipeline
#' 
#' @seealso [CytoPipeline::CytoPipelineClass], [CytoPipeline::CytoProcessingStep]
#'
#' @description
#'
#' `CytoPipeline` is a package that provides support for automation and 
#' visualization of flow cytometry data analysis pipelines. In the current 
#' state, the package focuses on the preprocessing and quality control part.   
#' 
#' The framework is based on two main S4 classes, i.e. `CytoPipeline` and 
#' `CytoProcessingStep`. The `CytoProcessingStep` defines the link between
#' pipeline step names and corresponding R functions that are either provided 
#' in the CytoPipeline package itself, or exported from a third party package,
#' or coded by the user her/himself. The processing steps need to be specified 
#' centrally and explicitly using either a json input file or through 
#' step by step creation of a `CytoPipeline` object with dedicated methods.  
#' 
#' After having run the pipeline, obtained results at all steps 
#' can be retrieved and visualized thanks to file caching 
#' (the running facility uses a `BiocFileCache` implementation). 
#' The package provides also specific visualization tools like 
#' pipeline workflow summary display, and 1D/2D comparison plots 
#' of obtained flowFrames at various steps of the pipeline.
#' 
#' For a step by step example using `CytoPipeline`, please have a look 
#' at the vignette!
NULL