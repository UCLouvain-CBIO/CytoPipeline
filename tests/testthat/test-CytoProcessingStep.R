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


# obtain OMIP021UTSamples, light-weight version used specifically for these 
# unit tests
path <- system.file("scripts",
                    package = "CytoPipeline"
)

source(file.path(path,"MakeOMIP021UTSamples.R"))

test_that("CytoProcessingStep basics works", {
    ps <- CytoProcessingStep("summing step", sum)
    
    psName <- getCPSName(ps)
    expect_equal(psName, "summing step")
    
    psFUN <- getCPSFUN(ps)
    expect_true(is.primitive(psFUN))
    
    psARGS <- getCPSARGS(ps)
    expect_identical(psARGS, list())
    
    expect_error(show(ps), NA)
    
    res <- executeProcessingStep(ps, 1:10)
    expect_equal(res, 55)
})

test_that("CytoProcessingStep wrong function works", {
    ps <- CytoProcessingStep("dummy step", "mistake_fun")
    
    psName <- getCPSName(ps)
    expect_equal(psName, "dummy step")
    
    expect_error(executeProcessingStep(ps, 1:10),
                 regexp = "not found")
})

test_that("CytoProcessingStep exports and imports work", {
    # case of a primitive
    ps <- CytoProcessingStep("summing step", sum)
    js_str <- as.json.CytoProcessingStep(ps)
    
    ps2 <- from.json.CytoProcessingStep(js_str)
    
    res <- executeProcessingStep(ps2, 1:10)
    expect_equal(res, 55)
    
    # case of a generic function
    ps <- CytoProcessingStep("median step", stats::median)
    
    js_str <- as.json.CytoProcessingStep(ps)

    ps2 <- from.json.CytoProcessingStep(js_str)
    
    res <- executeProcessingStep(ps2, 1:10)
    expect_equal(res, 5.5)
    
    # other case
    ps <- CytoProcessingStep("compensate step", "compensateFromMatrix")
    
    ff <- executeProcessingStep(ps, OMIP021UTSamples[[1]])
    res <- sum(flowCore::exprs(ff)[,"FSC-A"])
    expect_equal(res, 12553542.8)
    
    js_str <- as.json.CytoProcessingStep(ps)
    ps2 <- from.json.CytoProcessingStep(js_str)
     
    ff <- executeProcessingStep(ps2, OMIP021UTSamples[[1]])
    res <- sum(flowCore::exprs(ff)[,"FSC-A"])
    expect_equal(res, 12553542.8)
    
    # not yet implemented case (non generic, non primitive function as object)
    ps <- CytoProcessingStep("compensate step", compensateFromMatrix)
    expect_error(as.json.CytoProcessingStep(ps), 
                 regexp = "does not work")
                                        
})
