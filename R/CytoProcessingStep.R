# CytoPipeline - Copyright (C) <2022> <Université catholique de Louvain (UCLouvain), Belgique>
#   
#   Description and complete License: see LICENSE file.
# 
# This program (CytoPipeline) is free software: 
#   you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).

setClassUnion("characterOrFunction", c("character", "function"))

#' @title Cyto Processing step
#'
#' @aliases CytoProcessingStep-class characterOrFunction-class
#'
#' @description
#'
#' Class containing the function and arguments to be applied in a lazy-execution
#' framework.
#'
#' Objects of this class are created using the `CytoProcessingStep()` function.
#' The processing step is executed with the `executeProcessingStep()` function.
#'
#' @details
#'
#' This object contains all relevant information of a data analysis processing
#' step, i.e. the function and all of its arguments to be applied to the data.
#'
#' @return The `CytoProcessingStep` function returns and object of type
#'     `CytoProcessingStep`.
#'
#' @exportClass CytoProcessingStep
#'
#' @examples
#'
#' ## Create a simple processing step object
#' ps <- CytoProcessingStep("summing step", sum)
#'
#' getName(ps)
#'
#' executeProcessingStep(ps, seq_len(10))
#'
#'
#' @name CytoProcessingStep
NULL

setClass("CytoProcessingStep",
         representation = representation(
           name = "character",
           FUN = "characterOrFunction",
           ARGS = "list"
         ),
         prototype = prototype(
           name = character(),
           ARGS = list(),
           FUN = character()
         ),
         validity = function(object) {
           msg <- character()
           ## Fails with un-exported functions.
           if (length(object@FUN)) {
             if (!is.function(object@FUN)) {
               res <- try(match.fun(object@FUN), silent = TRUE)
               if (methods::is(res, "try-error"))
                 msg <- c(msg, paste0("Function '", object@FUN,
                                      "' not found."))
             }
           }
           if (length(msg))
             msg
           else TRUE
         })

#' @param name `character` denoting a name to the step, which can be different
#' from the function name
#'
#' @param FUN `function` or `character` representing a function name.
#'
#' @param ARGS `list` of arguments to be passed along to `FUN`.
#'
#' @rdname CytoProcessingStep
#'
#' @export
CytoProcessingStep <- function(name = character(), FUN = character(),
                               ARGS = list())  {
  if (missing(FUN))
    FUN <- character()
  methods::new("CytoProcessingStep", name = name, FUN = FUN, ARGS = ARGS)
}

.cat_fun <- function(x) {
  if (is.function(x))
    "user-provided function"
  else paste(x, collapse = ", ")
}

#' @rdname CytoProcessingStep
#' 
#' @param object a `CytoProcessingStep` object.
#'
#' @exportMethod show
setMethod("show", "CytoProcessingStep", function(object) {
  cat("Object of class \"", class(object), "\"\n", sep = "")
  cat(" Name: ", object@name, "\n", sep = "")
  cat(" Function: ", .cat_fun(object@FUN), "\n", sep = "")
  args <- object@ARGS
  if (length(args) > 0) {
    cat(" Arguments:\n")
    for (i in seq_along(args)) {
      cat("  o ", names(args)[i], " = ", .cat_fun(args[[i]]),
          "\n", sep = "")
    }
  }
})

#' @param x a `CytoProcessingStep` object.
#'
#' @param ... optional additional arguments to be passed along.
#'
#' @rdname CytoProcessingStep
#'
#' @md
#'
#' @export
executeProcessingStep <- function(x, ...) {
  if (!methods::is(x, "CytoProcessingStep"))
    stop("'x' should be a 'CytoProcessingStep' object!")
  do.call(x@FUN, args = c(list(...), x@ARGS))
}

#' @param x a `CytoProcessingStep` object.
#'
#' @rdname CytoProcessingStep
#'
#' @export
getName <- function(x) {
  if (!methods::is(x, "CytoProcessingStep"))
    stop("'x' should be a 'CytoProcessingStep' object!")
  return(x@name)
}

#' @param x a `CytoProcessingStep` object.
#' @param ... other arguments (not used)
#'
#' @rdname CytoProcessingStep
#'
#' @export
as.list.CytoProcessingStep <- function(x, ...) {
  if (!methods::is(x, "CytoProcessingStep"))
    stop("'x' should be a 'CytoProcessingStep' object!")
  ll <- list()
  ll$name <- x@name
  ll$FUN <- x@FUN
  ll$ARGS <- x@ARGS
  return(ll)
}

#' @param x a `CytoProcessingStep` object.
#' @param pretty formatting set-up (see jsonlite::toJSON doc)
#'
#' @rdname CytoProcessingStep
#'
#' @export
as.json.CytoProcessingStep <- function(x, pretty = FALSE) {
  myList <- as.list.CytoProcessingStep(x)
  ret <- jsonlite::toJSON(myList, pretty = pretty, null = "null")
  return(as.character(ret))
}

#' @param jsonString a `character()` containing a JSON string.
#'
#' @rdname CytoProcessingStep
#'
#' @export
from.json.CytoProcessingStep <- function(jsonString) {
  myList <- jsonlite::fromJSON(jsonString, simplifyDataFrame = FALSE)
  if (is.null(myList$name)) stop("name not found in Processing Step json string")
  if (is.null(myList$FUN)) stop("FUN not found in Processing Step json string")
  if (is.null(myList$ARGS)) stop("FUN not found in Processing Step json string")
  object <- CytoProcessingStep(myList$name, myList$FUN, myList$ARGS)
  return(object)
}

