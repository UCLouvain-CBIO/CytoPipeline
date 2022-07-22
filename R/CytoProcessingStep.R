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
#' executeProcessingStep(ps, 1:10)
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
               if (is(res, "try-error"))
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
  new("CytoProcessingStep", name = name, FUN = FUN, ARGS = ARGS)
}

.cat_fun <- function(x) {
  if (is.function(x))
    "user-provided function"
  else paste(x, collapse = ", ")
}

#' @rdname CytoProcessingStep
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

#' @param object `CytoProcessingStep` object.
#'
#' @param ... optional additional arguments to be passed along.
#'
#' @rdname CytoProcessingStep
#'
#' @md
#'
#' @export
executeProcessingStep <- function(object, ...) {
  if (!is(object, "CytoProcessingStep"))
    stop("'object' should be a 'CytoProcessingStep' object!")
  do.call(object@FUN, args = c(list(...), object@ARGS))
}

#' @param object `CytoProcessingStep` object.
#'
#' @rdname CytoProcessingStep
#'
#' @export
getName <- function(object) {
  if (!is(object, "CytoProcessingStep"))
    stop("'object' should be a 'CytoProcessingStep' object!")
  return(object@name)
}

#' @param object `CytoProcessingStep` object.
#'
#' @rdname CytoProcessingStep
#'
#' @export
as.list.CytoProcessingStep <- function(x, ...) {
  if (!is(x, "CytoProcessingStep"))
    stop("'object' should be a 'CytoProcessingStep' object!")
  ll <- list()
  ll$name <- x@name
  ll$FUN <- x@FUN
  ll$ARGS <- x@ARGS
  return(ll)
}

#' @param object `CytoProcessingStep` object.
#'
#' @rdname CytoProcessingStep
#'
#' @export
as.json.CytoProcessingStep <- function(object, pretty = FALSE) {
  myList <- as.list.CytoProcessingStep(object)
  ret <- jsonlite::toJSON(myList, pretty = pretty, null = "null")
  return(as.character(ret))
}

#' @param object `CytoProcessingStep` object.
#'
#' @rdname CytoProcessingStep
#'
#' @export
from.json.CytoProcessingStep <- function(jsonString) {
  myList <- jsonlite::fromJSON(jsonString, simplifyDataFrame = FALSE)
  if(is.null(myList$name)) stop("name not found in Processing Step json string")
  if(is.null(myList$FUN)) stop("FUN not found in Processing Step json string")
  if(is.null(myList$ARGS)) stop("FUN not found in Processing Step json string")
  object <- CytoProcessingStep(myList$name, myList$FUN, myList$ARGS)
  return(object)
}

