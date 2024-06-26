##' getInfo
##'
##' Retrieve information from genealogy process simulation
##'
##' @name getInfo
##'
##' @param object \code{gpsim} object.
##' @param prune logical; prune the genealogy?
##' @param obscure logical; obscure the demes?
##' @param hide logical; hide reassortment events?
##' @param time logical; return the current time?
##' @param t0 logical; return the zero-time?
##' @param tree logical; return the tree?
##' @param compact logical; return the tree in compact representation?
##' @param description logical; return the description?
##' @param yaml logical; return the structure in YAML format?
##' @param structure logical; return the structure in \R list format?
##' @param lineages logical; return the lineage-count function?
##'
##' @include package.R
##' @importFrom dplyr bind_cols
##' @importFrom tibble as_tibble
##'
##' @return
##' A list containing the requested elements, including any or all of:
##' \describe{
##'   \item{t0}{the initial time}
##'   \item{time}{the current time}
##'   \item{tree}{the genealogical tree, in Newick format}
##'   \item{description}{a human readable description of the state of the genealogy process}
##'   \item{yaml}{the state of the genealogy process in YAML format}
##'   \item{structure}{the state of the genealogy process in \R list format}
##'   \item{lineages}{a \code{\link[tibble]{tibble}} containing the lineage count function through time}
##' }
##'
##' @example examples/lbdpwr.R
##'
##' @rdname getinfo
##' @export
getInfo <- function (
    object, prune  = TRUE, obscure = TRUE, hide = FALSE,
    t0 = FALSE, time = FALSE,
    description = FALSE,
    structure = FALSE, yaml = FALSE,
    lineages = FALSE,
    tree = FALSE, compact = TRUE)
{
  x <- switch(
    paste0("model",as.character(attr(object,"model"))),
    modelLBDPwr = .Call(P_infoLBDPwr2,object,prune,obscure,hide,t0,time,
                        description,yaml,structure,lineages,tree,compact),
    model = stop("no model specified",call.=FALSE),
    stop("unrecognized model ",sQuote(attr(object,"model")),call.=FALSE)
  )
  if (!is.null(x$tree)) {
    x$tree <- gsub("nan","NA",x$tree)
  }

  if (!is.null(x$lineages)) {
    lapply(x$lineages, function(y) {
      n <- length(y$time)
      m <- length(y$count)/n
      if (m > 1L) {
        dig <- ceiling(log10(m))
        nm <- sprintf(paste0("deme%0",dig,"d"),seq_len(m))
      } else {
        nm <- "lineages"
      }
      bind_cols(
        time=y$time,
        y$count |>
          as.integer() |>
          array(dim=c(m,n),dimnames=list(nm,NULL)) |>
          t() |>
          as_tibble()
      )
    }) -> x$lineages
  }
  x
}
