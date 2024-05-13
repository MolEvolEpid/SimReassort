##' internals
##'
##'
##' @name internals
##'
##' @param x \code{gpsim} object
##' @param prune logical
##' @param init.include logical
##' @param tr1 character, newick
##' @param tr2 character, newick
##' @param phy \code{phylo} object
##' @param p node label
##' @param df data frame
##'
##'
##'
##' @importFrom stringr str_extract_all
##' @importFrom ggtree fortify
##' @importFrom dplyr arrange filter
##'
##'
##' @keywords internals


##' @rdname internals
current_popsize <- function (x) {
  x |>
    getInfo(prune=FALSE, lineages=TRUE) |>
    getElement("lineages") |>
    getElement(1) -> lineages

  return (lineages$lineages[nrow(lineages)-1])
}

##' @rdname internals
get_roottime <- function (x, prune=FALSE, init.include=FALSE) {
  if (!prune && current_popsize(x) < 2) {
    root.time <- 0
  } else {
    trs <- getInfo(x,prune=prune, tree=TRUE, hide=TRUE)$tree
    if (!init.include)  trs <- trs[-3]
    lapply(trs, function (tr) {
      read.tree(text=tr) |>
        fortify(ladderize=TRUE) |>
        arrange(x) |>
        filter(grepl("g_",label)) |>
        getElement("x") |> getElement(1)
    }) |> unlist() |> min() -> root.time
  }
  return (root.time)
}

##' @rdname internals
get_tiptimes <- function (x, prune=FALSE) {
  if (!prune && current_popsize(x) < 2) {
    tip.times <- NA
  } else {
    tr <- getInfo(x, prune=prune, tree=TRUE, hide=TRUE)$tree[3]
    read.tree(text=tr) |>
      fortify(ladderize=TRUE) |>
      arrange(x) |>
      filter(grepl("r",label)) |>
      getElement("x") -> tip.times
  }
  return (tip.times)
}

##' @rdname internals
get_discrepancies <- function (tr1, tr2) {
  lapply(c(tr1,tr2), function (tr) {
    str_extract_all(tr, "g_\\d_\\d+") |> unlist()
  }) -> coalset

  coal.common <- do.call(intersect, coalset)
  coal.all <- do.call(union, coalset)

  setdiff(coal.all, coal.common)
}

##' @rdname internals
drop_node <- function (phy, node) {
  lapply(node, function(nd) {
    sum(grepl(paste0(nd,"$"),phy$node.label))
  }) |> unlist() -> node.ind
  if (any(node.ind<1)) {
    warning("some nodes were not in the tree: skip them.")
    node <- node[node.ind>0]
  }
  if (length(node)==phy$Nnode) {
    warning("All nodes are dropped. Return NULL.")
    return (NULL)
  }
  df <- phy |> fortify(ladderize=TRUE)
  for (nnode in node) {
    node.num <- df |> filter(grepl(paste0(nnode,"$"),label)) |> getElement("node")
    child.num <- df |> filter(parent==node.num) |> getElement("node")
    if (length(child.num) > 1) {
      warning(paste0("the specified node ", nnode, " has more than 1 children: it cannot be removed, skip."))
      next
    } else {
      parent.num <- df |> filter(node==node.num) |> getElement("parent")
      df[df$node==child.num,"parent"] <- parent.num
      df[df$node==child.num,"branch.length"] <- sum(df |> filter(node%in%c(node.num,child.num))|>getElement("branch.length"))
      df <- df |> filter(node!=node.num)
    }
  }
  return (read.tree(text=df |> filter(!label %in% c("","i_NA_NA")) |> write_nwk()))
}

##' @rdname internals
nwk_recur <- function(p, df) {
  node <- df[df$node == p,] |> getElement("label")
  childs <- df[df$parent == p,] |> getElement("node")
  len <- df[df$node == p,] |> getElement("branch.length")
  paste0("(",
         ifelse(length(childs) > 0, paste0(lapply(childs,function(c) {nwk_recur(c,df)}),collapse=","),""),
         ")",
         node,":",len)
}

##' @rdname internals
write_nwk <- function(df) {
  p <- df |> arrange(x) |> getElement("parent") |> getElement(1)
  paste0(gsub("\\(\\)","",nwk_recur(p,df)),";") -> out
  gsub(":;",";",out)
}
