##' Track (a set of) components from a tree
##'
##'
##' @name track
##'
##' @param parent.lab character; the branching node
##' @param child.lab character;
##' @param coal.lab character;
##' @param Re.lab character;
##' @param node.lab1 character;
##' @param node.lab2 character;
##' @param nwktree character; tree in newick formate
##' @param parent.numset numeric vector;
##' @param desc.numset numeric vector;
##' @param nodecol character;
##' @param incongr.nodes character vector;
##' @param reversed.nodes character vector;
##'
##' \describe{
##'   \item{parent.lab}{the specified parent branching node}
##'   \item{child.lab}{the specified child node}
##'   \item{coal.lab}{the node label of the coalescent node}
##'   \item{Re.lab}{the node label of the reassortment node}
##'   \item{node.lab1}{a node label}
##'   \item{node.lab2}{another node label}
##'   \item{nwktree}{the tree string in newick format}
##'   \item{parent.numset}{a set of numbers indicating the candicate ancestors}
##'   \item{desc.numset}{a set of numbers indicating the candicate descendants}
##'   \item{nodecol}{the specified type/color (i.e. 'r' for samples, 'g' for branching nodes, '#Ha1' and '#Ha2' for reassortment nodes) of the requested nodes}
##'   \item{incongr.nodes}{the discording coalescent nodes between trees}
##'   \item{reversed.nodes}{the reversed coalescent nodes between trees}
##' }
##'
##'
##' @importFrom ape read.tree
##' @importFrom ggtree fortify
##' @importFrom dplyr filter
##'
NULL


##' @rdname track
##' @return A string for the node/tip
first.descendant <- function (parent.lab, nwktree, desc.numset, nodecol) {
  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    filter(label!="") -> df
  df |> filter(grepl(paste0("^",parent.lab,"$"), label)) |> getElement("node") -> parent.node
  if (length(parent.node) < 1)  return (NULL)
  df |> filter(parent==parent.node) -> children
  if (grepl("(g|#)",parent.lab) & nrow(children) > 0) {
    if (any(grepl(nodecol, children |>
                  filter(gsub(paste0(nodecol,"(\\d)?_\\d_"),"", label) %in% desc.numset) |>
                  getElement("label")))
        ) {
      children |> filter(grepl(nodecol,label)) |>
        filter(gsub(paste0(nodecol,"(\\d)?_\\d_"),"",label) %in% desc.numset) |>
        filter(x==min(x)) |>
        getElement("label") -> desc.node
    } else if (any(grepl("(g|#)", children |> getElement("label")))) {
      lapply(children |> filter(grepl("g|#",label)) |> getElement("label"),
             function (par) {
               first.descendant(par, nwktree, desc.numset, nodecol)
             }) |>
        unlist() -> desc.nodes
      desc.node <- desc.nodes[which(as.numeric(gsub(paste0(nodecol,"(\\d)?_\\d_"),"",desc.nodes))==
                                min(as.numeric(gsub(paste0(nodecol,"(\\d)?_\\d_"),"",desc.nodes))))]
    } else {
      warning(paste0("no node/tip corresponds in the clade of node ", parent.lab,"."))
      desc.node <- NULL
    }
  } else {
    warning("The node is not a branching/reassortment node!")
    return (NULL)
  }
  return (desc.node)
}

##' @rdname track
##' @return A string for the node/tip
last.descendant <- function (parent.lab, nwktree, desc.numset, nodecol) {
  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    filter(label!="") -> df
  df |> filter(grepl(paste0("^",parent.lab,"$"), label)) |> getElement("node") -> parent.node
  if (length(parent.node) < 1)  return (NULL)
  df |> filter(parent==parent.node & !grepl("r", label)) |> getElement("label") -> children.labs

  if (length(children.labs) > 0) {
    lapply(children.labs, function (par) {
      last.descendant(par, nwktree, desc.numset, nodecol)
    }) |> unlist() -> desc.nodes

    if (length(desc.nodes) > 0) {
      desc.node <- desc.nodes[which(as.numeric(gsub(paste0(nodecol,"(\\d)?_\\d_"),"",desc.nodes))==
                                      max(as.numeric(gsub(paste0(nodecol,"(\\d)?_\\d_"),"",desc.nodes))))]
    } else {
      parentNum <- as.numeric(gsub(paste0(nodecol,"_0_"), "", parent.lab))
      if (parentNum %in% desc.numset){
        desc.node <- parent.lab
      } else {
        desc.node <- NULL
      }
    }
  } else if (grepl(nodecol, parent.lab)) {
    parentNum <- as.numeric(gsub(paste0(nodecol,"_0_"), "", parent.lab))
    if (parentNum %in% desc.numset){
      desc.node <- parent.lab
    } else {
      desc.node <- NULL
    }
  } else {
    desc.node <- NULL
  }
  return (desc.node)
}

##' @rdname track
##' @return A string for the node/tip
first.parent <- function (child.lab, nwktree, parent.numset, nodecol) {
  child.num <- as.numeric(gsub(".+_.+_", "", child.lab))
  if (child.num <= min(parent.numset))  return (NULL)
  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    filter(label!="") -> df
  df |> filter(grepl(paste0(child.lab,"$"),label)) |> getElement("parent") -> parent.node
  if (length(parent.node) < 1) {
    warning("The child was not in the tree.")
    return (NULL)
  }
  df |> filter(node %in% parent.node) |> getElement("label") -> parent.lab
  as.numeric(gsub(paste0(".+_0_"), "", parent.lab)) -> parent.num
  if (grepl(nodecol, parent.lab) && parent.num %in% parent.numset) {
    return (parent.lab)
  } else {
    return (first.parent(parent.lab,nwktree, parent.numset, nodecol))
  }
}

##' @rdname track
##' @return A vector of characters
all.descendants <- function (parent.lab, nwktree, desc.col) {
  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    arrange(x) |>
    filter(label!="") -> df

  df |> filter(grepl(paste0("^",parent.lab,"$"), label)) |> getElement("node") -> parent.node
  df |> filter(parent==parent.node) |> getElement("node") -> childs.node

  if (length(childs.node) > 0) {
    lapply(childs.node, function (child.node) {
      df |> filter(node==child.node) |> getElement("label") -> child.label
      all.descendants(child.label, nwktree, desc.col)
    }) |> unlist() -> res
    if (grepl(desc.col, parent.lab))  return (c(parent.lab, res))
    return (res)
  } else {
    return (NULL)
  }
}

##' @rdname track
inbetween.nodes <- function (parent.lab, child.lab, nwktree, nodecol) {
  first.parent(
    child.lab, nwktree,
    as.numeric(gsub(".+_.+_","",parent.lab)),
    gsub("_.+_.+","",parent.lab)
  ) -> ind

  if (is.null(ind)) return(NULL)

  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    arrange(x) |>
    filter(label!="") -> df

  df |> filter(grepl(nodecol,label)) |> getElement("label") -> parent.labs
  c(0,sort(as.numeric(gsub(".+_.+_", "", parent.labs)))) -> parent.numset

  df |> filter(grepl(parent.lab, label)) |> getElement("node") -> start.node
  df |> filter(grepl(child.lab, label)) |> getElement("node") -> end.node

  parent.node <- Inf
  lab.seq <- NULL
  # if (grepl(nodecol, child.lab)) lab.seq <- child.lab
  while (any(parent.node > start.node)) {
    parent.lab <- first.parent(child.lab, nwktree, parent.numset, nodecol)
    if (is.null(parent.lab))  break
    parent.node <- df |> filter(grepl(paste0("^",parent.lab,"$"), label)) |> getElement("node")
    parent.node <- parent.node[parent.node>start.node]
    if (length(parent.node) < 1)  break
    lab <- df |> filter(node==parent.node) |> getElement("label")
    lab.seq <- c(lab,lab.seq)
    child.lab <- parent.lab
  }

  return (lab.seq)
}

##' @rdname track
##' @return An vector of characters
notate.branches <- function (coal.lab, nwktree) {
  ReLabs <- unlist(str_extract_all(nwktree, "#Ha\\d_\\d+_\\d+"))
  read.tree(text=nwktree) |>
    drop_node(ReLabs) |>
    fortify(ladderize=TRUE) |>
    arrange(x) |>
    filter(label!="") -> df

  df |> filter(grepl(paste0("^",coal.lab,"$"),label)) |> getElement("node") -> coal.node
  df |> filter(parent==coal.node) |> getElement("node") -> child.nodes
  df |> filter(node %in% child.nodes & grepl("g", label)) -> tmp

  if (nrow(tmp) < 1) {
    child1 <- df |>
      filter(node %in% child.nodes & grepl("r", label)) |>
      getElement("label") |> getElement(1)
  } else {
    child1 <- tmp |> getElement("label") |> getElement(1)
  }

  df |> filter(node %in% child.nodes & label!=child1) |>
    getElement("label") |> getElement(1) -> child2

  lapply(c(child1, child2), function (child) {
    if (grepl("g", child)) {
      br.state <- notate.branches(child, nwktree)
      if (br.state) return (TRUE)
      inbetween.nodes(coal.lab, child, nwktree, "#Ha") -> ReLab
      if (!is.null(ReLab))  return (TRUE)
      return (FALSE)
    } else {
      inbetween.nodes(coal.lab, child, nwktree, "#Ha") -> ReLab
      if (!is.null(ReLab))  return (TRUE)
      return (FALSE)
    }
  }) |> unlist() |>
    (\(x) {
      all(x)
    })()
}

##' @rdname track
##' @return logical
invisible.reassortments <- function (Re.lab, nwktree) {
  coal.labs <- unlist(str_extract_all(nwktree,"g_\\d+_\\d+"))
  coal.numset <- sort(as.numeric(gsub("g_\\d_","",coal.labs)))
  first.coal.desc <- first.descendant(Re.lab, nwktree, coal.numset, "g")
  if (!is.null(first.coal.desc) && notate.branches(first.coal.desc,nwktree))
    return (TRUE)
  return (FALSE)
}

##' @rdname track
##' @return An string
get.mrca <- function (node.lab1, node.lab2, nwktree) {
  root.lab <- "m_0_0"
  do.call(intersect,lapply(c(node.lab1, node.lab2), function (node) {
    inbetween.nodes(root.lab, node, nwktree, "(m|g)")
  })) -> common.anc
  return (common.anc[order(as.numeric(gsub("g_\\d_","",common.anc)),decreasing=TRUE)][1])
}

##' @rdname track
##' @return A string
track.theother <- function (nodelab, nwktree, nodecol) {
  read.tree(text=nwktree) |>
    fortify(ladderize=TRUE) |>
    arrange(x) |>
    filter(label!="") -> df

  CoalLabs <- unlist(str_extract_all(nwktree, "g_\\d_\\d+"))
  CoalNums <- sort(as.numeric(gsub("g_\\d_", "", CoalLabs)))
  parlab <- first.parent(nodelab, nwktree, CoalNums, "g")
  parnode <- df |> filter(grepl(paste0("^",parlab,"$"), label)) |> getElement("node")
  df |> filter(parent==parnode) |>
    filter(!grepl(paste0("^",nodelab,"$"),label)) |>
    filter(grepl(nodecol,label)) |>
    getElement("label")
}

##' @rdname track
##' @return An vector of characters
fixed.reassortments <- function (parent.lab, nwktree, incongr.nodes, reversed.nodes, fixed.renodes=NULL) {
  # check whether the parent.lab is in the tree
  if (str_extract_all(nwktree, parent.lab) |> unlist() |> length() < 1)
    stop("The node label is not in the tree.")

  # 1. get the first incongruent / reversed node descendant, node.lab
  cand.incongr.reversed.set <- c(incongr.nodes, reversed.nodes)
  cand.incongr.reversed.nums <- sort(as.numeric(gsub("g_\\d_","",cand.incongr.reversed.set)))
  cand.incongr.reversed.set <- paste0("g_0_",cand.incongr.reversed.nums)
  parent.num <- as.numeric(gsub("g_\\d_","",parent.lab))
  desc.coal.nodes <- all.descendants(parent.lab, nwktree, "g")
  cand.coal.nodes <- desc.coal.nodes[desc.coal.nodes %in% cand.incongr.reversed.set]
  cand.coal.nums <- as.numeric(gsub("g_\\d_","",cand.coal.nodes))
  node.lab <- cand.coal.nodes[min(which(cand.coal.nums >= parent.num))]
  node.num <- as.numeric(gsub("g_\\d_","",node.lab))

  if (node.num < max(cand.incongr.reversed.nums)) {
    # 2. for node.lab, get all incongr/reversed descendants
    #    and all reassortment node descendants
    ddesc.coal.nodes <- all.descendants(node.lab, nwktree, "g")
    ccand.coal.nodes <- ddesc.coal.nodes[ddesc.coal.nodes %in% cand.incongr.reversed.set]
    # ccand.coal.nodes <- ccand.coal.nodes[ccand.coal.nodes != node.lab]
    desc.renodes <- all.descendants(node.lab, nwktree, "#Ha\\d")
    lapply(ccand.coal.nodes[ccand.coal.nodes != node.lab], function (node) {
      fixed.reassortments(node, nwktree, incongr.nodes, reversed.nodes)
    }) |> unlist() |> unique() -> pre.fixed.renodes
    desc.cands <- ccand.coal.nodes[ccand.coal.nodes != node.lab]
    lapply(desc.cands, function (node) {
      set <- all.descendants(node, nwktree, "g")
      set <- set[set %in% cand.incongr.reversed.set]
      set.nums <- as.numeric(gsub("g_0_","",set))
      lapply(desc.cands, function (nnode) {
        nnode.num <- as.numeric(gsub("g_0_","",nnode))
        if (nnode.num %in% set.nums && nnode.num != min(set.nums))
          return (nnode)
      }) |> unlist()
    }) |> unlist() -> ignored.nodes
    if (!node.lab %in% reversed.nodes && length(desc.renodes) <= length(ccand.coal.nodes)) {
      lapply(desc.cands[!desc.cands %in% ignored.nodes], function (node) {
        all.descendants(node, nwktree, "#Ha\\d")
      }) |> unlist() |> unique() -> pre.renodes
      fixed <- setdiff(desc.renodes, pre.renodes)
      return (c(pre.fixed.renodes, fixed))
    }
    return (pre.fixed.renodes)
  } else {
    # the last incongruent / reversed node
    desc.renodes <- all.descendants(node.lab, nwktree, "#Ha\\d")
    if ((length(desc.renodes) <= 1) && (!parent.lab %in% reversed.nodes))
      return (c(fixed.renodes,desc.renodes))
    return (NULL)
  }
}

##' @rdname track
##' @return An vector of characters
tree.fixed <- function (nwktree, incongr.nodes, reversed.nodes) {
  # get the root of the tree
  br.nodes <- str_extract_all(nwktree,"g_\\d_\\d+") |> unlist()
  br.nums <- sort(as.numeric(gsub("g_0_","",br.nodes)))
  root.lab <- paste0("g_0_",min(br.nums))

  return (fixed.reassortments(root.lab, nwktree, incongr.nodes, reversed.nodes))
}

