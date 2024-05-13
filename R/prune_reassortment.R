## Function to get the right number of reassortments
## and the reassortment times
## updated Jul 18, 2023
## 1. change the removal of extra reassortments with reversed coal nodes
drop_before_root <- function (init.nwk, seg.nwk) {
  nwks <- c(init.nwk, seg.nwk)
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |>
      fortify(ladderize=TRUE) |>
      arrange(x) |>
      filter(grepl("g_", label)) |>
      getElement("x") |>
      min()
  }) |> unlist() |> max() -> latest.root.time
  
  lapply(nwks, function (nwk) {
    #### reassortment node labels before latest root time
    read.tree(text=nwk) |>
      fortify(ladderize=TRUE) |>
      arrange(x) |>
      filter((x < latest.root.time) & (abs(x-latest.root.time) > 1e-5) & grepl("#Ha", label)) |>
      getElement("label") -> ReLabs
    
    read.tree(text=nwk) |> drop_node(ReLabs) |> write.tree()
  }) |> unlist()
}

keep_shared_re <- function (init.nwk, seg.nwk) {
  nwks <- c(init.nwk, seg.nwk)
  
  ### (0) drop #Ha2 on seg tree and corresponding #Ha2 on init infect tree
  unrelated.re <- unlist(str_extract_all(nwks[2], "#Ha2_\\d+_\\d+"))
  unrelated.nums <- sort(as.numeric(gsub("#Ha2_\\d+_","",unrelated.re)))
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |>
      drop_node(paste0("#Ha\\d_\\d+_", unrelated.nums)) |>
      write.tree()
  }) |> unlist() -> nwks
  
  ### (1) get the labels of common reassortments, #Ha2 + #Ha1
  lapply(seq_along(nwks), function (i) {
    nwk <- nwks[i]
    pattern1 <- "#Ha\\d_\\d+_\\d+"
    pattern2 <- "#Ha\\d_\\d+_"
    str_extract_all(nwk, pattern1) |> 
      unlist() |>
      (\(x) {gsub(pattern2,"",x)})() |>
      as.numeric() |>
      sort()
  }) |> unlist() |>
    table() -> Re.count.table
  
  # keep #Ha2 in init infect tree, and only #Ha1 in segment tree
  kept.ReNum <- as.numeric(names(Re.count.table)[Re.count.table>1])
  
  ### (2) drop the unshared reassortment nodes
  lapply(seq_along(nwks), function (i) {
    nwk <- nwks[i]
    str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+") |> 
      unlist() |>
      (\(x) {gsub("#Ha\\d+_\\d+_","",x)})() |>
      as.numeric() |> sort() -> all.nodenum
    
    read.tree(text=nwk) |> 
      drop_node(paste0("#Ha\\d_0_",all.nodenum[!all.nodenum%in%kept.ReNum])) |>
      write.tree()
  }) |> unlist()
}

drop_parent_re <- function (init.nwk, seg.nwk) {
  nwks <- c(init.nwk, seg.nwk)
  ## Scenario 1:  #Ha ---- #Ha
  lapply(nwks, function (nwk) {
    lapply(unlist(str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+")), function (ReLab) {
      read.tree(text=nwk) |>
        fortify(ladderize=TRUE) |>
        arrange(x) -> df
      
      df |> filter(label==ReLab) |> getElement("parent") -> parnode
      df |> filter(node==parnode) |> getElement("label") -> parlab
      
      if (grepl("#", parlab)) return (gsub("#Ha\\d+_\\d+_", "", parlab))
      return (NULL)
    }) |> unlist() |>
      as.numeric() 
  }) |> unlist() |> unique() -> parReNum1
  
  # drop nodes
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |> 
      drop_node(paste0("#Ha\\d_\\d+_",parReNum1)) |>
      write.tree() -> nwk
  }) |> unlist() -> nwks
  
  ## Scenario 2:
  ##               |----#Ha-------
  ##    #Ha --> g -|
  ##               |--------#Ha---
  lapply(nwks, function (nwk) {
    lapply(unlist(str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+")), function (ReLab) {
      if (invisible.reassortments(ReLab, nwk))
        return (as.numeric(gsub("#Ha\\d_\\d+_","",ReLab)))
      return (NULL)
    }) |> unlist()
  }) |> unlist() |> unique() -> parReNum2
  
  # drop nodes
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |> 
      drop_node(paste0("#Ha\\d_\\d+_",parReNum2)) |>
      write.tree()
  }) |> unlist()
}

get_reversions <- function (init.nwk, seg.nwk) {
  nwks <- c(init.nwk, seg.nwk)
  lapply(nwks, function (nwk) {
    CoalLabs <- str_extract_all(nwk, "g_\\d+_\\d+") |> unlist()
    parnums <- as.numeric(gsub("g_0_","",CoalLabs))
    relabs <- str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+") |> unlist()
    seg.unReNums <- sort(as.numeric(gsub("#Ha\\d+_\\d+_", "", relabs[grep("#Ha\\d+", relabs)])))
    common.set <- NULL
    for (renum in rev(seg.unReNums)) {
      ReLab <- paste0("#Ha\\d+_\\d+_", renum)
      first.parent(ReLab, nwk, parnums, "g") -> parlab
      if (!is.null(parlab)) {
        parnums <- setdiff(parnums, as.numeric(gsub("g_0_","",parlab)))
        common.set <- c(common.set, renum, parlab)
      }
    }
    return (common.set)
  }) |> unlist() -> reversionRe
  
  if (!is.null(reversionRe)) {
    reversionRe |> 
      matrix(ncol=2, byrow=TRUE) |>
      as.data.frame() |>
      arrange(V1) |>
      group_by(V2) |>
      filter(n()>1) |>
      (\(x) {
        colnames(x) <- c("ReLab", "CoalLab")
        return (x)
      })() -> reversionRe
  }
  
  return (reversionRe)
}

compute_preRRtab <- function (init.nwk, seg.nwk, incongr.init, incongr.seg, reversions) {
  trees <- c(init.nwk, seg.nwk)
  
  reversed.set <- unique(reversions$CoalLab)[order(as.numeric(gsub("g_0_","",unique(reversions$CoalLab))))]
  reversed.Numset <- sort(as.numeric(gsub("g_0_","",reversed.set)))
  
  read.tree(text=trees[1])$node.lab -> labs
  
  ## if the # of the par set < # of unre set
  ## remove the first few unsampled reassortment nodes
  reNums <- sort(as.numeric(gsub("#Ha2_\\d+_", "", labs[grep("#Ha2", labs)])))
  
  ## the rejoined nodes
  incongr.set <- incongr.seg
  incongr.Numset <- sort(as.numeric(gsub("g_0_","",incongr.set)))
  seg.reNums <- reNums
  seg.set <- NULL
  seg.parnums <- sort(c(incongr.Numset, reversed.Numset))
  candidates <- paste0("g_0_", seg.parnums)
  for (dispar in rev(candidates)) {
    last.descendant(dispar, trees[2], seg.reNums, "#Ha1") -> corres.re
    if (is.null(corres.re))
      first.parent(dispar, trees[2], seg.reNums, "#Ha1") -> corres.re
    # print(c(dispar, corres.re))
    if (is.null(corres.re)) next
    seg.reNums <- setdiff(seg.reNums, as.numeric(gsub("#Ha1_0_","",corres.re)))
    seg.set <- c(seg.set, gsub("#Ha1_0_", "", corres.re), dispar)
  }
  matrix(seg.set, ncol=2, byrow=TRUE) |>
    as.data.frame() |>
    (\(x) {
      colnames(x) <- c("ReLab", "RejoinNode")
      return (x)
    })() |>
    arrange(ReLab) -> seg.join
  
  # the removed nodes
  incongr.set <- incongr.init
  incongr.Numset <- sort(as.numeric(gsub("g_0_","",incongr.set)))
  seg.set <- NULL
  seg.parnums <- sort(c(incongr.Numset, reversed.Numset))
  for (unre in sort(as.numeric(seg.join$ReLab))) {
    first.parent(paste0("#Ha2_0_", unre), trees[1], seg.parnums, "(g)|(#Ha2)|(m)") -> corres.par
    seg.parnums <- setdiff(seg.parnums, as.numeric(gsub("((g)|(#Ha2))_0_", "", corres.par)))
    seg.set <- c(seg.set, unre, corres.par)
  }
  matrix(seg.set, ncol=2, byrow=TRUE) |>
    as.data.frame() |>
    (\(x) {
      colnames(x) <- c("ReLab", "RemoveNode")
      return (x)
    })() -> seg.leave
  
  merge(seg.leave, seg.join, all=TRUE, by="ReLab",incomparables="") |>
    arrange(as.numeric(ReLab)) |> 
    filter(RemoveNode!=RejoinNode)
}

refine_RRtab <- function (init.nwk, seg.nwk, incongr.init, incongr.seg, RRtab) {
  
  intersect(RRtab$RemoveNode, RRtab$RejoinNode) |> 
    (\(x) {
      x[order(as.numeric(gsub("g_\\d+_","",x)))]
    })() -> RevNodes
  
  trees <- c(init.nwk, seg.nwk)
  
  lapply(trees, function (nwk) {
    rev(sort(as.numeric(gsub("g_\\d_", "", unlist(str_extract_all(nwk, "g_\\d_\\d+"))))))
  }) -> CoalNums
  
  ## Condition: some incongruent nodes can not be interpreted
  ## if the corresponding reassortment nodes are removed
  incongr.nodes <- c(incongr.init, incongr.seg)
  incongr.nodes <- incongr.nodes[order(as.numeric(gsub("g_0_","",incongr.nodes)))]
  renums.init <- as.numeric(RRtab$ReLab)
  renums.seg <- renums.init
  ReCoalnodes.init <- RRtab$RemoveNode
  ReCoalnodes.seg <- RRtab$RejoinNode
  lapply(trees, function (tree) {
    tree.fixed(tree, incongr.nodes, RevNodes)
  }) |> unlist() |>
    (\(x) {
      as.numeric(gsub("#Ha\\d_\\d_","",x))
    })() |> unique() |> sort() -> fixed.renodes
  
  lapply(RevNodes, function (node) {
    RRtab |> filter((RemoveNode==node) | (RejoinNode==node)) -> Rev
    
    ## Format
    #    Remove Rejoin
    # m       b      a
    # n       a      c
    a <- node
    b <- Rev$RemoveNode[Rev$RemoveNode!=node]
    c <- Rev$RejoinNode[Rev$RejoinNode!=node]
    m <- as.numeric(Rev$ReLab[Rev$RemoveNode==b])
    n <- as.numeric(Rev$ReLab[Rev$RejoinNode==c])
    
    ## Scenario 1: in init infect tree
    ## and remove m
    if (!m %in% fixed.renodes) {
      get.mrca(paste0("#Ha\\d_\\d_",m), paste0("#Ha\\d_\\d_",n), trees[1]) -> mrca.init
      node.mrca.desc.init <- first.parent(
        mrca.init, trees[1], 
        as.numeric(gsub("g_\\d_","",node)), "g")
      
      if (!is.null(node.mrca.desc.init) && mrca.init==b) return(m)
      
      btwn.nodes <- inbetween.nodes(a,paste0("#Ha2_0_",m),trees[1],"g")
      if (all(lapply(btwn.nodes, function (node) {
        notate.branches(node, trees[1])
      }) |> unlist()) && !is.null(btwn.nodes)) 
        return(m)
    }
    
    ## Scenario 2: in the seg tree
    ## and remove n
    if (!n %in% fixed.renodes) {
      get.mrca(paste0("#Ha\\d_\\d_",m), paste0("#Ha\\d_\\d_",n), trees[2]) -> mrca.seg
      node.mrca.desc.seg <- first.parent(
        mrca.seg, trees[2], 
        as.numeric(gsub("g_\\d_","",node)), "g")
      
      if (!is.null(mrca.seg) && mrca.seg %in% unname(unlist(RRtab[,3])) && 
          !is.null(node.mrca.desc.seg))   return(n)
      
      AllReLabs <- unlist(str_extract_all(trees[2],"#Ha\\d_\\d_\\d+"))
      AllReNums <- sort(as.numeric(gsub("#Ha\\d_\\d_","",AllReLabs)))
      c.par.seg.re <- first.parent(c,trees[2],c(CoalNums[[2]],AllReNums),"#Ha1")
      
      if (!is.null(c.par.seg.re) && c.par.seg.re==m) return(n)
    }
    
    ## Scenario 3: in both trees, neither m nor n
    ## but descendants of the reversed node
    ReDescs.init.num <- as.numeric(gsub("#Ha\\d_\\d_","",all.descendants(a,trees[1],"#Ha2")))
    ReDescs.seg.num <- as.numeric(gsub("#Ha\\d_\\d_","",all.descendants(a,trees[2],"#Ha1")))
    CommonReDesc.num <- sort(ReDescs.init.num[ReDescs.init.num %in% ReDescs.seg.num])
    CommonReDesc.num <- CommonReDesc.num[!CommonReDesc.num %in% fixed.renodes]
    if (length(CommonReDesc.num) < 1) return (NULL)
    lapply(CommonReDesc.num, function (renode.num) {
      tryCatch({
        par1.init <- inbetween.nodes(a,paste0("#Ha2_0_", renode.num),trees[1],"g")[1]
        descs.init.par1.br <- all.descendants(par1.init,trees[1],"g")
        descs.init.par1.br.incongr <- descs.init.par1.br[descs.init.par1.br %in% c(incongr.init, incongr.seg, RevNodes)]
        descs.init.par1.re <- all.descendants(par1.init,trees[1],"#Ha2")
        btwn.init.re <- inbetween.nodes(par1.init,paste0("#Ha2_0_", renode.num),trees[1],"#Ha2")
        par1.seg <- inbetween.nodes(a,paste0("#Ha1_0_", renode.num),trees[2],"g")[1]
        descs.seg.par1.br <- all.descendants(par1.seg,trees[2],"g")
        descs.seg.par1.br.incongr <- descs.seg.par1.br[descs.seg.par1.br %in% c(incongr.init, incongr.seg, RevNodes)]
        descs.seg.par1.re <- all.descendants(par1.seg,trees[2],"#Ha1")
        if (length(descs.init.par1.re) > length(descs.init.par1.br.incongr) && 
            length(btwn.init.re) < 1 &&
            length(descs.seg.par1.re) > length(descs.seg.par1.br.incongr))
          return (renode.num)
      },error = function(e) {return(NULL)})
    }) |> unlist() -> remove.cands
    if (length(remove.cands) > 0) return (min(remove.cands))
    return (NULL)
  }) |> unlist() |> unique() -> Re.removed
  
  if (is.null(Re.removed)) {
    apply(RRtab, 1, function (l) {
      ## scenario 1: wrong assignment at the first place
      re.par.init <- first.parent(paste0("#Ha\\d_\\d_",l[1]), trees[1], CoalNums[[1]], "g")
      re.par.seg <- first.parent(paste0("#Ha\\d_\\d_",l[1]), trees[2], CoalNums[[2]], "g")
      
      if (!all(re.par.init %in% c(incongr.init, RevNodes), re.par.seg %in% c(incongr.seg, RevNodes)))
        return (as.numeric(l[1]))
      
      ## scenario 2
      remove.mark <- notate.branches(l[2], trees[1])
      rejoin.mark <- notate.branches(l[3], trees[2])
      
      if (!all(remove.mark,rejoin.mark))  return (NULL)
      
      first.parent(
        unname(unlist(l[2])),trees[1],
        CoalNums[[1]],"g"
      ) -> par.init
      first.parent(
        unlist(unname(l[3])),trees[2],
        CoalNums[[2]],"g"
      ) -> par.seg
      
      if (is.null(par.seg) || is.null(par.init))  return (NULL)
      
      first.parent(
        par.seg, trees[2],
        as.numeric(gsub("g_0_","",par.init)), "g"
      ) -> par.seg.par
      
      if (is.null(par.seg.par)) return (NULL)
      
      if (par.init %in% RevNodes) {
        btwn.nodes <- c(inbetween.nodes(par.init, par.seg, trees[2], "g"), par.seg)
        lapply(btwn.nodes, function (node) {
          notate.branches(node, trees[2])
        }) |> unlist() -> ind
        if (ind)  return(as.numeric(unname(l[1])))
      }
    }) |> unlist() -> Re.removed
  }
  
  if (length(Re.removed) < 1) return (RRtab)
  
  lapply(trees, function (nwk) {
    read.tree(text=nwk) |>
      drop_node(paste0("#Ha\\d_\\d+_",Re.removed)) |>
      write.tree()
  }) |> unlist() -> new.trees
  
  return (func_LeaveJoin(new.trees[1], new.trees[2]))
}

func_LeaveJoin <- function (init.nwk, seg.nwk) {
  
  ## ------------------------------- ##
  #   Part I: Data Preprocessing      #
  ## ------------------------------- ##
  nwks <- c(init.nwk, seg.nwk)
  str_extract_all(nwks[1], "g_0_\\d+") |> unlist() -> Init.CoalLabs
  str_extract_all(nwks[2], "g_0_\\d+") |> unlist() -> Seg.CoalLabs
  
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |>
      fortify(ladderize=TRUE) |>
      arrange(x) |>
      filter(grepl("g_", label)) |>
      getElement("label")
  }) |> unlist() |>
    table() -> coallab.count.table
  
  CommonCoalLabs <- names(coallab.count.table)[coallab.count.table>1]
  IncongrCoalLabs <- names(coallab.count.table)[coallab.count.table<2]
  lapply(nwks, function (nwk) {
    read.tree(text=nwk) |>
      fortify(ladderize=TRUE) |>
      arrange(x) |>
      filter(label %in% IncongrCoalLabs) |> 
      getElement("label")
  }) -> tmp
  IncongrInit <- tmp[[1]]
  IncongrSeg <- tmp[[2]]
  
  if (length(IncongrCoalLabs) < 1)  return (NULL)
  
  ## ------------------------------------------ ##
  #   Part II: Invisible Reassortment Removal    #
  ## ------------------------------------------ ##
  ## Step 1. Drop reassortments before the most earliest root time
  drop_before_root(nwks[1], nwks[2]) -> seg.nwks.init
  
  ## Step 2. Keep shared reassortments between init infect trees and segment tree
  keep_shared_re(seg.nwks.init[1], seg.nwks.init[2]) -> seg.nwks.init_2
  
  ## Step 3. drop reassortment nodes (#Ha1) in segment tree
  drop_parent_re(seg.nwks.init_2[1], seg.nwks.init_2[2]) -> seg.nwks.init_3
  
  ## Step 4. Find the reversion coalescent nodes by reassortments
  get_reversions(seg.nwks.init_3[1], seg.nwks.init_3[2]) -> reversionRe
  
  ## Step 5. Compute the preliminary Remove-Rejoin table
  compute_preRRtab(
    seg.nwks.init_3[1],  seg.nwks.init_3[2], 
    IncongrInit, IncongrSeg, reversionRe
  ) -> RRtab
  
  if (sum(!(unique(unlist(unname(RRtab[,-1]))) %in% IncongrCoalLabs)) < 1) return (RRtab)
  
  lapply(seg.nwks.init_3, function (nwk) {
    str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+") |>
      unlist() |>
      (\(x) {
        gsub("#Ha\\d+_\\d+_", "", x)
      })() |>
      as.numeric() |>
      sort() -> res
    remove <- res[!res %in% as.numeric(RRtab$ReLab)]
    read.tree(text=nwk) |>
      drop_node(paste0("#Ha\\d_\\d+_",remove)) |>
      write.tree()
  }) |> unlist() -> seg.nwks.init_4
  
  ## Step 6. Refine the Remove-Rejoin table by removing fake reversions
  refine_RRtab(
    seg.nwks.init_4[1],  seg.nwks.init_4[2], 
    IncongrInit, IncongrSeg, RRtab
  ) -> RRtab
  
  return (RRtab)
}



##' dismiss all invisible reassortments given the first-infection tree and segment trees
Pruned_Re <- function (gp) {
  init.nwk <- getInfo(gp,prune=TRUE,hide=FALSE,tree=TRUE)$tree[3]
  seg.nwks <- getInfo(gp,prune=TRUE,hide=FALSE,tree=TRUE)$tree[1:2]
  
  # get the prune-and-regraph (leave-and-rejoin) data frame
  lapply(seg.nwks, function (seg.nwk) {
    func_LeaveJoin(init.nwk, seg.nwk)
  }) -> Reassortments
  
  # Get the final nwk
  lapply(Reassortments, function (list) {
    return (list$ReLab)
  }) |> unlist() |>
    as.numeric() |>
    sort() -> AllReNums
  
  all.nwks <- c(init.nwk, seg.nwks)
  
  lapply(seq_along(all.nwks), function (i) {
    nwk <- all.nwks[i]
    str_extract_all(nwk, "#Ha\\d+_\\d+_\\d+") |>
      unlist() |>
      (\(x) {
        gsub("#Ha\\d+_\\d+_", "", x)
      })() |>
      as.numeric() |>
      sort() -> res
    remove <- res[!res %in% AllReNums]
    
    if (i == 1) {
      read.tree(text=nwk) |>
        drop_node(paste0("#Ha2_\\d+_", remove)) |>
        write.tree() -> output
    } else {
      read.tree(text=nwk) |>
        drop_node(paste0("#Ha2_\\d+_", res)) |>
        drop_node(paste0("#Ha1_\\d+_", remove)) |>
        write.tree() -> output
    }
    return (output)
  }) |> unlist() -> nwksfinal
  
  return(list(nwks=nwksfinal, reassortments=Reassortments))
}

