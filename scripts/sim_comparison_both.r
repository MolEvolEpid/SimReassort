##########################################################################
##  infer reassortments using CoalRe and TreeKnit, Version 2024-01-18   ##
##########################################################################

#------------------------------------#
# Scenario 2: Both segment reassort  #
#------------------------------------#

## Procedure Summary:
### 1. simulate sets of trees, fixed parameters, constant expected pop size (lambda = mu + psi)
### 2. simulate alignments: seg A, B: len=1e3, FIT: len=1e4, evol rate=1e-5
### 3. One inference schemes
#### 3.1 only seg trees imported
#### 3.2 seg trees + FIT imported
####     - pairwise comparison in TreeKnit
####     - all imported to CoalRe, fixing tree structure of FIT approximately

rm(list=ls())
path <- "~/your_path/"
if (!dir.exists(path))	dir.create(path)
setwd(path)
library(SimReassort)

args <- commandArgs(trailingOnly=TRUE)
re.num <- args[1]
re.dir <- paste0("renum",re.num)
current.path <- paste0(path,re.dir,"/")
if (!dir.exists(current.path)) dir.create(current.path)
setwd(current.path)

nsamp <- 50
nrep <- 100
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
print(paste0("Number of cpus used: ",n.cpus))

params <- c(lambda=5,mu=4.75,psi=.25,rhoA=.12,rhoB=.08,rhoAB=0,n0=1)

registerDoParallel(cores=n.cpus)
foreach(task.no=seq_len(nrep), .combine=rbind, .packages=c("SimReassort","dplyr","stringr","ape"),
        .errorhandling="remove") %dopar% {
        
  ##################################################################
  ##  1. simulate a set of trees that has `re.num` reassortments  ##
  ##################################################################
	
	while (TRUE) {
	  gc()
	  t <- .01
	  y <- simulate("LBDPwr",lambda=200,mu=0,psi=0,rhoA=0,rhoB=0,rhoAB=0,
	                n0=1,time=t,t0=0,cont=FALSE)
	  
	  while (current_popsize(y) < 100) {
	    t <- t + .00001
	    y <- y |> simulate(time=t)
	  }
	  
	  # current_popsize(y)
	  if (current_popsize(y)!=100)  next
	  
	  t100 <- t
	  # burn-in period, constant popsize without sample. The root is younger than t100
	  withTimeout({
	    while (current_popsize(y) > 0 && get_roottime(y, init.include=TRUE) <= t100 && current_popsize(y) < 200) {
	      t <- t + .1
	      y |> simulate(lambda=params[1], mu=params[1], psi=0, 
	                    rhoA=params[4], rhoB=params[5], rhoAB=params[6],
	                    time=t, cont=FALSE) -> y
	    }
	  }, timeout=600, onTimeout="warning") -> res
	  
	  if (!is.null(res) || current_popsize(y) < 1 ||  get_roottime(y, init.include=TRUE) <= t100)  next
	  
	  # getting `nsamp` samples
	  withTimeout({
	    while (current_popsize(y) > 0) {
	      t <- t + .001
	      y <- y |> simulate(lambda=params[1],mu=params[2],psi=params[3],
	                         rhoA=params[4],rhoB=params[5],rhoAB=params[6],
	                         time=t,cont=FALSE)
	      tr <- read.tree(text=getInfo(y,prune=TRUE,tree=TRUE)$tree[1])
	      if (length(tr$tip.label) >= nsamp + 2) break
	    }
	  }, timeout=60, onTimeout="warning") -> res
	  
	  if (!is.null(res) || (length(tr$tip.label) != nsamp + 2) || current_popsize(y) < 1) next
	  
	  Pruned_Re(y) -> pruned
	  Reassorts <- pruned$reassortments
	  trees <- pruned$nwks
	  reassorts <- unlist(str_extract_all(trees[1], "#Ha2_\\d+_\\d+"))
	  
	  if (length(reassorts) == re.num)  break
	}
          
  dir.name <- paste0("run",task.no)
  dir.path <- paste0(current.path,dir.name,"/")
  
          
	if (!dir.exists(dir.path))  dir.create(dir.path)
  
  setwd(dir.path)
  
  saveRDS(list(params,y), file=paste0(dir.path,"gpsim_output",task.no,".rds"))
  
  lapply(y |> getInfo(hide=TRUE,tree=TRUE) |> getElement("tree"), function(tr) {
    read.tree(text=tr) |>
      fortify(ladderize=TRUE) |>
      filter(!grepl("(i)|(m)", label)) |>
      arrange(x) |>
      mutate(label=if_else(!grepl("r",label),"",label)) |>
      group_by(isTip) |>
      mutate(label=if_else(isTip, paste0("seq",row_number()), "")) |>
      ungroup()
  }) -> treedfs
  
  tipdates <- treedfs[[1]] |> filter(isTip) |> select(label, x)
  write.table(as.matrix(tipdates), paste0(dir.path,"tipdates",task.no,".txt"),
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  evol.rate <- 5e-3
  lapply(treedfs[1:3], function(df) {
    df |> write_nwk() -> nwk
    gsub("^(\\()(.+)(\\))(;)$","\\2\\4", nwk) -> nwk
    tmp <- read.tree(text=nwk)
    tmp$edge.length <- tmp$edge.length*evol.rate
    write.tree(tmp)
  }) |> unlist() -> trees_nwk
      

	##############################
  ##  2. simulate alignments  ##
  ##############################

  nseq <- 1e3
  lapply(1:3, function (k) {
    tr <- trees_nwk[k]
    if (k > 2)  nseq <- 1e4
    write.table(trees_nwk[k],file=paste0(dir.path,"tre",task.no,"_",k,".nwk"),col.names=FALSE,row.names=FALSE,quote=FALSE)
    system(paste0("iqtree2 --alisim ", dir.path, "tre",task.no,"_",k," -m JC -t ", dir.path, "tre",task.no,"_",k,".nwk -af fasta -seed 1234 --length ",nseq), intern=TRUE) -> fasta.out
  })

  	
	##############################
  ##  3. Infer reassortments  ##
  ##############################

  #------------------------------------------#
  #  3.1 Scheme 1: only using segment trees  #
  #------------------------------------------#

  doc <- read_xml("template.xml", package="xml2")
  for (k in seq_along(trees_nwk)) {
    read.fasta(paste0(dir.path,"tre",task.no,"_",k,".fa")) -> fa.file
    # building the xml
    data <- xml_find_first(doc, paste0("//data[@id='tre",k,"']"))
    lapply(seq_len(nrow(tipdates)), function(j) {
      xml_add_child(.x=data, .value="sequence",
                    id=paste0("seq_",k,"_", tipdates$label[j]),
                    spec="Sequence", taxon=tipdates$label[j], totalcount="4", value=toupper(paste0(fa.file[[tipdates$label[j]]],collapse="")))
    })
    
    trait <- xml_find_first(doc, paste0("//trait[@id='dateTrait.t:tre",k,"']"))
    xml_attr(trait, "value") <- paste0(apply(tipdates,1,function(td) paste0(td,collapse="=")), collapse=",")
  }
  
  write_xml(doc, paste0(dir.path,"simulator",task.no,".xml"))
  
  system(paste0("beast -beagle_sse -threads 4 -overwrite ",dir.path, "simulator",task.no,".xml"))
  
  system(paste0("applauncher ReassortmentNetworkSummarizer ", dir.path, "simulator",task.no,".network.trees ", dir.path, "simulator",task.no,".mcc.network.tree"))
  
  system(paste0("loganalyser ",dir.path, "tre2.log"), intern=TRUE) -> logoutput
  write.table(logoutput, paste0(dir.path, "logoutput", task.no, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  lines <- c(
    "using TreeTools", "using TreeKnit", "",
    paste0("cd(\"/target_dir/",re.dir,"/", dir.name, "/\")"), "",
    paste0("t1 = read_tree(\"tre", task.no, "_1.nwk\", label=\"tree1\")"),
    paste0("t2 = read_tree(\"tre", task.no, "_2.nwk\", label=\"tree2\")"),
    "treelist = Any[t1, t2]",
    "mcc = naive_mccs(treelist)",
    "",
    "MCCs = run_treeknit!(t1, t2, OptArgs(;pre_resolve=true, resolve=true, strict=true));",
    "",
    "arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t1, MCCs.mccs[Set([\"tree1\", \"tree2\"])]);",
    "",
    paste0("write(\"arg", task.no, ".nwk\", arg)"),
    paste0("write_mccs(\"MCCs", task.no, ".json\", MCCs)"),
    "",
    "exit()"
  )
  
  writeLines(lines, paste0(dir.path,"treeknit_script",task.no,".jl"))
  
  system(paste0("julia ", dir.path, "treeknit_script",task.no,".jl"))


	#-------------------------------------#
  #  3.2 Scheme 2: segment trees + FIT  #
  #-------------------------------------#

	doc <- read_xml("template_fit.xml", package="xml2")
	for (k in seq_along(trees_nwk[1:3])) {
	  read.fasta(paste0(dir.path,"tre",task.no,"_",k,".fa")) -> fa.file
	  
	  # building the xml
	  data <- xml_find_first(doc, paste0("//data[@id='tre",k,"']"))
	  lapply(seq_len(nrow(tipdates)), function(j) {
	    xml_add_child(.x=data, .value="sequence",
	                  id=paste0("seq_",k,"_", tipdates$label[j]),
	                  spec="Sequence", taxon=tipdates$label[j], totalcount="4", value=toupper(paste0(fa.file[[tipdates$label[j]]],collapse="")))
	  })
	  
	  trait <- xml_find_first(doc, paste0("//trait[@id='dateTrait.t:tre",k,"']"))
	  xml_attr(trait, "value") <- paste0(apply(tipdates,1,function(td) paste0(td,collapse="=")), collapse=",")
  }
	write_xml(doc, paste0(dir.path,"fit_simulator",task.no,".xml"))

	system(paste0("beast -beagle_sse -threads 4 -overwrite ",dir.path, "fit_simulator",task.no,".xml"))

	system(paste0("applauncher ReassortmentNetworkSummarizer ", dir.path, "fit_simulator",task.no,".network.trees ", dir.path, "fit_simulator",task.no,".mcc.network.tree"))
	
	system(paste0("loganalyser ",dir.path, "tre3.log"), intern=TRUE) -> logoutput
	write.table(logoutput, paste0(dir.path, "fit_logoutput", task.no, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)


	# TreeKnit in two steps:
  ## Step 1: input the first-infection tree and the seg A tree
	lines <- c(
	  "using TreeTools", "using TreeKnit", "",
	  paste0("cd(\"/home/qianying/cchfv/Comparison/twostrain/",re.dir,"/", dir.name, "/\")"), "",
	  paste0("t1 = read_tree(\"tre", task.no, "_3.nwk\", label=\"tree1\")"),
	  paste0("t2 = read_tree(\"tre", task.no, "_1.nwk\", label=\"tree2\")"),
	  "treelist = Any[t1, t2]",
	  "mcc = naive_mccs(treelist)",
	  "",
	  "MCCs = run_treeknit!(t1, t2, OptArgs(;pre_resolve=true, resolve=true, strict=true));",
	  "",
	  "arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t1, MCCs.mccs[Set([\"tree1\", \"tree2\"])]);",
	  "",
	  paste0("write(\"fit1_arg", task.no, ".nwk\", arg)"),
	  paste0("write_mccs(\"fit1_MCCs", task.no, ".json\", MCCs)"),
	  "",
	  "exit()"
	)
	
	writeLines(lines, paste0(dir.path,"fit1_treeknit_script",task.no,".jl"))
	
	system(paste0("julia ", dir.path, "fit1_treeknit_script",task.no,".jl"))
	
	## Step 2: input the first-infection tree and the seg B tree
	lines <- c(
	  "using TreeTools", "using TreeKnit", "",
	  paste0("cd(\"/target_dir/",re.dir,"/", dir.name, "/\")"), "",
	  paste0("t1 = read_tree(\"tre", task.no, "_3.nwk\", label=\"tree1\")"),
	  paste0("t2 = read_tree(\"tre", task.no, "_2.nwk\", label=\"tree2\")"),
	  "treelist = Any[t1, t2]",
	  "mcc = naive_mccs(treelist)",
	  "",
	  "MCCs = run_treeknit!(t1, t2, OptArgs(;pre_resolve=true, resolve=true, strict=true));",
	  "",
	  "arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t1, MCCs.mccs[Set([\"tree1\", \"tree2\"])]);",
	  "",
	  paste0("write(\"fit2_arg", task.no, ".nwk\", arg)"),
	  paste0("write_mccs(\"fit2_MCCs", task.no, ".json\", MCCs)"),
	  "",
	  "exit()"
	)
	
	writeLines(lines, paste0(dir.path,"fit2_treeknit_script",task.no,".jl"))
	
	system(paste0("julia ", dir.path, "fit2_treeknit_script",task.no,".jl"))
}
gc()
q(save="no")
