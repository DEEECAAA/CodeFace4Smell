#! /usr/bin/env Rscript
## Analyse the developer connections

## This file is part of Codeface. Codeface is free software: you can
## redistribute it and/or modify it under the terms of the GNU General Public
## License as published by the Free Software Foundation, version 2.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
## Copyright 2010, 2011 by Wolfgang Mauerer <wm@linux-kernel.net>
## Copyright 2012, 2013, Siemens AG, Wolfgang Mauerer <wolfgang.mauerer@siemens.com>
## All Rights Reserved.

## TODO: Can we deploy the page rank idea to patches as well?
## For instance, we could interpret that patches are "touched" by
## files, so patches involving more heaviliy patched files would
## receive higher scores. Or we could also determine how many
## people touched a particular piece of code (to simplify things,
## this computation could be done on the file level in a first stage),
## and use this as basis for the patch evaluation.
## NOTE: To compare the results to traditional metric, we could
## also use anonymised names to not insult anyone.
## TODO: Could we use graph cohesion as a more mature (and
## mathematically sound) replacement to OOP cohesion?
## TODO: Another idea is also categorisation by subsystems.
##       (i.e., use a known/suspected distribution into subsystems
##        for learning, and then apply the inferred clustering on
##        the complete data set to infer subsystems in a quasi-Bayesian
##        way)
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(logging))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
source("../utils.r", chdir=TRUE)
source("../config.r", chdir=TRUE)
source("../db.r", chdir=TRUE)
source("../query.r", chdir=TRUE)
source("community_metrics.r")
source("network_visualization.r")

#######################################################################
##		Lower Level Functions
#######################################################################

##========================================================================
##		Utility
##========================================================================
read.oslom <- function(input.file){
  ## reads a file created by the oslom clustering program
  con   	 <- file(input.file, open='r')
  comms 	 <- list() 
  membership <- c()
  csize      <- c()
  res		 <- list()
  while(length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
    if ( "#" == substring(line,1,1)) {
	  next
	}
	else{
      comm.str <- (strsplit(line, " "))
      comm.int <- list(as.numeric(comm.str[[1]]))
      comms    <- c(comms,comm.int)
	}
  }
  close(con)  

  ##TODO: change this to support overlapping communities, ATM nodes split
  ##			between two clusters end up arbitrarily in one
  ## create a igraph communities style object
  for (i in 1:length(comms)){
    verts <- comms[[i]]
    membership[verts] <- i
  }
  membership.conseq <- remapConsecSeq(membership)
  for (i in 1:length(unique(membership.conseq))) {
    csize[i] <- length(which(membership.conseq == i))
  }
  class(res)     <- "communities"
  res$membership <- membership.conseq 
  res$csize	     <- csize
  res$algorithm  <- "OSLOM"
  return(res)
}


## Given an eMail address like "Name N. Surname <name.surname@domain.com>",
## extract the person name without the electronic address
name.of.email <- function(str) {
  return(str_sub(str, 1, str_locate(str, "<")[1]-2))
}

## NOTE: We rely that the ids are already sorted numerically. To ensure
## this is really the case, a check could do no harm
ID.to.name <- function(.iddb, .id) {
  return(.iddb[which(.iddb$ID==.id),]$Name)
}

IDs.to.names <- function(.iddb, .ids) {
  return(sapply(.ids,
                function(.id) { .iddb[which(.iddb$ID==.id),]$Name }))
}

name.to.ID <- function(.iddb, name) {
  return(which(.iddb$Name==name))
}

rotate.label <- function(x) {
  return(paste("\\rotatebox{60}{", x, "}", sep=""))
}

rotate.label.30 <- function(x) {
  return(paste("\\rotatebox{30}{", x, "}", sep=""))
}

txt.comm.subsys <- function(.comm, .id.subsys, i) {
  idx <- which(.comm$membership==i)

  ##  return(summary(.id.subsys[idx,2:dim(.id.subsys)[2]]))
  res <- apply(.id.subsys[idx,2:dim(.id.subsys)[2]], 2, fivenum)
  rownames(res) <- c("lw", "lh", "med", "uh", "uw")
  return(res)
}

get.rank.by.field <- function(.iddb, .field, N=dim(.iddb)[1]) {
  res <- .iddb[c("ID", "Name", .field)]
  res <- res[order(res[[.field]], decreasing=TRUE),]
  s <- sum(res[,3])
  res <- cbind(res, data.frame(percent=res[,3]/s*100))
  res <- cbind(res, data.frame(norm=scale.data(res[,3], 0, 1)))

  return(res[1:N,])
}

int.to.hex <- function(n, fill=2) {
  loc <- n
  class(loc) <- "hexmode"
  loc <- as.character(loc)
  if (nchar(loc) != fill) {
    loc <- paste(rep("0", fill - nchar(loc)), loc, sep="")
  }
  return(as.character(loc))
}

four.digit <- function(n) {
  loc <- as.character(n)

  if (nchar(loc) != 4) {
    loc <- paste(paste(rep("0", 4-nchar(loc)), collapse=''), loc, sep="")
  }

  return(loc)
}

col.to.hex <- function(prefix="0x", r,g,b) {
  return (paste(prefix, int.to.hex(r), int.to.hex(g), int.to.hex(b), sep=""))
}


## Return indices of vertices that are in the largest connected subgraph
largest.subgraph.idx <- function(graph) {
  ## Find all connected subgraphs
  g.clust <- clusters(graph)

  ## Select the id of the largest cluster, and find all mathing indices
  clusters <- data.frame(id=1:length(g.clust$csize), size=g.clust$csize)
  clusters <- clusters[sort(clusters$size, index.return=TRUE, decreasing=TRUE)$ix,]
  id.largest <- clusters$id[1]

  ## Get all indices of connected developers for the largest cluster
  idx <- which(g.clust$membership==id.largest)

  return(idx)
}


largest.subgraph <- function(graph) {
  ##############################################################################
  ## Returns a graph composed only of the largest connected component
  ## - Input -
  ## graph: igraph object
  ## - Ouput -
  ## graph.connected: igraph object, composed of the largest connected component
  ##                  provided by the input graph
  ##############################################################################
  ## Get all vertices that exist in the largest connected component
  idx <- largest.subgraph.idx(graph)
  ## Get vertices to remove
  idx.rmv <- setdiff(V(graph),idx)

  ## Remove all vertices that are not in the largest connected component
  graph.connected <- delete.vertices(graph, idx.rmv)

  return(graph.connected)
}
##========================================================================
##    						Edge Functons
##========================================================================


## These two function return the absolute number of tags received and given.
## If two persons exchange tags multiple times, this is counted repeatedly
tags.received <- function(.id, .tags) {
  return(sum(.tags[,.id]))
}

tags.given <- function(.id, .tags) {
  return(sum(.tags[.id,]))
}

## These functions, in turn, compute from/to how many _different_ developers
## a tag was given to/received from
tags.received.norep <- function(.id, .tags) {
  return(length(which(.tags[,.id]>0)))
}

tags.given.norep <- function(.id, .tags) {
  return(length(which(.tags[.id,]>0)))
}



##========================================================================
##		Community Detection
##========================================================================
oslom.community <- function(g) {
  ## uses the OSLOM progam to generate an igraph-like communities object
  ## Args:
  ##  g: igraph graph object
  ## Returns:
  ##  community: igraph-like communities object
  ##TODO: we need to figure out how to set this up better w.r.t. where the
  ##      OSLOM program files should be stored
  prog.loc <- getwd()
  file.name <- paste(prog.loc, "/oslom.dat", sep="")
  ## write graph to file
  g.frame <- get.data.frame(g, what="edges")
  g.frame["weight"] <- E(g)$weight
  write.table(g.frame, file.name, sep="\t", row.names=FALSE, col.names=FALSE)

  ## make system call to oslom
  oslom.prog <- paste(prog.loc, "/oslom_undir -w -t 1.0 -cp 1.0 -copra 10 -infomap 10 -f", sep="")
  cmd <- paste(oslom.prog, file.name, sep=" ")
  system(cmd, ignore.stdout=TRUE)

  ## read output file
  community <- read.oslom(paste(file.name, "_oslo_files/tp_without_singletons", sep=""))
  return (community)
}


link.community <- function(g){
  #########################################
  ## Description:
  ##   Utilize linkcomm package to perform network decomposition
  ##   provided with an igraph graph object
  #########################################
  ## construct a directed and weighted edge list
  edge.list.dir.we <- data.frame( cbind(get.edgelist(g), E(g)$weight))
  ## perform decomposition
  link.communities <- getLinkCommunities(edge.list.dir.we, hcmethod="single",
										 directed=TRUE, plot=FALSE, verbose=FALSE)
  ## create an igraph community object to return results in accordance with
  ## existing infrastructure in handeling communities
  node.clust <- link.communities$nodeclusters
  num.clust <- link.communities$numbers[3]
  membership <- list()
  csize      <- c()
  for (i in 1:num.clust) {
    membership[[i]] <- as.numeric(as.character(node.clust[node.clust$cluster == i,1]))
	csize[i] <- length(membership[[i]])
  }
  membership$csize <- csize
  class(membership) <- "overlapComm"
  return(membership)
}


## Clique percolation. Stolen from http://igraph.wikidot.com/community-detection-in-
## Does not work on our graph. Maybe it works in an undirected version
clique.community <- function(graph, k) {
  clq <- cliques(graph, min=k, max=k)
  edges <- c()
  for (i in seq_along(clq)) {
    for (j in seq_along(clq)) {
      if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
        edges <- c(edges, c(i,j)-1)
      }
    }
  }
  clq.graph <- simplify(graph(edges))
  V(clq.graph)$name <- seq_len(vcount(clq.graph))
  comps <- decompose.graph(clq.graph)

  lapply(comps, function(x) {
    unique(unlist(clq[ V(x)$name ]))
  })
}

## Select communities with more than .min members
select.communities.more <- function(.comm, .min) {

	if (class(.comm) == "communities"){
		N <- length(unique(.comm$membership))
		num.members <- sapply(1:(N),
				function(x) { return(length(which(.comm$membership==x))) })
		elems <- which(num.members > .min)
	}
	else if (class(.comm) == "overlapComm"){
		N <- length(.comm$csize)
		elems <- which(sapply(1:(N),
						function(x) { return (length(.comm[[x]]) > .min)  }) == TRUE)
	}

  return(elems)
}


## For large projects like the Linux kernel, the default setting of
## 25 communities necessarily leads to some very large contributions.
## Compute a more reasonable upper bound that gives the algorithm a
## change to detect reasonably small communities
compute.num.spins <- function(g) {
  AVG.SIZE <- 5

  num.spins <- as.integer(vcount(g)/AVG.SIZE)

  if (num.spins < 25) {
    num.spins <- 25
  }

  if (num.spins > 500) {
    num.spins <- 500
  }

  return(num.spins)
}


######
## Community detection algorithms require the input graph to be connected,
## here the graph is decomposed into connected components and the cluster
## algorithm is applied to each connected component, the results are then
## collected into a single membership vector
## Args:
##  g: igraph graph object
##  cluster.algo: clustering algorithm with an interface equivalent to igraph
##                clustering algorithms
## Returns:
##  community: igraph communities object, includes membership vector, size of
##             communities and number of communities
community.detection.disconnected <- function(g, cluster.algo) {
  ## Global community
  membership.global <- c()
  csize.global      <- c()

  ## Find all connected subgraphs
  g.conn.clust    <- clusters(g)
  g.conn.mem      <- g.conn.clust$membership
  g.conn.clust.no <- g.conn.clust$no ## number of clusters

  ## Perform clustering on each connected subgraph and aggregate results
  ## into single data structure
  global.idx.start <- global.idx.end <- 0
  algorithm.name <- "None"
  for (sub.g.id in 1:g.conn.clust.no) {
    ## Get all vertices for subgraph
    sub.g.verts <- which(g.conn.mem == sub.g.id)

    ## Computed connected graph
    g.conn <- induced.subgraph(g, sub.g.verts)

    if (vcount(g.conn) != 1) {
      ## Perform clustering
      g.comm <- cluster.algo(g.conn)
      comm.membership <- g.comm$membership
      algorithm.name  <- g.comm$algorithm
      csize           <- g.comm$csize
    }
    else {
      ## Singleton clusters don't require clustering
      comm.membership <- c(1)
      csize <- c(1)
    }
    comm.no <- length(unique(comm.membership))

    ## Map local cluster index system to global index system
    global.idx.start <- global.idx.end + 1
    global.idx.end   <- global.idx.end + comm.no
    comm.global.map  <- global.idx.start:global.idx.end

    ## Map local membership ids to global ids
    membership.global[sub.g.verts] <- comm.global.map[comm.membership]

    ## compute community sizes
    for (comm.id in 1:comm.no) {
      csize.global[comm.global.map[comm.id]] <- length(which(comm.membership ==
                                                             comm.id))
    }
  }
  ## Calculate other communities class attributes 
  community            <- list(membership=c(), csize=c(), modularity=0, no=0)
  community$membership <- membership.global
  community$no         <- global.idx.end
  community$csize      <- csize.global
  community$algorithm  <- algorithm.name
  class(community)     <- "communities"

  return(community)
 }


spinglass.community.connected <- function(graph, spins=compute.num.spins(graph)) {
	## wrapper for spinglass clustering algorithm

	## Description:
	## 	Spinglass will detect clusters that are disjoint and place them in one
	## 	cluster. Our definition of a community requires a connected graph. This
	## 	wrapper splits a disjoint cluster into multiple communities
	## Args:
	## 	graph: igraph graph object
	## Returns:
	## 	comms.new: an igraph communities object with connected communities
	comms.new <- list(membership=c(), csize=c(), modularity=0, no=0,
			algorithm="spinglass")
	class(comms.new) <- "communities"

	## perform normal spinglass clustering
	comms <- spinglass.community(graph, spins=spins, update.rule="config")
	## construct new communities instance
	comms.new$modularity <- comms$modularity
	## check if any communities are not disjoint
	numComms <- length(comms$csize)
	for (comm.indx in 1:numComms){
		## get vertex set corresponding to the comm.indx-th cluster
		vert.set  <- which(comms$membership == comm.indx)
		g.induced <- induced.subgraph(graph, vert.set)
		clust     <- clusters(g.induced, mode="weak")

		## if more than one cluster is found we need to split into two communities
		for (i in 1:clust$no) {
			comms.new$no <- comms.new$no + 1
			vert.set.sub <- which(clust$membership == i)
			comms.new$membership[vert.set[vert.set.sub]] <- comms.new$no
			comms.new$csize[comms.new$no] <- length(vert.set.sub)
		}
	}
	return(comms.new)
}


minCommGraph <- function(graph, comm, min=10){
	## create a graph and the associated communities object only consisting of
	## vertices that belong to communities large than a minimum size

	## Args:
	## 	graph: igraph graph object
	## 	comm:  igraph communities object
	## 	min: minimum number of vertices that must exist in a community to be
	##      kept
	##
	## Returns:
	## 	res$graph: resulting igraph graph object onces insignificant vertices
	##						are removed
	## 	res$community: resulting igraph communities object after small communities
	##								have been removed
	comm.idx <- which(comm$csize > min)
  if(length(comm.idx) == 0) {
    res <- list(graph=NULL, community=NULL)
  }
  else{
    verts <- rle(unlist(as.vector(sapply(comm.idx,
           function(x) { return(which(comm$membership==x)) }))))$values

	  V(graph)$key <- 1:vcount(graph)
	  graph.comm <- induced.subgraph(graph, verts)
	  ## use the unique key to determine the mapping of community membership to the
	  ## new graph index
	  comm$membership <- comm$membership[V(graph.comm)$key]
	  comm$csize <- sapply(1:length(comm.idx),
			  function(x) {return(length(which(comm$membership == comm.idx[x])))})
	  comm$membership <- remap.consec.seq(comm$membership)
	  res <- list(graph=graph.comm, community=comm)
  }
  return(res)
}


remap.consec.seq <- function(values) {
  ## Map an arbitrary number sequence to increase from 1
  ## For instance, map the sequence 2,2,2,3,3,4,3 to 1,1,1,2,2,3,2

  return(mapvalues(values, unique(values), 1:length(unique(values))))
}


## Select communities with less or equal than .max members
select.communities.less.equal <- function(.comm, .max) {
  N <- length(unique(.comm$membership))
  num.members <- sapply(1:(N),
			function(x) { return(length(which(.comm$membership==x))) })

  elems <- which(num.members <= .max)

  return(elems)
}

## Select communities with size between a certain range inclusive
select.communitiy.size.range <- function(.comm, bound1, bound2) {

  ## Determine which bound is the upper which is the lower
  if (bound1 <= bound2){
    lowBound = bound1
    upBound   = bound2
  }
  else{
    lowBound = bound2
    upBound  = bound1
  }

  ## Locate elements that suit the appropriate range
  N <- length(unique(.comm$membership))
  num.members <- sapply(1:(N),
			function(x) { return(length(which(.comm$membership==x))) })

  elems <- which(num.members <= upBound & num.members >= lowBound)

  return(elems)

}

## Helper function for select.communities below: Determine the largest possible
## community size that can be removed so that min.size contributors remain
select.threshold <- function(cmts, min.size) {
  cut.size <- 0

  for (i in sort(unique(cmts$size))) {
    tmp <- cmts[cmts$size > i,]$size.csum
    if (length(tmp) > 0) {
      remaining <- tmp[length(tmp)]

      if (remaining < min.size) {
        return(cut.size)
      } else {
        cut.size <- i
      }
    }
    else {
      ## The current cut would leave no contributors
      return(cut.size)
    }
  }

  return(cut.size)
}

## Remove small communities, but make sure that the fraction of
## surviving contributors on the reduced set is as least min.fract.
## Optionally, an upper bound on the community size that is deemed appropriate
## to be removed can be specified.
select.communities <- function(comm, min.fract=0.95, upper.bound=NA) {
  ## This function needs to work with community objects generated
  ## by walktrap and spinglass (the former does not produce a vsize
  ## member, so we cannot use it)
  min.size <- round(min.fract*length(comm$membership))

  comm.idx <- sort(unique(comm$membership))
  ## Provide a mapping between community labels and their size
  cmts <- data.frame(id=comm.idx,
                     size=sapply(comm.idx, function(i) {
                       sum(comm$membership==i)}))
  ## ... and sort the communities by size (they are still identifiable
  ## by their id)
  cmts <- cmts[sort(cmts$size, index.return=TRUE, decreasing=TRUE)$ix,]

  cmts$size.csum <- cumsum(cmts$size)

  cut.size <- select.threshold(cmts, min.size)

  if(!is.na(upper.bound)) {
    if (cut.size > upper.bound)
      cut.size <- upper.bound
  }

  return(cmts[cmts$size > cut.size,]$id)
}


## Summarise in which subsystems the authors of a given community are active
comm.subsys <- function(.comm, .id.subsys, N) {
  idx <- which(.comm$membership==N)
  suppressMessages(molten <- melt(.id.subsys[idx,2:dim(.id.subsys)[2]]))
  colnames(molten) <- c("Subsystem", "Fraction")

  return(data.frame(Community=N, molten))
}



## This function allows to do a sanity check for the plot
##N <- 3
##summary(id.subsys.connected[which(g.spin.community$membership==N), 2:dim(id.subsys.connected)[2]])
plot.comm.subsys <- function(.comm, .id.subsys, filename, .alg,
                             elems=1:(length(unique(.comm$membership))),
                             .height=8, .width=14) {
  comb <- vector("list", length(elems))
  for (i in 1:length(elems)) {
    comb[[i]] <- comm.subsys(.comm, .id.subsys, elems[i])
  }

  comb <- do.call(rbind, comb)

  ggsave(filename,
         ggplot(comb, aes(Subsystem, Fraction)) + geom_boxplot(outlier.colour="blue",
                                                               outlier.size=1.5, alpha=0.5) +
         facet_wrap(~Community) + geom_jitter(size=1, alpha=0.5) +
         theme(axis.text.x=element_text(angle=-90, hjust=0)) +
         labs(title=paste("Subsystem distribution for community clusters (algorithm: ", .alg,
                ")", sep="")),
         height=.height, width=.width)
}



## TODO: Dunno if that tells us something. Most likely not.
##densityplot(~pranks|group, data=construct.pr.info(g.spin.community, pr.for.all),
##            plot.points="rug")

## NOTE: When raw IDs are used as labels, they are zero-based, not
## 1-based as in the ID mapping table
## NOTE: This works for both, communities derived from walktrap and
## spinglass.
plot.group <- function(N, .tags, .iddb, .comm) {
  s <- which(.comm$membership==N)
  g <- graph_from_adjacency_matrix(.tags[s, s], mode = "directed")
  V(g)$name <- IDs.to.names(.iddb, V(g)$name)
  plot(g, vertex.label=IDs.to.names(.iddb, s))
}


## Given a single cluster of persons, construct an igraph object,
## compute some attributes for proper visualisation, and export the
## result as a graphviz dot format if a filename is provided.
save.group <- function(conf, .tags, idx, .prank, .filename = NULL, label, iddb) {
  subset <- .tags[idx, idx]

  # 🩹 Forza la struttura corretta di iddb
  iddb <- as.data.frame(iddb)
  iddb$ID <- as.numeric(iddb$ID)
  iddb$total <- as.numeric(iddb$total)
  iddb$numcommits <- as.numeric(iddb$numcommits)

  if (length(class(subset)) == 1 && class(subset) == "matrix") {
    rownames(subset) <- 1:ncol(subset)
    colnames(subset) <- 1:ncol(subset)
  }

  g <- graph_from_adjacency_matrix(subset, mode = "directed", weighted = TRUE)

  V(g)$label <- as.character(IDs.to.names(iddb, idx))

  ## 🛠️ Robust handling of .prank
  prank_vector <- if (is.list(.prank) && !is.null(.prank$vector)) {
    .prank$vector
  } else {
    .prank  # assume it's already a numeric vector
  }

  V(g)$fontsize <- scale.data(prank_vector, 15, 50)[idx]

  fc <- as.character(as.integer(100 - scale.data(log(iddb$total + 1), 0, 50)[idx]))
  V(g)$fillcolor <- paste("grey", fc, sep = "")
  V(g)$style <- "filled"
  V(g)$penwidth <- as.character(scale.data(log(iddb$numcommits + 1), 1, 5)[idx])

  if (length(label) == 1 && !is.na(label)) {
    g$label <- label
    g$fontsize <- 30
  }

  if (!is.null(.filename)) {
    g.scaled <- g
    if (ecount(g.scaled) > 0 && length(E(g.scaled)$weight) > 0) {
      scaled_weights <- log(E(g.scaled)$weight + 1)
      if (all(is.finite(scaled_weights))) {
        E(g.scaled)$weights <- scale.data(scaled_weights, 0, 100)
      } else {
        logwarn("Non-finite weights encountered during scaling")
        E(g.scaled)$weights <- rep(0, ecount(g.scaled))
      }
    } else {
      loginfo("Graph has no edges or weights; skipping weight scaling.")
    }
    write.graph(g.scaled, .filename, format = "dot")
  }

  return(g)
}


## Prepare graph data for database and insert
store.graph.db <- function(conf, baselabel, idx, .iddb, g.reg, g.tr, clusterNumber, releaseRangeId = conf$range.id) {
  message("📥 Storing cluster to DB: ", baselabel, " [cluster=", clusterNumber, "]")

  ## Create and insert cluster entry
  cluster.entry <- data.frame(
    projectId = conf$pid,
    releaseRangeId = releaseRangeId,
    clusterNumber = clusterNumber,
    label = baselabel
  )
  message("📝 cluster.entry:")
  print(cluster.entry)
  dbWriteTable(conf$con, "cluster", cluster.entry, append = TRUE, row.names = FALSE)
  message("✅ cluster entry inserted")

  ## Skip if graph is missing or invalid
  if (is.null(g.reg) || inherits(g.reg, "try-error")) {
    message("⚠️ Graph is NULL or errored — skipping edge list generation")
    return(invisible(NULL))
  }

  ## Attempt edge extraction
  edges <- tryCatch(
    get.data.frame(g.reg, what = "edges"),
    error = function(e) {
      message("⚠️ get.data.frame failed: ", e$message)
      return(data.frame(fromId = integer(0), toId = integer(0)))
    }
  )

  ## Only proceed if edges were found
  if (nrow(edges) > 0 && ncol(edges) >= 2) {
    colnames(edges)[1:2] <- c("fromId", "toId")

    ## Defensive check to ensure index bounds are valid
    max.idx <- max(idx, na.rm = TRUE)
    if (any(edges$fromId > length(idx)) || any(edges$toId > length(idx))) {
      message("⚠️ Edge indices out of bounds — skipping edge writing")
      return(invisible(NULL))
    }

    edges$fromId <- .iddb[idx[edges$fromId], "ID.orig"]
    edges$toId   <- .iddb[idx[edges$toId], "ID.orig"]

    edges <- gen.weighted.edgelist(edges)
    write.graph.db(conf, releaseRangeId, baselabel, edges, clusterNumber)
    message("✅ edge list written")
  } else {
    message("ℹ️ No edges to write for cluster ", clusterNumber)
  }
}

## Iterate over all clusters in a community decomposition, create
## a graphviz input file (via save.groups), and also store the
## cluster into the database.
## NOTE: The decomposition into clusters is identical for the
## regular and transposed pagerank, only the vertex _attributes_ (but not
## the vertices as such) differ. Edges are identical in all respects.
## Consequently, we only need to write into tables cluster_user_mapping
## and edgelist once.
save.groups <- function(conf, .tags, .iddb, .comm, .prank.list, .basedir,
                        .prefix, comm.quality, label) {
  baselabel <- label
  label.tr <- NA
  j <- 0

  for (i in unique(.comm$membership)) {
    filename.reg <- paste(.basedir, "/", .prefix, "reg_", "group_", four.digit(i),
                          ".dot", sep="")
    filename.tr <- paste(.basedir, "/", .prefix, "tr_", "group_", four.digit(i),
                          ".dot", sep="")

    if (class(.comm) == "communities") {
      idx <- as.vector(which(.comm$membership==i))

      ## Do not store clusters of size one
      if(length(idx) < 2) {
        next
      }
    }
    else if(class(.comm) == "overlapComm") {
      idx <- as.vector(.comm[[i]])
    }

    if (!is.na(baselabel)) {
      label <- paste(baselabel, i, "Community Quality = ", comm.quality[i],
			         sep=" ")
      label.tr <- paste(baselabel, "(tr)", i,  "Community Quality = ",
			            comm.quality[i], sep=" ")
    }
    g.reg <- save.group(conf, .tags, .iddb, idx, .prank.list$reg,
                        filename.reg, label)
    g.tr <- save.group(conf, .tags, .iddb, idx, .prank.list$tr,
                       filename.tr, label.tr)

    ## Store the cluster content into the database
    if (!is.na(baselabel)) {
      store.graph.db(conf, baselabel, idx, .iddb, g.reg, g.tr, j)
      j <- j + 1
    }
  }
}


## Determine which persons are in a community detection algorithm
## determined sub-community N
persons.in.group <- function(N, .comm, .iddb) {
  idlist <- which(.comm$membership==N)
  return(data.frame(ID=idlist, name=IDs.to.names(.iddb, idlist)))
}

##print("Persons in community a as determined by the spinglass algorithm")
##print(persons.in.group(7, g.spin.community, ids.connected))

## Determine the page rank factors of persons in community N,
## and merge them with the metrics information
pr.in.group <- function(N, .comm, .iddb, .pr) {
  ## Select all IDs that are contained in cluster N. Then, select all
  ## the page ranks, and combine them with the traditional metrics
  idx <- .iddb$ID %in% which(.comm$membership==N)
  return(cbind(.iddb[idx,], group=N, prank=.pr$vector[idx],
               num.members=length(which(.comm$membership==N))))
}

## Collect page rank and developer metrics for each cluster
construct.group.info <- function(.comm, .pr, .iddb, .elems) {
  res <- vector("list", length(.elems))

  ## Clusters can be empty, skip processing in this case
  if (length(.elems) == 0)
    return (res)

  for (i in 1:length(.elems)){
    res[[i]] <- pr.in.group(.elems[[i]], .comm, .iddb, .pr)
  }

  res <- do.call(rbind, res)

  res$group=as.factor(res$group)
  return (res)
}

save.cluster.stats.subsys <- function(.comm, .id.subsys, .elems,
                                      .outdir, .basename) {
  for (i in .elems) {
    print(xtable(txt.comm.subsys(.comm, .id.subsys, i)), type="latex",
          floating=FALSE, file=paste(.outdir, "/", .basename, four.digit(i), ".tex", sep=""))
  }
}

save.all <- function(conf, .tags, .iddb, .prank.list, .comm,
                     .filename.base = NULL, label, idx = NULL,
                     releaseRangeId = conf$range.id) {

  # 🛠️ Fix columns to avoid "invalid subscript type 'list'"
  .iddb <- as.data.frame(.iddb)
  .iddb$ID <- as.numeric(.iddb$ID)
  .iddb$total <- as.numeric(.iddb$total)
  .iddb$numcommits <- as.numeric(.iddb$numcommits)

  if (is.null(idx)) {
    idx <- seq_along(.iddb$ID)
  }

  # Wrap raw vectors into list with $vector field if needed
  safe.prank <- function(p) {
    if (is.list(p) && !is.null(p$vector)) return(p)
    list(vector = p)
  }

  # Build graphs for both PageRank variants
  g.all.reg <- save.group(conf, .tags, idx = as.numeric(.iddb$ID), .prank = safe.prank(.prank.list$reg),
                        .filename = NULL, label = NA, iddb = .iddb)
  g.all.tr  <- save.group(conf, .tags, idx = as.numeric(.iddb$ID), .prank = safe.prank(.prank.list$tr),
                        .filename = NULL, label = NA, iddb = .iddb)

  # Standard labels and colors
  V(g.all.reg)$label    <- .iddb$ID
  V(g.all.tr)$label     <- .iddb$ID
  V(g.all.reg)$pencolor <- V(g.all.reg)$fillcolor
  V(g.all.tr)$pencolor  <- V(g.all.reg)$fillcolor

  # Color by membership if available
  if (!is.null(.comm$membership)) {
    elems <- unique(.comm$membership)
    red <- as.integer(scale.data(0:(length(elems) + 1), 0, 255))
    for (i in elems) {
      membership_idx <- which(.comm$membership == i)
      col <- col.to.hex("#", red[i + 1], 0, 0)
      V(g.all.reg)[membership_idx]$fillcolor <- col
      V(g.all.tr)[membership_idx]$fillcolor <- col
    }
  }

  # Label fallback
  if (is.na(label) || is.null(label)) {
    label <- "Global Cluster"
  }
  g.all.reg$label <- label
  g.all.tr$label  <- label

  # Save to DB
  message("💾 Storing cluster -1: ", label, " [releaseRangeId=", releaseRangeId, "]")
  store.graph.db(conf, label, idx, .iddb, g.all.reg, g.all.tr,
                 clusterNumber = -1,
                 releaseRangeId = releaseRangeId)

  # Write Graphviz .dot files
  if (!is.null(.filename.base)) {
    write.graph(g.all.reg, paste0(.filename.base, "reg_all.ldot"), format = "dot")
    write.graph(g.all.tr,  paste0(.filename.base, "tr_all.ldot"),  format = "dot")

    # Only save .ldot for community if it's not a dummy cluster
    if (!is.null(.comm$membership) &&
        length(unique(.comm$membership)) > 1 &&
        is.character(label) &&
        !is.na(label) &&
        nzchar(label)) {

      message("✅ Saving community graph: ", label)
      save.graph.graphviz(conf$con, conf$pid, releaseRangeId, label,
                          paste0(.filename.base, "community.ldot"))
    } else {
      message("ℹ️ Skipping save.graph.graphviz(): dummy or trivial cluster [", label, "]")
    }
  }
}


## TODO: Investigate the results of inner.links and outer.links (maybe compute
## averages and min/max for all elements of a given community to see how
## well "closed" the communities are)
## TODO: Together with a subsystem distribution of the authors, this should be
## a good basis for a nice ggplot graph.
compute.community.links <- function(g, .comm, N) {
  ## TODO: Continue working here. Compute averages for inner.link and
  ## outer.link over all vertices of community N
  idx <- which(.comm$membership==N)
  function(i) {
    subspin <- spinglass.community(g.connected, vertex=V(g)[idx[i]])
    return(c(subspin$inner.links, subspin$outer.links))
  }
}
##========================================================================
##     						Page Rank
##========================================================================
compute.pagerank <- function(.tags, .damping=0.85, transpose=FALSE, weights=NULL) {
  if (transpose) {
    g <- graph_from_adjacency_matrix(t(.tags), mode="directed", weighted=weights)
  } else {
    g <- graph_from_adjacency_matrix(.tags, mode="directed", weighted=weights)
  }
  ranks <- page.rank(g, directed=TRUE, damping=.damping)

  return(ranks)
}




## Determine the N most important developers (as per the
## PageRank measure). This returns a list ordered by pagerank.
## (Note that the raw pagerank data are given by compute.pagerank()
## give the ranks ordered by ID)
influential.developers <- function(N, .ranks, .tags, .iddb) {
  if (is.na(N)) {
    N <- length(.ranks$vector)
  }

  idx = order(.ranks$vector, decreasing=TRUE)
  idlst = seq(1,length(.ranks$vector))[idx[1:N]]

  res = data.frame(name=as.character(IDs.to.names(.iddb, idlst)), ID=idlst,
    TagsGiven=sapply(idlst, function(.id) { tags.given(.id, .tags)}),
    TagsReceived=sapply(idlst, function(.id) { tags.received(.id, .tags)}),
    TagsGivenNoRep=sapply(idlst, function(.id) { tags.given.norep(.id, .tags)}),
    TagsReceiveNoRep=sapply(idlst, function(.id) { tags.received.norep(.id, .tags)}),
    rank = .ranks$vector[idx[1:N]])

  return(res)
}

store.pageranks <- function(conf, .iddb, devs.by.pr, range.id, technique) {
  ## First, create an entry in table pagerank to get the DB internal
  ## id for the releaseRangeId/technique tuple
  prank.id <- get.clear.pagerank.id(conf, range.id, technique)

  dat <- devs.by.pr[,c("ID", "rank")]
  colnames(dat) <- c("personId", "rankValue")
  ## Convert personId to in-DB values from the local indices
  dat$personId <- .iddb[dat$personId,]$ID.orig
  dat$pageRankId <- prank.id

  res <- dbWriteTable(conf$con, "pagerank_matrix", dat, append=TRUE, row.names=FALSE)
  if (!res) {
    stop("Internal error: Could not write pagerank matrix into database!")
  }
}

writePageRankData <- function(conf, outdir, .iddb, devs.by.pr, devs.by.pr.tr) {
  ## Top 20 page rank (focus on giving tags)
  print(xtable(devs.by.pr[1:20,]), type="latex", floating=FALSE,
        file=paste(outdir, "/top20.pr.tex", sep=""),
        sanitize.colnames.function=rotate.label.30)

  ## Top 20 page rank (focus on being tagged)
  print(xtable(devs.by.pr.tr[1:20,]), type="latex", floating=FALSE,
        file=paste(outdir, "/top20.pr.tr.tex", sep=""),
        sanitize.colnames.function=rotate.label.30)

  ## Emit the results into the database
  store.pageranks(conf, .iddb, devs.by.pr, conf$range.id, 0)
  store.pageranks(conf, .iddb, devs.by.pr.tr, conf$range.id, 1)
}

#########################################################################
##     					 Main Functions
#########################################################################

performAnalysis <- function(outdir, conf) {
  ################## Preprocess: force correct range.id #################
  message("📌 DEBUG: conf$pid = ", conf$pid)
  message("📌 DEBUG: conf$range.id = ", conf$range.id)
  args <- commandArgs(trailingOnly = TRUE)
  original_range_id <- as.numeric(args[length(args)])
  conf$range.id <- original_range_id
  message("🔥 conf$range.id set to ", conf$range.id)

  ################## Process the data #################
  logdevinfo("Reading files", logger = "cluster.persons")
  mat.file <- paste(outdir, "/adjacencyMatrix.txt", sep = "")
  adjMatrix <- read.table(mat.file, sep = "\t", header = TRUE)
  adjMatrix.ids <- unlist(strsplit(readLines(mat.file, n = 1), "\t"))

  colnames(adjMatrix) <- rownames(adjMatrix)
  adjMatrix <- t(adjMatrix)  # Transpose for R convention

  ids.db <- get.range.stats(conf$con, conf$range.id)
  assign(".iddb", ids.db, envir = .GlobalEnv)  # 🔥 Required for save.group()

  remapping <- unlist(lapply(adjMatrix.ids, function(id) which(id == ids.db$ID)))
  ids <- ids.db[remapping, ]
  if (!all(ids$ID == adjMatrix.ids)) {
    logerror("Id mismatch", logger = "cluster.persons")
  }

  id.subsys <- read.csv(paste(outdir, "/id_subsys.txt", sep = ""), sep = "\t", header = TRUE)
  id.subsys$ID <- id.subsys$ID + 1
  if (length(colnames(id.subsys)) == 2) {
    id.subsys <- NULL
  }

  ## Check if adjacency matrix is empty (excluding diagonal)
  mat.diag <- diag(adjMatrix)
  diag(adjMatrix) <- 0
  non.diag.sum <- sum(adjMatrix)
  diag(adjMatrix) <- mat.diag

  if (non.diag.sum == 0) {
    loginfo("⚠️ Adjacency matrix empty — forcing dummy cluster", logger = "cluster.persons")

    ids$ID.orig <- ids$ID
    ids$ID <- seq_len(nrow(ids))
    ids$Name <- as.character(ids$Name)

    pr.for.all <- rep(0, nrow(ids))
    pr.for.all.tr <- rep(0, nrow(ids))
    prank.list <- list(reg = pr.for.all, tr = pr.for.all.tr)
    dummy_comm <- list(membership = rep(1, nrow(ids)))
    idx.global <- ids$ID.orig

    # Save dummy cluster
    save.all(conf, adjMatrix, ids, prank.list, dummy_comm,
             paste(outdir, "/dummy_", sep = ""),
             label = "Dummy Cluster", idx = idx.global,
             releaseRangeId = original_range_id)

    # Save global cluster (-1)
    save.all(conf, adjMatrix, ids, prank.list,
             list(membership = rep(1, nrow(ids))),
             paste(outdir, "/global_", sep = ""),
             label = "Global (Forced)", idx = idx.global,
             releaseRangeId = original_range_id)

    return(0)
  }

  ################## Standard analysis #################
  performGraphAnalysis(conf, adjMatrix, ids, outdir, id.subsys)

  ## 🔥 Fallback: ensure global cluster is stored also in the standard path
  idx.global <- ids$ID
  prank.list <- list(reg = rep(1, length(ids$ID)), tr = rep(1, length(ids$ID)))
  save.all(conf, adjMatrix, ids, prank.list,
           list(membership = rep(1, length(ids$ID))),
           paste(outdir, "/global_", sep = ""),
           label = "Global (Forced)", idx = idx.global,
           releaseRangeId = original_range_id)
}

writeClassicalStatistics <- function(outdir, ids.connected) {
  rank.by.total <- get.rank.by.field(ids.connected, "total", 20)
  rank.by.numcommits  <- get.rank.by.field(ids.connected, "numcommits", 20)

  write.table(rank.by.total, file=paste(outdir, "/top20.total.txt", sep=""),
              sep="\t", quote=FALSE)
  print(xtable(rank.by.total), type="latex", floating=FALSE,
        file=paste(outdir, "/top20.total.tex", sep=""),
        sanitize.colnames.function=rotate.label)

  write.table(rank.by.numcommits, file=paste(outdir, "/top20.numcommits.txt",
                                    sep=""), sep="\t", quote=FALSE)
  print(xtable(rank.by.numcommits), type="latex", floating=FALSE,
        file=paste(outdir, "/top20.numcommits.tex", sep=""),
        sanitize.colnames.function=rotate.label)
}

## Detect, visualise and save communities with clustering algorithm FUN
## (for instance spinglass.communities or walktrap.communities)
detect.communities <- function(g, ids, adjMatrix, prank.list, outdir,
                               prefix, label, min.fract, upper.bound, FUN) {
  set.seed(42)
  if (vcount(g) == 1) {
    ## If there is only one vertex in the graph (which can happen for
    ## very short bug-fix cycles when a single contributor interacts
    ## with the repo), it does not make sense to try detecting communities.
    ## (besids, spinglass community detection would run into an infinite
    ## loop in this case).
    ## Note: We don't construct a community with one member, but choose
    ## the interpretation that there are no communities in this case.
    ## This possibility needs to be supported anyway because the OSLOM
    ## clustering method, for instance, can fail to detect significant communities
    ## even for larger graphs if there aren't any.
    g.community <- NULL
    elems.selected <- logical(0)
  } else {
    g.community <- community.detection.disconnected(g, FUN)
    ## compute community quality
    comm.quality <- community.metric(g, g.community, "conductance")
  }

  logdevinfo(str_c("Writing community graph sources for algorithm ", label), logger="cluster.persons")
  ## NOTE: The cluster decomposition is independent of the page
  ## rank calculation technique -- only the edge strengths, but not the
  ## page rank values influence the decomposition.
  clear.all.clusters(conf, conf$range.id, label)
  save.groups(conf, adjMatrix, ids,
              g.community, prank.list, outdir, prefix,
			  comm.quality, label=label)
  return(g.community)
}


performGraphAnalysis <- function(conf, adjMatrix, ids, outdir, id.subsys = NULL) {
  ##====================================
  ##     Find Connected Subgraphs
  ##====================================
  n <- ncol(adjMatrix)
  adjMatrix <- adjMatrix * abs(diag(1, n, n) - 1)

  logdevinfo("Computing adjacency matrices", logger = "cluster.persons")
  g <- graph_from_adjacency_matrix(adjMatrix, mode = "directed", weighted = TRUE)
  idx <- V(g)

  message("🔢 Number of nodes: ", length(idx))
  message("📋 ID list:")
  print(ids)

  ids$ID.orig <- ids$ID
  ids$ID <- seq(1, length(idx))
  ids$Name <- as.character(ids$Name)

  if (!is.null(id.subsys)) {
    message("📌 Subsystem info provided")
    id.subsys <- id.subsys[idx, ]
    id.subsys$ID <- seq(1, length(idx))
  }

  V(g)$label <- as.character(ids$Name)

  ##========================
  ##  Page rank analysis
  ##========================
  logdevinfo("Computing page rank", logger = "cluster.persons")
  pr.for.all <- compute.pagerank(adjMatrix, transpose = TRUE, weights = TRUE)
  pr.for.all.tr <- compute.pagerank(adjMatrix, .damping = 0.3, weights = TRUE)

  if (!is.numeric(pr.for.all) || anyNA(pr.for.all) || length(pr.for.all) != nrow(ids)) {
    message("⚠️ Invalid page rank vector detected (reg) — fallback to 1")
    pr.for.all <- rep(1, nrow(ids))
  }
  if (!is.numeric(pr.for.all.tr) || anyNa(pr.for.all.tr) || length(pr.for.all.tr) != nrow(ids)) {
    message("⚠️ Invalid page rank vector detected (tr) — fallback to 1")
    pr.for.all.tr <- rep(1, nrow(ids))
  }

  .prank.list <- list(reg = pr.for.all, tr = pr.for.all.tr)
  message("📈 PageRank (reg):")
  print(pr.for.all)
  message("📈 PageRank (tr):")
  print(pr.for.all.tr)

  ##========================
  ## Save developer rankings
  ##========================
  devs.by.pr <- influential.developers(NA, pr.for.all, adjMatrix, ids)
  devs.by.pr.tr <- influential.developers(NA, pr.for.all.tr, adjMatrix, ids)
  writePageRankData(conf, outdir, ids, devs.by.pr, devs.by.pr.tr)

  logdevinfo("Computing classical statistics", logger = "cluster.persons")
  writeClassicalStatistics(outdir, ids)

  ##-----------------------
  ## Save graphs to DB
  ##-----------------------
  logdevinfo("Writing the all-developers graph sources", logger = "cluster.persons")
  idx.global <- ids$ID.orig
  rrid <- conf$range.id
  if (!is.null(attr(conf, "forced_range_id"))) {
    message("⚠️ Overriding releaseRangeId with forced value")
    rrid <- attr(conf, "forced_range_id")
  }
  message("🔥 DEBUG: Saving Global cluster for releaseRangeId = ", rrid)

  message("🧮 idx.global (original IDs):")
  print(idx.global)
  message("🗂️  Release range ID: ", rrid)

  if (nrow(adjMatrix) == 0 || sum(adjMatrix) == 0) {
    message("⚠️ Adjacency matrix empty — forcing dummy cluster entry")

    .prank.list <- list(
      reg = rep(1, nrow(ids)),
      tr  = rep(1, nrow(ids))
    )

    dummy_comm <- list(membership = rep(1, nrow(ids)))

    ## Salva dummy cluster
    message("💾 Saving dummy cluster...")
    save.all(conf, adjMatrix, ids, .prank.list,
             dummy_comm,
             paste(outdir, "/dummy_", sep = ""),
             label = "Dummy Cluster", idx = idx.global, releaseRangeId = rrid)

    ## 🔥 Salva cluster globale forzato anche nel ramo dummy
    message("✅ Forcing global cluster entry (-1) in dummy branch")
    save.all(conf, adjMatrix, ids, .prank.list,
             list(membership = rep(1, length(ids$ID))),
             paste(outdir, "/global_", sep = ""),
             label = "Global (Forced)", idx = idx.global, releaseRangeId = rrid)

    return(invisible(NULL))
  }

  ##=======================
  ## Find Communities
  ##=======================
  MIN.CUT.FRACTION <- 0.95
  MAX.CUT.SIZE <- 10

  logdevinfo("Inferring communities with spin glasses", logger = "cluster.persons")
  g.spin.community <- detect.communities(g, ids, adjMatrix,
                                         .prank.list,
                                         outdir, "sg_", "Spin Glass Community",
                                         MIN.CUT.FRACTION, MAX.CUT.SIZE,
                                         spinglass.community.connected)

  logdevinfo("Inferring communities with random walks", logger = "cluster.persons")
  g.walktrap.community <- detect.communities(g, ids, adjMatrix,
                                             .prank.list,
                                             outdir, "wt_", "Random Walk Community",
                                             MIN.CUT.FRACTION, MAX.CUT.SIZE,
                                             walktrap.community)

  message("💾 Saving spin glass community...")
  save.all(conf, adjMatrix, ids, .prank.list, g.spin.community,
           paste(outdir, "/sg_", sep = ""),
           label = "Spin Glass Community", idx = idx.global, releaseRangeId = rrid)

  message("💾 Saving walktrap community...")
  save.all(conf, adjMatrix, ids, .prank.list, g.walktrap.community,
           paste(outdir, "/wt_", sep = ""),
           label = "Random Walk Community", idx = idx.global, releaseRangeId = rrid)

  ## Cluster globale (-1) sempre salvato
  message("✅ Forcing global cluster entry (-1)")
  save.all(conf, adjMatrix, ids, .prank.list,
           list(membership = rep(1, length(ids$ID))),
           paste(outdir, "/global_", sep = ""),
           label = "Global (Forced)", idx = idx.global, releaseRangeId = rrid)
}


get.community.graph <- function(graph, community, prank, ids, outdir) {
  community.idx <- sort(unique(community$membership))
  influential.people <- sapply(community.idx,
                               function(comm.idx) {
                                 which(prank$vector == max(prank$vector[which(community$membership == comm.idx)]))[1]
                               })

  names <- ids$Name[influential.people]

  g.contracted <- contract.vertices(graph, membership(community))
  E(g.contracted)$weight <- 1
  g.simplified  <- simplify(g.contracted)
  V(g.simplified)$label <- names

  ## We also use the page rank to specify the font size of the vertex
  V(g.simplified)$fontsize <- scale.data(prank$vector, 15, 50)[influential.people]

  ## The amount of changed lines is visualised by the nodes background colour:
  ## The darker, the more changes.
  ##fc <- as.character(as.integer(100-scale.data(log(.iddb$total+1),0,50)[idx]))
  V(g.simplified)$fillcolor <- paste("grey", 50, sep="")
  V(g.simplified)$style="filled"
  write.graph(g.simplified, outdir, format="dot")
}

runRandCompare <- function(nonTagDir, tagDir, outfile) {
  ## Read files for ids and adjacency matrix
  nonTagAdjMatrix <- read.table(file=paste(nonTagDir, "/adjacencyMatrix.txt",
                                  sep=""), sep="\t", header=FALSE)
  ids.nonTag <- read.csv(file=paste(nonTagDir, "/ids.txt", sep=""),
                         sep="\t", header=TRUE, stringsAsFactors = FALSE)
  ids.nonTag$ID <- ids.nonTag$ID + 1


  tagAdjMatrix <- read.table(file=paste(tagDir, "/tags.txt", sep=""),
                             sep="\t", header=FALSE)
  ids.Tag <- read.csv(file=paste(tagDir, "/ids.txt", sep=""),
                      sep="\t", header=TRUE, stringsAsFactors = FALSE)
  ids.Tag$ID <- ids.Tag$ID + 1


  colnames(nonTagAdjMatrix) <- rownames(nonTagAdjMatrix)
  colnames(tagAdjMatrix)    <- rownames(tagAdjMatrix)


  ## The adjacency file format uses a different convention for edge direction
  ## than GNU R, so we need to transpose the matrix
  nonTagAdjMatrix <- t(nonTagAdjMatrix)
  tagAdjMatrix    <- t(tagAdjMatrix)

  ## Randomize Graph
  ## Get dimension of the adjacency matrix
  dimension.nonTag <- ncol(nonTagAdjMatrix)
  idx <- 1:dimension.nonTag
  ## Randomize the indecies
  idx.rand <- sample(idx, replace=FALSE)
  ## Randomize the adjacency matrix
  nonTagAdjMatrixRand <- nonTagAdjMatrix[idx.rand, idx.rand]

  ## Run comparison on randomized adjacency matrix
  graphComparison(nonTagAdjMatrixRand, ids.nonTag, tagAdjMatrix,
                  ids.Tag, outfile)
}


write.graph.2.file <- function(.filename, g, .iddb, idx) {
  V(g)$label <- as.character(IDs.to.names(.iddb, idx))

  write.graph(g, .filename, format="dot")
}

#########################################################################
                                        #     					 Experiments
#########################################################################
experiment <- function(g, g.connected){
                                        # Somce other generic graph measures. TODO: See if/how we can use them
  max(closeness(g))
  max(betweenness(g))

                                        # NOTE: Computing the adjacency graphs with weighted instead
                                        # of multiple edges could be done with
                                        # g <- graph.adjaceny(tags, mode="directed", weighted=TRUE)
                                        # but takes quite a long time (and is also not suitable for most
  o# community detection algorithms)

  ranks <- page.rank(g)
                                        # Map the page rank values to [0,100]
  ranks.norm <-  ranks
  ranks.norm$vector <- scale.data(ranks.norm$vector, 0, 100)

  test <-  clique.community(g.connected, 5)

}

#########################################################################
##     					 Testing Section
#########################################################################
test.community.quality <- function() {

  r.1 <- c(0,1,1,1,1,0,0,0)
  r.2 <- c(1,0,1,1,0,0,0,0)
  r.3 <- c(1,1,0,1,0,0,0,1)
  r.4 <- c(1,1,1,0,0,0,0,0)
  r.5 <- c(1,0,0,0,0,1,0,0)
  r.6 <- c(0,0,0,0,1,0,1,0)
  r.7 <- c(0,0,0,0,0,1,0,1)
  r.8 <- c(0,0,1,0,0,0,1,0)

  adj.matrix <- matrix(data = c(r.1,r.2,r.3,r.4,r.5,r.6,r.7,r.8), ncol = 8, nrow = 8)

  g <- graph_from_adjacency_matrix(adj.matrix)


  ## Test that modularity is correct
  g.spincommunity <- spinglass.community(g)
  igraph.modularity.result <- modularity(g, g.spincommunity$membership)
  modularity.result        <- sum(community.metric(g, g.spincommunity, "modularity"))
  if( !(igraph.modularity.result == modularity.result)){
    logerror("modularity test failed", logger="cluster.persons")
  }
  else{
    logdevinfo("Success: modularity test passed", logger="cluster.persons")
  }
}

test.community.quality.modularity <- function() {

  ##        1,2,3,4,5,6,7,8
  r.1 <- c(0,0,1,0,0,0,0,0)
  r.2 <- c(0,0,1,0,0,0,0,0)
  r.3 <- c(0,0,0,0,1,0,0,0)
  r.4 <- c(1,0,0,0,1,0,0,0)
  r.5 <- c(0,0,0,0,0,0,0,0)
  r.6 <- c(0,0,0,0,0,0,1,1)
  r.7 <- c(0,0,0,0,0,0,0,0)
  r.8 <- c(0,0,0,0,1,0,1,0)

  adj.matrix <- t(matrix(data = c(r.1,r.2,r.3,r.4,r.5,r.6,r.7,r.8), ncol = 8, nrow = 8))

  g <- graph_from_adjacency_matrix(adj.matrix)
  g.clust <- list()
  g.clust$membership <- c(1,1,1,2,2,3,3,3)

  quality <- community.metric(g, g.clust, "modularization")

}

#########################################################################
##     					 Executed Statements
#########################################################################
##----------------------------
## Parse commandline arguments
##----------------------------

config.script.run({
  conf <- config.from.args(positional.args=list("resdir", "range.id"),
                           require.project=TRUE)
  loginfo(str_c("DEBUG: conf$pid = ", conf$pid))
  performAnalysis(conf$resdir, conf)
})
