library(phytools)


get.descendant.spp <- function(node,ctree){
	#returns the species in location 1 that descend from any given node in the tree
	tree <- ctree$tree

	#get all species that descend from a node
	descendant.nodes <- getDescendants(tree, node)
	if(length(descendant.nodes)==0){descendant.nodes=node} #this happens if it's a tip
	descendant.tips <- descendant.nodes[which(descendant.nodes<=length(tree$tip.label))]
	descendant.spp <-tree$tip.label[descendant.tips]
	
	#now just keep the species in location 1
	loc1 <- rownames(ctree$tip.locations[ctree$tip.locations$locations==1,])	
	descendant.spp <- descendant.spp[descendant.spp %in% loc1]
	
	descendant.spp <- list(descendant.spp)
	names(descendant.spp) <- node
	return(descendant.spp)
}

reconstruct.colonisations <- function(ctree){
	
	#let's first get a matrix of node and tip locations
	tip.locs <- as.matrix(ctree$tip.locations$locations) 
	rownames(tip.locs) <- ctree$tip.locations$tip.number
	all.locations <- rbind(ctree$anc.node.ml, tip.locs)
		
	#this gives us a matrix of zeros and ones in place of the node numbers, see
	#http://stackoverflow.com/questions/15693340/replace-matrix-values-in-r-with-rownames-another-matrix
	states <- matrix(all.locations[match(ctree$tree$edge, rownames(all.locations) )], ncol=2)	
	
	# and this is TRUE/FALSE matrix with the same number of rows as edges, 
	colonisation.edges <- states[,1]==0 & states[,2]==1

	#these are the start/end times of each node in the tree
	times <- matrix(ctree$branching.times[match(ctree$tree$edge, rownames(ctree$branching.times) )], ncol=2)
	times[which(is.na(times))] <- 0

	# now we can build a data frame with all the information for each colonisation event
	colonisation.events <- as.data.frame(cbind(colonisation.edges,ctree$tree$edge,times))
	colnames(colonisation.events) <- c("colonisation","start.node","end.node","start.time","end.time")
	
	#now we thin out that matrix (this is all a bit inefficient...) to just the colonisation events
	colonisation.events <- colonisation.events[which(colonisation.events$colonisation==1),]	
	colonisation.events$colonisation.branch <- rownames(colonisation.events)
	colonisation.events <- colonisation.events[,-1]
	rownames(colonisation.events) <- NULL

	colonisation.spp <- apply(as.matrix(colonisation.events$end.node), 1, get.descendant.spp, ctree)	
	colonisation.spp <- unlist(colonisation.spp,recursive=FALSE)
	
	#now we add the total number of species in each colonisation event to colonistion.events
	colonisation.events$num.spp <- unlist(lapply(colonisation.spp, length))

	colonisations <- list(colonisation.events,colonisation.spp)
	names(colonisations) <- c('summary', 'events')

	return(colonisations)
}





reconstruct.speciations <- function(ctree){
	
	
	
	
	
}