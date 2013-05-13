
get.location.characters <- function(tree, dataset){

	#set up a character matrix which is '1' if in particularl location, and '0' if not 
	location.characters <-  tree$tip.label %in% dataset$tree.name
	location.characters <- replace(location.characters, which(location.characters==TRUE), 1)
	location.characters <- replace(location.characters, which(location.characters==FALSE), 0)
	names(location.characters) <- tree$tip.label 
	return(location.characters)
}


get.node.frame <- function(tree, dataset, ace.method="ER"){
	
	# take a tree and a data frame of tips, mindate, maxdate, and return a DF 
	# that has the ML designations of each node in the tree: 
	# '1' if it's reconstructed as within the focal location
	# '0' if it's reconstructed as somewhere else
	
	root.node <- length(tree$tip.label) + 1
	
	# 1. we build a simple data.frame of location characters, where a taxon has a '1' if it's 
	# in the focal location, and a zero if it's not
	location.characters <- get.location.characters(tree, dataset)
	
	# 2. reconstruct ancestral states on that tree using any method in ACE, specified with ace.method
	ACE.result <- ace(location.characters, tree, type="discrete", model=ace.method)
		
	# 3. Assign ML internal node states (location.one is the focal location)
	in.location.one <- ACE.result$lik.anc[,2]
	in.location.one <- replace(in.location.one, which(in.location.one<=0.5), 0)
	in.location.one <- replace(in.location.one, which(in.location.one>0.5), 1)
	
	# 4. make a DF of node numbers ordered as post-order tree-traversal, matched with phylo objects numbers ($int.nodes)
	nums <- 1:tree$Nnode
	node.frame <- as.data.frame(in.location.one, row.names=nums)
	node.frame$int.nodes <- root.node:(root.node+tree$Nnode-1)
	#add branching times to data frame
	node.frame$branching.times <- branching.times(tree)
	
	return(node.frame)
	
	}
	
reconstruct.spp.additions <- function(tree, dataset, node.frame, colonisation.placement='CENTRAL', plot.trees=FALSE){
	
	
	#identify as 'in situ speciation' events, those nodes (but not terminal nodes...) whose state is 1, 
	# and for whom both descendent also have a state of 1
	in.situ.times <- c()
	in.situ.nodes <- node.frame$int.nodes[node.frame$in.location.one==1]
	for(node in in.situ.nodes){	
		dec.nodes <- get.descendent.nodes(tree, node)
		dec.states <- get.states(node.frame, location.characters, tree, dec.nodes)
		if(sum(dec.states)==2){ in.situ.times <- append(in.situ.times, branching.times(tree)[which(names(branching.times(tree)) == node)]) }
	}
									   
	
	#identify 'colonisation events' as those branches which start as 0 and end up as 1
	colonisation.times <- c() #a list of dates, labelled with descendent nodes
	
	for(node in node.frame$int.nodes){
		# get the state of the ancestral node
		anc.state <- get.states(node.frame, location.characters, tree, node)
				
		# get descendent nodes as their internal node numbers, and their states too
		descendents <- tree$edge[,2][tree$edge[,1]==node]
	
		for(descendent in descendents){
			#get the descendent state
			dec.state <- get.states(node.frame, location.characters, tree, descendent)
	
			#cat("anc.state: ", anc.state)
			#cat("dec.state: ", dec.state)
			#cat("node: ", node)
	
			#if the ancestor is 0 and the descendent is 1:
			if(anc.state-dec.state == -1){
	
				if(descendent<root.node){is.tip<-TRUE}
				else{ is.tip<-FALSE }
				
				
				#get the node times of the ancestor and descendent
				anc.time <- node.frame$branching.times[node.frame$int.nodes==node]
				#if the descendent was a tip, then its time will be zero
				if(is.tip){dec.time<- 0}
				else { dec.time <- node.frame$branching.times[node.frame$int.nodes==descendent] }
				
				#tips may have additional data on arrival times in the data file, and so the date may be more constrained than just from the tree
				if(descendent<root.node){ 			# this defines a tip in the tree
	
					#get max date from data file
					max.date <- get.date("max.date", names(location.characters[descendent]), bird.data, tree)
					if(!is.na(max.date) && max.date<anc.time){anc.time<-max.date}
					
					#get min date from data file
					min.date <- get.date("min.date", names(location.characters[descendent]), bird.data, tree)
					if(!is.na(min.date) && min.date>dec.time){dec.time<-min.date}
	
				}
				
				#add the colonisation date to the list of dates, and name it with the ancestral node name, for checking
				colonisation.date <- get.colonisation.time(anc.time, dec.time, method=method)
				names(colonisation.date) <- which(tree$edge[,2]==descendent)          #name the colonisation date with the edge label
				colonisation.times <- append(colonisation.times, colonisation.date)
			}
		}
	}
	
	##############################################################################################
	############### USE GEO-TREE and DATA FILE TO ADD IN DATA FOR SPP NOT IN THE TREE ############
	
	#first, get a list of all the species which werent in the tree, so havent yet been accounted for
	additional.spp <- which(is.na(bird.data$tree.name))
	
	#the extra information comes in three possible forms: (i) NONE, (ii) a known date of addition, (iii) a sister group
	#deal with each species in turn using the apply() function
	additional.spp.times <- c()
	for(id in additional.spp){
		
		#for the max.date column, the date can be extracted with the get.date function
		#this will be either the numeric date, NA, or the MRCA date of the specified taxa
		max.date <- parse.date(bird.data[id, 'max.date'])
		if(is.na(max.date)) max.date <- -Inf 
		min.date <- parse.date(bird.data[id, 'min.date'])
		if(is.na(min.date)) min.date <- 0
	
		sister.spp <- as.character(bird.data[id, 'sister.species'])
	
		# for the sister spp column two things may happen:
		# (i) if its >1 spp, then we need to take the next node down from that (unless were already at the root)
		if(length(strsplit(sister.spp, ",")[[1]])>1){
			current.node <- which(branching.times(tree)==sister.spp.date)
			if(current.node != root.node){
				current.node <- tree$edge[,1][tree$edge[,2]==current.node]
				sister.spp.date <- node.frame$branching.times[node.frame$int.nodes==current.node] 
			}
			else sister.spp.date <- node.frame$branching.times[node.frame$int.nodes==root.node] #if the current node was root
		}
		# (ii) if its just a single spp, we can use the date from get.date:
		else sister.spp.date <- parse.date(sister.spp)
		
		if(is.na(sister.spp.date)){ sister.spp.date <- -Inf }
		
		#there can be only one max date, so fix this
		max.date <- max(max.date, sister.spp.date)
		   
		#now get the arrival time for that species using the min and max dates
		new.time <- get.colonisation.time(max.date, min.date, method=method)
		
		additional.spp.times <- append(additional.spp.times, new.time)
	}	
	
	all.arrival.times <- list(in.situ.times, colonisation.times, additional.spp.times)

	# plot the trees as you get them. This is useful if for checking that you have sensible answers, 
	# but avoid it in general (makes the script much slower than it has to be 
	if(plot.trees) plot.colonisation.tree(tree, in.situ.times, colonisation.times, additional.spp.times, location.characters, do.density=TRUE)
	
	return(all.arrival.times)

}
