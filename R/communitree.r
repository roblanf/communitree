#' Build a communitree object from an ape phylo object, and a binary location vector

require(phytools)

getStates<-function(x){
	# this function is taken from Liam Revell's blog here: 
	# http://blog.phytools.org/2013/03/a-little-more-on-ancestral-state.html
	y<-setNames(sapply(x$maps,function(x) names(x)[1]),x$edge[,1])
	y<-y[as.character(length(x$tip)+1:x$Nnode)]
	return(y)
}


check.tree<- function(tree){
	# this function performs basics checks on a phylogenetic tree
	# to make sure it'll work with the rest of the package
	if(!is.rooted(tree)){ 
		stop("Tree not rooted, please fix") 
	}
	if(!is.binary.tree(tree)){ 
		stop("Tree not binary, please fix")
	}
	if(!is.ultrametric(tree)){
		stop("Tree not ultrametric, please fix")
	}
}

communitree <- function(tree, locations,nsim=500,model="ER"){

	check.tree(tree)
	
	# 1. we estimate ancestral states the right way
	# http://blog.phytools.org/2013/03/conditional-scaled-likelihoods-in-ace.html
	mtrees<-make.simmap(tree,locations,nsim=nsim,model=model)
	AA<-sapply(mtrees,getStates)
	anc.node.freqs <- t(apply(AA,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=c(0,1), Nsim=length(mtrees)))
	
	# 2. calculate the ML states - resolve conflicts assuming that a node is in the focal location
	anc.node.ml <- matrix(nrow=nrow(anc.node.freqs), 1)
	anc.node.ml[which(anc.node.freqs[,1]>anc.node.freqs[,2])] = 0
	rownames(anc.node.ml) <- rownames(anc.node.freqs)	
	
	# 3. add in the observed states at tips, labelled according to the numbers in the tree$edge matrix
	tip.number <- match(tree$tip.label, names(locations))
	tip.locations <- as.data.frame(cbind(locations,tip.number))
	
	# 4. build communitree object
	ctree <- list()
	ctree$tree <- tree
	ctree$tip.locations <- tip.locations
	ctree$anc.node.freqs <- anc.node.freqs
	ctree$anc.node.ml    <- anc.node.ml
	ctree$branching.times <- as.matrix(branching.times(tree))
	
	class(ctree) <- "communitree"
	
	return(ctree)
}

plot.communitree <- function(ctree, method='freqs'){
	if(method=='freqs'){
		# this code is also from Liam Revell's blog:
		# http://blog.phytools.org/2013/03/a-little-more-on-ancestral-state.html
		plot(ctree$tree,no.margin=TRUE,show.tip.label=FALSE, edge.width=2)
		nodelabels(pie=ctree$anc.node.freqs,piecol=setNames(c("blue","red"),colnames(ctree$anc.node.freqs)),cex=0.6)
		tips<-sapply(c('0','1'),"==",ctree$tip.locations)*1
		tiplabels(pie=tips,piecol=setNames(c("blue","red"),colnames(ctree$anc.node.freqs)),cex=0.6)
	}	

	if(method=='ml'){
		plot(ctree$tree,no.margin=TRUE,show.tip.label=FALSE, edge.width=2)
		nodes <- sapply(c('0','1'),"==",ctree$anc.node.ml)*1
		nodelabels(pie=nodes,piecol=setNames(c("blue","red"),colnames(ctree$anc.node.freqs)),cex=0.6)
		tips<-sapply(c('0','1'),"==",ctree$tip.locations)*1
		tiplabels(pie=tips,piecol=setNames(c("blue","red"),colnames(ctree$anc.node.freqs)),cex=0.6)
	}


}