require(phytools)

context("ancstates")

test_that("I can successfully build a communitree object", {
	
	#setup	
	#simulate nice data from an equal-rates model (thanks to Liam Revell's blog for this)
	tree <- pbtree(n=50,scale=1)
	Q <- matrix(c(-1,1,1,-1),2,2)
	rownames(Q)<-colnames(Q)<-c(0,1)
	x <- sim.history(tree,Q)$states		
	ctree <- communitree(tree, x, nsim=5)
	
	non.ultrametric.tree <- rtree(50)
	multifurcating.tree<-di2multi(tree,tol = sort(tree$edge.length[tree$edge[,2]>length(tree$tip)])[3])
	unrooted.tree <- unroot(tree)
	
	#tests
	expect_that(ctree, is_a("communitree"))
	expect_that(ctree$tree, equals(tree))
	expect_that(as.vector(ctree$tip.locations$locations), equals(as.vector(x)))
	expect_that(names(ctree)[1], equals("tree"))
	expect_that(names(ctree)[2], equals("tip.locations"))
	expect_that(names(ctree)[3], equals("anc.node.freqs"))
	expect_that(names(ctree)[4], equals("anc.node.ml"))
	expect_that(names(ctree)[5], equals("branching.times"))
	
	expect_that(communitree(non.ultrametric.tree, x), throws_error("Tree not ultrametric, please fix"))
	expect_that(communitree(multifurcating.tree, x), throws_error("Tree not binary, please fix"))
	expect_that(communitree(unrooted.tree, x), throws_error("Tree not rooted, please fix"))
	

})