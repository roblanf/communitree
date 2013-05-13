context("ancstates")

test_that("I can reconstruct colonisation events", {
	
	#setup
	#simulate nice data from an equal-rates model (thanks to Liam Revell's blog for this)
	tree <- pbtree(n=50,scale=1)
	Q <- matrix(c(-1,1,1,-1),2,2)
	rownames(Q)<-colnames(Q)<-c(0,1)
	x <- sim.history(tree,Q,anc='0')$states #make sure the ancestral state is zero		
	ctree <- communitree(tree, x, nsim=10)
	species.in.loc1 = length(which(x==1))
	ctree <- reconstruct.colonisations(ctree)
	colonisations <- ctree$colonisations
	
	#tests
	#all times should be positive, e.g. 5mya
	expect_that(length(which(colonisations$summary$start.time<0)), equals(0))	
	expect_that(length(which(colonisations$summary$end.time<0)), equals(0))	

	#the total number of species resulting from colonistaions should equal the number of species in location 1
	expect_that(species.in.loc1, equals(sum(colonisations$summary$num.spp)))
	
	
})