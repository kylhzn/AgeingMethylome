#Profiling drift
# input betavalues, vec_age
# output as table row = probeID, columns = spearman rho, p-value 


#1) probes that methylate with age
#a. Run spearmans rho across a data set of probes


#1.
# input a matrix of beta_values, row = probes, columns = samples, require a vector of ages for samples in same order
# > need age meta data
	install.packages("MonoPoly")
	install.packages("polynom")
	library(MonoPoly)
	library("polynom")

profile = function(vec_age, betavalues, alpha = 0.05){
	#vec_age 
	#betavalues
	#alpha
	
	#define m value function
	m.values = function(x){
	  log2(x/(1-x))
	}
	#2.
	# convert to m values
	mvalues <- m.values(betavalues)

	#which probes dirft with age x-axis = age of all samples, y-axis = m value for a samples probe
	# i will run on each column
	d <- dim(mvalues)
	m1 <- matrix(nrow = d[1],ncol = 2) 
	for(i in 1:d[1]){
	out <- cor.test(vec_age,mvalues[i,],method="spearman")
	m1[i,1] <- out$estimate
	m1[i,2] <- out$p.value
	}
	rm(out)
	#set rownames of m1 to probe IDs
	row.names(m1) <- row.names(betavalues)

	# bonferroni correction
	b <- (alpha/d[2])

	#sort by p-value and keep only significant probes
	m1 <- m1[sort.list(m1[,2]), ]
	m1 <- m1[m1[,2] < b,]

	#rn is a object of the significant probe ID, subset mvalue by rn
	rn <- rownames(m1)
	sigmvalues <- mvalues[rn,]



	#b. If spearman rho is statistically significant, fit a regression line to define direction LMLIST
	LMLIST <- list()
	for(i in 1:dim(sigmvalues)[1]){
	LMLIST(i) <- lm(sigmvalues[i,]~vec_age)
	}

	# and Monotonic Polonomial LSR package, includes monpol$RSS and monpol$scly which is difference in max and min of y MPLIST
	
	MPLIST <- list()
	for(i in 1:dim(sigmvalues)[1]){
	MPLIST[i] <- monpol(sigmvalues[i,]~vec_age)
	}

	#saving the monotone function, POLYNOMIALS
	
	mplistbeta <- list()
	for(i in 1:length(mplist){
	tmp <- mplist[i]
	mplistbeta[[i]] <- tmp$beta.raw
	}
	POLYNOMIALS <- list()
	for(i in 1:length(mplist){
	tmp <- mplistbeta[i]
	POLYNOMIALS[[i]] <- polynomial(tmp[[1]])
	}
	names(POLYNOMIALS) <- rn

	#plotting the monotome function for a probe, cpg######plot.jpg
	for(i in 1:length(POLYNOMIALS)){
	jpeg(paste0(rn[i],"plot.jpg")
	plot(sigmvalues[i,]~vec_age)
	lines(POLYNOMIALS[[i]])
	dev.off()
	}
return(LMLIST)
return(MPLIST)
return(POLYNOMIALS)
	#c. Desired probes possess positive correlation

	#2) probe that demethylate with age
	#a. Following from 2b, desired probes possess a negative correlation
	#3) probes that are static with age (low variance among all beta values)
	#a. Folling from 2b if not significant, run a variance analysis for each probe across ages
	#b. low variance represent static probes
	#4) probes that have a high variance but not correlated
	#a. Follow from 3a, if high variance these are not correlated represent highly variable probes
	#5) plot gene regions and cpg regions for 1 2 3 4 and compare

	}
