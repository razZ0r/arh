# --------------------------------------------------------------------
ARH = function(x = "numeric", y = "numeric", f = "character", na.rm = FALSE) {
# INPUT
#	x, y: exon expression; numeric vectors, x,y >= 0
#	f : gene IDs corresponding to x,y; groups the exon expressions to genes, character
#	x, y, f all of the same length
#	na.rm: how do you want to treat NA in the expressions? R standard is here implemented, but we use 'na.rm = TRUE'.
# OUTPUT
#	ARH values, numeric vector, length(ARH) == length(unique(f)), names(ARH) == names(unique(f))

	# negative exon expressions not interpretable
	if(any(x < 0, na.rm = TRUE) | any(y < 0, na.rm = TRUE)) stop("What are negative expressions?")
	
	# avoid division by zero and arbitrarily big numbers
	y[ y < 0.0001 ] = 0.0001

	# compute splicing probabilities, each exon compared to gene
	splicingDeviations = log2(x / y)
	splicingDeviationsMedian = split(splicingDeviations, f)
	splicingDeviationsMedian = sapply(X = splicingDeviationsMedian, FUN = median, na.rm = na.rm)
	splicingDeviationsMedian = splicingDeviationsMedian[ match(f, names(splicingDeviationsMedian)) ]
	splicingDeviations = 2^abs( splicingDeviations - splicingDeviationsMedian )
	rm(splicingDeviationsMedian)
	splicingProbabilitiesSum = split(splicingDeviations, f)
	splicingProbabilitiesSum = sapply(splicingProbabilitiesSum, sum, na.rm = na.rm)
	splicingProbabilitiesSum = splicingProbabilitiesSum[ match(f, names(splicingProbabilitiesSum)) ]
	splicingProbabilities = splicingDeviations / splicingProbabilitiesSum
	rm(splicingDeviations, splicingProbabilitiesSum)
	
	# compute ARH, each gene	
	entropy = split(splicingProbabilities, f)
	entropy = entropy[ match(unique(f), names(entropy)) ]
	entropy = sapply(X = entropy, FUN = function(X) return( -sum(X * log2(X), na.rm = na.rm) ))
	
	iqrQuotient = x / y
	iqrQuotient = split(iqrQuotient, f)
	iqrQuotient = iqrQuotient[ match(unique(f), names(iqrQuotient)) ]
	iqrQuotient = sapply(X = iqrQuotient, FUN = quantile, probs = c(0.25, 0.75), na.rm = TRUE)
	iqrQuotient = iqrQuotient[ "75%" ,  ] / iqrQuotient[ "25%" ,  ]

	geneLength = table(f)
	geneLength = geneLength[ match(unique(f), names(geneLength)) ]

	arh = as.numeric(iqrQuotient * (log2(geneLength) - entropy))
	names(arh) = unique(f)
	
	# filter for at least two finite expression values
	# (1) no splicing in genes with only one exon
	# (2) infinite values lead to extreme splicing indication, hard to interprete
	good = split(is.finite(x) & is.finite(y) & is.finite(splicingProbabilities), f = f)
	good = sapply(good, function(X) return(sum(X) >= 2))
	good = good[ match(unique(f), names(good)) ]
	arh[ !good ] = NA
	
	return(arh)
}
# EXAMPLE
#x = c(1, 1, 1, 2, 2, Inf, 3, 4, 4, 5)
#y = c(5, 6, 7, 2, 0, 5, 5, NA, 5, 2)
#f = c("b", "b", "b", "a", "a", "c", "c", "d", "d", "e")
#ARH(x = x, y = y, f = f)
#ARH(x = x, y = y, f = f, na.rm = TRUE)
# EXAMPLE, compare to OpenOfficeCalc
#x = c(
#	1.5, 1.2, 0.85, 2, 1.9, 1.8, 0.9,
#	1.5, 1.2, 0.2, 2, 1.9, 1.8, 0.9, 
#	0.7, 0.9, 0.75, 0.8, 1.9, 1.8, 2.1
#)
#y = rep(1, 3*7)
#f = c(rep("Sheet1", 7), rep("Sheet2", 7), rep("Sheet3", 7))
#ARH(x = x, y = y, f = f)

# --------------------------------------------------------------------
ARH_sd = function(x = "numeric", y = "numeric", f = "character", na.rm = FALSE) {
# INPUT
#	x, y: exon expression; numeric vectors, x,y >= 0
#	f : gene IDs corresponding to x,y; groups the exon expressions to genes, character
#	x, y, f all of the same length
#	na.rm: how do you want to treat NA in the expressions? R standard is here implemented.
# OUTPUT
#	numeric vector, length(ARH_sd) == length(f), names(ARH_sd) == names(f)

	# negative exon expressions not interpretable
	if(any(x < 0, na.rm = TRUE) | any(y < 0, na.rm = TRUE)) stop("What are negative expressions?")
	
	# avoid division by zero and arbitrarily big numbers
	y[ y < 0.0001 ] = 0.0001

	# compute splicing probabilities, each exon compared to gene
	splicingDeviations = log2(x / y)
	splicingDeviationsMedian = split(splicingDeviations, f)
	splicingDeviationsMedian = sapply(X = splicingDeviationsMedian, FUN = median, na.rm = na.rm)
	splicingDeviationsMedian = splicingDeviationsMedian[ match(f, names(splicingDeviationsMedian)) ]
	splicingDeviations = splicingDeviations - splicingDeviationsMedian
	rm(splicingDeviationsMedian)
	
	# filter for at least two finite expression values
	# (1) no splicing in genes with only one exon
	# (2) infinite values lead to extreme splicing indication, hard to interprete
	good = split(is.finite(x) & is.finite(y) & is.finite(splicingDeviations), f = f)
	good = sapply(good, function(X) return(sum(X) >= 2))
	good = good[ match(f, names(good)) ]
	splicingDeviations[ !good ] = NA
	
	return(splicingDeviations)
}
# EXAMPLE
#x = c(1, 1, 1, 2, 2, Inf, 3, 4, 4, 5)
#y = c(5, 6, 7, 2, 0, 5, 5, NA, 5, 2)
#f = c("b", "b", "b", "a", "a", "c", "c", "d", "d", "e")
#ARH_sd(x = x, y = y, f = f)
#ARH_sd(x = x, y = y, f = f, na.rm = TRUE)
# EXAMPLE, compare to OpenOfficeCalc
#x = c(
#	1.5, 1.2, 0.85, 2, NA, NA, 0.9,
#	1.5, 1.2, 0.2, 2, 1.9, 1.8, 0.9, 
#	0.7, 0.9, 0.75, 0.8, 1.9, 1.8, 2.1
#)
#y = rep(1, 3*7)
#f = c(rep("Sheet1", 7), rep("Sheet2", 7), rep("Sheet3", 7))
#ARH_sd(x = x, y = y, f = f)
#ARH_sd(x = x, y = y, f = f, na.rm = TRUE)
# EXAMPLE, compare to OpenOfficeCalc
#x = c(
#	1.5, 1.2, 0.85, 2, 1.9, 1.8, 0.9,
#	1.5, 1.2, 0.2, 2, 1.9, 1.8, 0.9, 
#	0.7, 0.9, 0.75, 0.8, 1.9, 1.8, 2.1
#)
#y = rep(1, 3*7)
#f = c(rep("Sheet1", 7), rep("Sheet2", 7), rep("Sheet3", 7))
#ARH_sd(x = x, y = y, f = f)

# --------------------------------------------------------------------
ARH_p = function(x = "numeric", y = "numeric", f = "character", na.rm = FALSE) {
# requires package "evd"
# INPUT
#	x, y: exon expression; numeric vectors, x,y >= 0
#	f : gene IDs corresponding to x,y; groups the exon expressions to genes, character
#	x, y, f all of the same length
#	na.rm: how do you want to treat NA in the expressions? R standard is here implemented, but we use 'na.rm = TRUE'.
# OUTPUT
#	numeric vector, 0 <= ARH_p <= 1, p-values of the generalised extreme value distribution fit for the biologic ARH background, 
#	length(ARH) == length(unique(f)), names(ARH) == names(unique(f))

	# negative exon expressions not interpretable
	if(any(x < 0, na.rm = TRUE) | any(y < 0, na.rm = TRUE)) stop("What are negative expressions?")
	
	# avoid division by zero and arbitrarily big numbers
	y[ y < 0.0001 ] = 0.0001

	# compute splicing probabilities, each exon compared to gene
	splicingDeviations = log2(x / y)
	splicingDeviationsMedian = split(splicingDeviations, f)
	splicingDeviationsMedian = sapply(X = splicingDeviationsMedian, FUN = median, na.rm = na.rm)
	splicingDeviationsMedian = splicingDeviationsMedian[ match(f, names(splicingDeviationsMedian)) ]
	splicingDeviations = 2^abs( splicingDeviations - splicingDeviationsMedian )
	rm(splicingDeviationsMedian)
	splicingProbabilitiesSum = split(splicingDeviations, f)
	splicingProbabilitiesSum = sapply(splicingProbabilitiesSum, sum, na.rm = na.rm)
	splicingProbabilitiesSum = splicingProbabilitiesSum[ match(f, names(splicingProbabilitiesSum)) ]
	splicingProbabilities = splicingDeviations / splicingProbabilitiesSum
	rm(splicingDeviations, splicingProbabilitiesSum)
	
	# compute ARH, each gene	
	entropy = split(splicingProbabilities, f)
	entropy = entropy[ match(unique(f), names(entropy)) ]
	entropy = sapply(X = entropy, FUN = function(X) return( -sum(X * log2(X), na.rm = na.rm) ))
	
	iqrQuotient = x / y
	iqrQuotient = split(iqrQuotient, f)
	iqrQuotient = iqrQuotient[ match(unique(f), names(iqrQuotient)) ]
	iqrQuotient = sapply(X = iqrQuotient, FUN = quantile, probs = c(0.25, 0.75), na.rm = TRUE)
	iqrQuotient = iqrQuotient[ "75%" ,  ] / iqrQuotient[ "25%" ,  ]

	geneLength = table(f)
	geneLength = geneLength[ match(unique(f), names(geneLength)) ]

	arh = as.numeric(iqrQuotient * (log2(geneLength) - entropy))
	names(arh) = unique(f)
	
	# filter for at least two finite expression values
	# (1) no splicing in genes with only one exon
	# (2) infinite values lead to extreme splicing indication, hard to interprete
	good = split(is.finite(x) & is.finite(y) & is.finite(splicingProbabilities), f = f)
	good = sapply(good, function(X) return(sum(X) >= 2))
	good = good[ match(unique(f), names(good)) ]
	arh[ !good ] = NA
	
	require("evd")
	arh[ !is.na(arh) ] = pgev(arh[ !is.na(arh) ], loc = 0.00633807, scale = 0.00550684, shape = 0.332936, lower.tail = FALSE)

	return(arh)
}
# EXAMPLE
#	p1 = pgev(c(0.02285041, 0.03072437, 0.05717684, 0.1308201), loc = 0.004776, scale = 0.005359, shape = 0.3938, lower.tail = FALSE)
#	p2 = pexp(c(0.02285041, 0.03072437, 0.05717684, 0.1308201), rate = 94.1699136, lower.tail = FALSE)
