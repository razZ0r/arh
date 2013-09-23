# author: Axel Rasche
# 
# This script is a minimal implementation of the method comparison in 
# "ARH: Predicting Splice Variants from Genome-wide Data with Modified 
# Entropy"; Rasche and Herwig (2009). All methods are implemented in 
# this script with subsequent evaluation. The implementation is 
# focussed on the generation of figure 2.B in the liver vs. pancreas 
# test case.
# 
# Input:  load preprocessed data and annotation as matrices
# Output: Performance evaluation as ROC (figure) and AUC (table)


library(ROCR)



######################################################################
### Input
######################################################################

chip_treat = c("huex_wta_liver_A.CEL", "huex_wta_liver_B.CEL", "huex_wta_liver_C.CEL")
chip_ctrl = c("huex_wta_pancreas_A.CEL", "huex_wta_pancreas_B.CEL", "huex_wta_pancreas_C.CEL")

load(file = "preproc.RData")
	# load "probe_int", "probe_map"
	# probe_map is a simple matrix for the mapping of probe ID to exon ID or gene ID
	# In probe_int are the normalised probe intensities, 
	# Rows correspond to probes and requires to correspond to the rows of probe_map.
	# Columns correspond to chips. Column names have to correspond to chip_treat and chip_ctrl.
#probe_map = probe_map[ 1:10000 ,  ]
#probe_int = probe_int[ 1:10000 ,  ]
#> head(probe_int)
#        huex_wta_liver_A.CEL huex_wta_liver_B.CEL huex_wta_liver_C.CEL
#3426751             69.05053             69.08707             80.88705
#5039551             77.00563             55.61893             69.14864
#2172317             95.35908             80.93277             82.14248
#2179997            106.81056             88.52317            104.02855
#5160184             33.95972             35.32264             27.53149
#5316344             33.76084             37.36830             36.46284
#        huex_wta_pancreas_A.CEL huex_wta_pancreas_B.CEL huex_wta_pancreas_C.CEL
#3426751                52.12755                49.41366                39.10893
#5039551                46.33757                48.81233                59.63601
#2172317                59.40152                77.64930               101.17854
#2179997                72.10185                69.18507                60.44979
#5160184                17.46553                21.71740                37.07784
#5316344                20.69318                17.14619                17.30234
#> dim(probe_int)
#[1] 1434190       6
#> class(probe_int)
#[1] "matrix"
#> mode(probe_int)
#[1] "numeric"
#> head(probe_map)
#        probe_id  exon_id           gene_id          
#3426751 "3426751" "ENSE00000328462" "ENSG00000092377"
#5039551 "5039551" "ENSE00000328462" "ENSG00000092377"
#2172317 "2172317" "ENSE00000328462" "ENSG00000092377"
#2179997 "2179997" "ENSE00000328462" "ENSG00000092377"
#5160184 "5160184" "ENSE00000328462" "ENSG00000092377"
#5316344 "5316344" "ENSE00000328462" "ENSG00000092377"
#> dim(probe_map)
#[1] 1434190       3
#> class(probe_map)
#[1] "matrix"
#> mode(probe_map)
#[1] "character"
#


######################################################################
## Functions

#---------------------------------------------------------------------
lbind = function(x = "list", y = "list") {
	if(length(x) != length(y)) { stop("ERROR in lbind, x and y must have same length") }
	 	Llist = lapply(X = 1:length(x), FUN = function(PartNr, x, y) {
		return(list(x[[ PartNr ]], y[[ PartNr ]]))
	}, x = x, y = y)
	if(identical(names(x), names(y))) { names(Llist) = names(x) }
 	return(Llist)
}

#---------------------------------------------------------------------
stretch = function(x = "ANY") {
	new = x[ match(x = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "gene_id" ], table = names(x)) ]
	names(new) = unique(probe_map[  , "exon_id" ])
	return(new)
}

# --------------------------------------------------------------------
setMethod("split", signature(x = "matrix", f = "ANY", drop = "missing"), definition = function(x = "matrix", f = "ANY", drop = FALSE) {
	Split = split(x = 1:nrow(x), f = f, drop = drop)
	L = lapply(X = Split, FUN = function(X) {
		return(x[ X ,  ])
	})
	return(L)
})

# --------------------------------------------------------------------
Score_prepare = function(X) {
	# X is a filter
	if(class(X) == "logical") { 
		add = rep(NA, sum(!(names(Reference) %in% names(X))))
		names(add) = names(Reference)[ !(names(Reference) %in% names(X)) ]
		X = c(X, add)
		rm(add)
		X = X[ match(x = names(Reference), table = names(X)) ]
		Fine = !is.na(X) & X & is.finite(Reference)
		X[ !is.na(X) & X ] = Reference[ !is.na(X) & X ]
		X[ !(!is.na(X) & X) ] = NA
	} else {
		Fine = !is.na(X) & is.finite(X)
	}
	
	# NA, Inf, missing predictions on behind of X
	X[ !Fine ] = max(X[ Fine ]) + 1
	if(all(names(X)[ 1:5 ] %in% probe_map[  , "exon_id" ])) {
		X[ !Fine & names(X) %in% TP_exon ] = max(X[ Fine ]) + 2
	} else if(all(names(X)[ 1:5 ] %in% probe_map[  , "gene_id" ])) {
		X[ !Fine & names(X) %in% TP_gene ] = max(X[ Fine ]) + 2
	} else { stop("predictions for exons or genes!") }
	X = sort(X)

	return(X)
}



######################################################################
### stat evaluation
######################################################################

Methods = vector("list")
# Collect results of different methods, numeric vectors of positive 
# values. Small values indicate higher splicing indication.


######################################################################
## ARH, gene level splicing prediction, outer turbo (IQR-Q for corrected entropy)
# --------------------------------------------------------------------
ARH_f = function(x = "numeric", y = "numeric", f = "character", na.rm = FALSE) {
# INPUT
#	x, y: exon expression; numeric vectors, x,y >= 0
#	f : gene IDs corresponding to x,y; groups the exon expressions to genes, character
#	x, y, f all of the same length
#	na.rm: how do you want to treat NA in the expressions? R standard is here implemented, but 'na.rm = TRUE' also is reasonable.
# OUTPUT
#	ARH values, numeric vector, length(ARH) == length(unique(f)), names(ARH) == names(unique(f))

	# negative exon expressions not interpretable
	if(any(x < 0, na.rm = TRUE) | any(y < 0, na.rm = TRUE)) stop("What are negative expressions?")
	
	# avoid division by zero and arbitrarily big numbers
	y[ y < 0.0001 ] = 0.0001

	# compute splicing probabilities, each exon compared to gene
	splicingProbabilities = log2(x / y)
	splicingProbabilitiesMedian = split(splicingProbabilities, f)
	splicingProbabilitiesMedian = sapply(X = splicingProbabilitiesMedian, FUN = median, na.rm = na.rm)
	splicingProbabilitiesMedian = splicingProbabilitiesMedian[ match(f, names(splicingProbabilitiesMedian)) ]
	splicingProbabilities = 2^abs( splicingProbabilities - splicingProbabilitiesMedian )
	rm(splicingProbabilitiesMedian)
	splicingProbabilitiesSum = split(splicingProbabilities, f)
	splicingProbabilitiesSum = sapply(splicingProbabilitiesSum, sum, na.rm = na.rm)
	splicingProbabilitiesSum = splicingProbabilitiesSum[ match(f, names(splicingProbabilitiesSum)) ]
	splicingProbabilities = splicingProbabilities / splicingProbabilitiesSum
	rm(splicingProbabilitiesSum)
	
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
Entropy_treat = split(x = probe_int[  , chip_treat ], f = probe_map[  , "exon_id" ])
Entropy_treat = Entropy_treat[ match(x = unique(probe_map[  , "exon_id" ]), table = names(Entropy_treat)) ]
Entropy_treat = sapply(Entropy_treat, median, na.rm = TRUE)
Entropy_ctrl = split(x = probe_int[  , chip_ctrl ], f = probe_map[  , "exon_id" ])
Entropy_ctrl = Entropy_ctrl[ match(x = unique(probe_map[  , "exon_id" ]), table = names(Entropy_ctrl)) ]
Entropy_ctrl = sapply(Entropy_ctrl, median, na.rm = TRUE)
Entropy = ARH_f(x = Entropy_treat, y = Entropy_ctrl, f = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "gene_id" ], na.rm = TRUE)
rm(Entropy_treat, Entropy_ctrl)

ARH = 1/Entropy
Methods = c(Methods, ARH = list(ARH))
rm(Entropy, ARH_f, ARH)


######################################################################
## Splicing Index
gene_xpr_treat = split(probe_int[  , chip_treat ], f = probe_map[  , "gene_id" ])
gene_xpr_treat = gene_xpr_treat[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_treat)) ]
gene_xpr_treat = sapply(gene_xpr_treat, median)
gene_xpr_ctrl = split(probe_int[  , chip_ctrl ], f = probe_map[  , "gene_id" ])
gene_xpr_ctrl = gene_xpr_ctrl[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_ctrl)) ]
gene_xpr_ctrl = sapply(gene_xpr_ctrl, median)
exon_xpr_treat = split(probe_int[  , chip_treat ], f = probe_map[  , "exon_id" ])
exon_xpr_treat = exon_xpr_treat[ match(x = unique(probe_map[  , "exon_id" ]), table = names(exon_xpr_treat)) ]
exon_xpr_treat = sapply(exon_xpr_treat, median)
exon_xpr_ctrl = split(probe_int[  , chip_ctrl ], f = probe_map[  , "exon_id" ])
exon_xpr_ctrl = exon_xpr_ctrl[ match(x = unique(probe_map[  , "exon_id" ]), table = names(exon_xpr_ctrl)) ]
exon_xpr_ctrl = sapply(exon_xpr_ctrl, median)

SplicingIndex = log2(exon_xpr_treat) - log2(exon_xpr_ctrl) + stretch(log2(gene_xpr_ctrl)) - stretch(log2(gene_xpr_treat))
SplicingIndex = 1/abs(SplicingIndex)
Methods = c(Methods, SplicingIndex = list(SplicingIndex))
rm(SplicingIndex, gene_xpr_treat, gene_xpr_ctrl, exon_xpr_treat, exon_xpr_ctrl)


######################################################################
## SPLICE
SPLICE_probe_treat = probe_int[  , chip_treat ]
if(NCOL(SPLICE_probe_treat) > 1) { SPLICE_probe_treat = rowMeans(SPLICE_probe_treat) }
SPLICE_gene_xpr_treat = split(SPLICE_probe_treat, f = probe_map[  , "gene_id" ])
SPLICE_gene_xpr_treat = SPLICE_gene_xpr_treat[ match(x = unique(probe_map[  , "gene_id" ]), table = names(SPLICE_gene_xpr_treat)) ]
SPLICE_gene_xpr_treat = sapply(SPLICE_gene_xpr_treat, median)
SPLICE_gene_xpr_treat = SPLICE_gene_xpr_treat[ match(x = probe_map[  , "gene_id" ], table = names(SPLICE_gene_xpr_treat)) ]
SPLICE_probe_ctrl = probe_int[  , chip_ctrl ]
if(NCOL(SPLICE_probe_ctrl) > 1) { SPLICE_probe_ctrl = rowMeans(SPLICE_probe_ctrl) }
SPLICE_gene_xpr_ctrl = split(SPLICE_probe_ctrl, f = probe_map[  , "gene_id" ])
SPLICE_gene_xpr_ctrl = SPLICE_gene_xpr_ctrl[ match(x = unique(probe_map[  , "gene_id" ]), table = names(SPLICE_gene_xpr_ctrl)) ]
SPLICE_gene_xpr_ctrl = sapply(SPLICE_gene_xpr_ctrl, median)
SPLICE_gene_xpr_ctrl = SPLICE_gene_xpr_ctrl[ match(x = probe_map[  , "gene_id" ], table = names(SPLICE_gene_xpr_ctrl)) ]

SPLICE = log2(SPLICE_probe_treat) - log2(SPLICE_probe_ctrl) + log2(SPLICE_gene_xpr_ctrl) - log2(SPLICE_gene_xpr_treat)
rm(SPLICE_probe_treat, SPLICE_gene_xpr_treat, SPLICE_probe_ctrl, SPLICE_gene_xpr_ctrl)
SPLICE = split(SPLICE, f = probe_map[  , "exon_id" ])
SPLICE = SPLICE[ match(x = unique(probe_map[  , "exon_id" ]), table = names(SPLICE)) ]
SPLICE = 1/abs(sapply(SPLICE, median))

Methods = c(Methods, SPLICE = list(SPLICE))
rm(SPLICE)


######################################################################
## PAC
gene_xpr_treat = split(probe_int[  , chip_treat ], f = probe_map[  , "gene_id" ])
gene_xpr_treat = gene_xpr_treat[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_treat)) ]
gene_xpr_treat = sapply(gene_xpr_treat, median)
exon_xpr_treat = split(probe_int[  , chip_treat ], f = probe_map[  , "exon_id" ])
exon_xpr_treat = exon_xpr_treat[ match(x = unique(probe_map[  , "exon_id" ]), table = names(exon_xpr_treat)) ]
exon_xpr_treat = sapply(exon_xpr_treat, median)
gene_xpr_PAC = split(probe_int, f = probe_map[  , "gene_id" ])
gene_xpr_PAC = gene_xpr_PAC[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_PAC)) ]
gene_xpr_PAC = sapply(gene_xpr_PAC, median)
exon_xpr_PAC = split(probe_int, f = probe_map[  , "exon_id" ])
exon_xpr_PAC = exon_xpr_PAC[ match(x = unique(probe_map[  , "exon_id" ]), table = names(exon_xpr_PAC)) ]
exon_xpr_PAC = sapply(exon_xpr_PAC, median)

PAC = exon_xpr_treat - stretch(gene_xpr_treat) * (exon_xpr_PAC / stretch(gene_xpr_PAC))
PAC = 1/abs(PAC)
Methods = c(Methods, PAC = list(PAC))
rm(PAC, gene_xpr_treat, exon_xpr_treat, gene_xpr_PAC, exon_xpr_PAC)


######################################################################
## ANOSVA
Intensity = split(x = log2(probe_int), f = probe_map[  , "gene_id" ])
Intensity = Intensity[ match(x = unique(probe_map[  , "gene_id" ]), table = names(Intensity)) ]
Intensity = lapply(X = Intensity, as.numeric)
#table(sapply(Intensity, length))

Exon = split(x = rep(probe_map[  , "exon_id" ], times = ncol(probe_int)), f = probe_map[  , "gene_id" ])
Exon = Exon[ match(x = unique(probe_map[  , "gene_id" ]), table = names(Exon)) ]
Exon = lapply(Exon, factor)
ExonNr = sapply(Exon, function(X) length(unique(X)))

Condition = split(x = c(rep("treat", times = length(chip_treat) * nrow(probe_int)), rep("ctrl", times = length(chip_ctrl) * nrow(probe_int))), f = probe_map[  , "gene_id" ])
Condition = Condition[ match(x = unique(probe_map[  , "gene_id" ]), table = names(Condition)) ]
Condition = lapply(Condition, factor)

Data = lbind(Exon, Condition)
Data = lapply(X = 1:length(Data), FUN = function(X, Intensity, Data) {
	O = c(Intensity[ X ], Data[[ X ]])
	names(O) = c("y", "exon", "condition")
	return(O)
}, Intensity = Intensity, Data = Data)
names(Data) = unique(probe_map[  , "gene_id" ])
rm(Intensity, Exon, Condition)

ANOSVA = rep(NA, times = length(Data))
ANOSVA[ ExonNr > 1 ] = sapply(X = Data[ ExonNr > 1 ], FUN = function(X) {
	return(anova(lm(y ~ 1 + exon + condition + exon:condition, data = X))[ "exon:condition" , "Pr(>F)" ])
})
names(ANOSVA) = names(Data)
rm(Data, ExonNr)

ANOSVA_log2 = ANOSVA
Methods = c(Methods, ANOSVA = list(ANOSVA_log2))
rm(ANOSVA, ANOSVA_log2)


######################################################################
## MiDAS
# processed by Affymetrix Power Tools (apt), install from Affymetrix-homepage
dir.create(path = "MiDAS", showWarnings = FALSE, recursive = FALSE)

Norm.exon.matrix = apply(X = probe_int, MARGIN = 2, FUN = function(X) {
	S = split(X, probe_map[  , "exon_id" ])
	S = S[ match(x = unique(probe_map[  , "exon_id" ]), table = names(S)) ]
	return(sapply(S, median))
})
Norm.gene.matrix = apply(X = probe_int, MARGIN = 2, FUN = function(X) {
	S = split(X, probe_map[  , "gene_id" ])
	S = S[ match(x = unique(probe_map[  , "gene_id" ]), table = names(S)) ]
	return(sapply(S, median))
})

# write files for apt
exonR = 1:length(unique(probe_map[  , "exon_id" ])) + 1000000
names(exonR) = unique(probe_map[  , "exon_id" ])
geneR = 1:length(unique(probe_map[  , "gene_id" ])) + 2000000
names(geneR) = unique(probe_map[  , "gene_id" ])
probe_mapR = probe_map
probe_mapR[  , "exon_id" ] = exonR[ probe_mapR[  , "exon_id" ] ]
probe_mapR[  , "gene_id" ] = geneR[ probe_mapR[  , "gene_id" ] ]

Chips.m = cbind(c(chip_treat, chip_ctrl), c(rep("treat", length(chip_treat)), rep("ctrl", length(chip_ctrl))))
colnames(Chips.m) = c("cel_files", "group_id")
write.table(Chips.m, file = "MiDAS/cels.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

Extract = probe_mapR[ !duplicated(probe_mapR[  , "exon_id" ]) , c("exon_id", "gene_id") ]
Extract = split(Extract[  , "exon_id" ], Extract[  , "gene_id" ])
Extract = Extract[ match(x = unique(probe_mapR[  , "gene_id" ]), table = names(Extract)) ]
Extract = sapply(Extract, paste, collapse = " ")
Extract = cbind(names(Extract), Extract)
colnames(Extract) = c("probeset_id", "probeset_list")
write.table(Extract, file = "MiDAS/metaprobeset.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(Extract)

N = cbind(geneR[ rownames(Norm.gene.matrix) ], Norm.gene.matrix)
colnames(N)[ 1 ] = "probeset_id"
write.table(N, file = "MiDAS/gene.summary.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
N = cbind(exonR[ rownames(Norm.exon.matrix) ], Norm.exon.matrix)
colnames(N)[ 1 ] = "probeset_id"
write.table(N, file = "MiDAS/exon.summary.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(N)

# call apt
WD = getwd()
setwd("MiDAS")
system("/project/hfa_work/Rasche/Programme/apt-1.8.0/bin/apt-midas --cel-files cels.txt -g gene.summary.txt -e exon.summary.txt -m metaprobeset.txt")
setwd(WD)
rm(WD)

# read results from apt
MiDAS.r = read.table(file = "MiDAS/midas.pvalues.txt", quote = "", header = TRUE)
MiDAS.r = as.matrix(MiDAS.r)[  , c("probeset_list_id", "pvalue") ]
MiDAS = MiDAS.r[  , c("pvalue") ]
exonRR = names(exonR)
names(exonRR) = exonR
names(MiDAS) = exonRR[ as.character(MiDAS.r[  , c("probeset_list_id") ]) ]
MiDAS = MiDAS[ match(x = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "exon_id" ], table = names(MiDAS)) ]
Methods = c(Methods, MiDAS = list(MiDAS))

rm(MiDAS, MiDAS.r, Norm.exon.matrix, Norm.gene.matrix, exonR, geneR, exonRR)
system(paste("rm -rf ", "MiDAS", sep = ""))


######################################################################
## FIRMA
# FIRMA uses the aroma.affymetrix package with some own annotation files. 
# Please refer to the aroma.affymetrix package and manuals for a local installation.
library(aroma.affymetrix)
path_CEL = "/project/altsplice/projekte/tissue_human_AEdb/cel/"
	# path to the respective CEL files for liver and pancreas
dir.create(path = "FIRMA", showWarnings = FALSE, recursive = FALSE)
WD = getwd()
setwd("FIRMA")
	system(paste("cp -r /project/altsplice/etc/FIRMA/HuEx-1_0-st-v2,U-Ensembl49,G-Affy/annotationData .", sep = ""))
	chipType <- "HuEx-1_0-st-v2"
	cdf <- AffymetrixCdfFile$byChipType(chipType, tags="U-Ensembl49,G-Affy")
dir.create(path = "rawData/HuEx-1_0-st-v2/", showWarnings = FALSE, recursive = TRUE)
sapply(X = paste("ln -s ", path_CEL, c(chip_treat, chip_ctrl), " rawData/HuEx-1_0-st-v2/", c(chip_treat, chip_ctrl), sep = ""), FUN = system)
cs <- AffymetrixCelSet$fromFiles("rawData/HuEx-1_0-st-v2", cdf=cdf)
setCdf(cs,cdf)
	bc <- RmaBackgroundCorrection(cs)
	#bc <- RmaBackgroundCorrection(cs, tag="coreR3")
csBC <- process(bc)
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn)
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
fit(plmTr)
rs <- calculateResiduals(plmTr)
cesTr <- getChipEffectSet(plmTr)
cesEx <- getChipEffectSet(plmEx)
firma <- FirmaModel(plmTr)
fit(firma)
fs <- getFirmaScores(firma)
out = extractDataFrame(fs, addNames = TRUE)
setwd(WD)
rm(WD)
rownames(out) = gsub(pattern = "_at", replacement = "", x = out[  , "groupName" ])
colnames(out)[ 6:length(colnames(out)) ] = paste(colnames(out)[ 6:length(colnames(out)) ], ".CEL", sep = "")
F = out[  , "unitName" ]
out = out[  , 6:length(colnames(out)) ]
out = as.matrix(out)
out = split(out[  , c(chip_treat, chip_ctrl) ], F)
rm(F)
out = lapply(X = out, FUN = function(X) {
	if(NCOL(X) == 1) return(X)
	ret = apply(X = X, MARGIN = 2, FUN = function(Y) {
		return(Y[ which(abs(log2(Y)) == max(abs(log2(Y)))) ])
	})
	return(ret)
})
out = do.call("rbind", out)
M = unique(probe_map[  , "gene_id" ])[ !(unique(probe_map[  , "gene_id" ]) %in% rownames(out)) ]
N = matrix(NA, nrow = length(M), ncol = ncol(out))
mode(N) = "numeric"
rownames(N) = M
colnames(N) = colnames(out)
out = rbind(out, N)
rm(M, N)
out = out[ match(x = unique(probe_map[  , "gene_id" ]), table = rownames(out)) ,  ]
FIRMA = out[  , c(chip_treat, chip_ctrl) ]
rm(chipType, cdf, cs, bc, csBC, qn, csN, plmTr, rs, cesTr, firma, fs, out, path_CEL)
system(paste("rm -rf FIRMA", sep = ""))

FIRMA_log =      1/abs(log2(rowMeans(FIRMA[  , chip_treat ], na.rm = TRUE)))
Methods = c(Methods, FIRMA = list(FIRMA_log))
rm(FIRMA, FIRMA_lin, FIRMA_log, FIRMA_diff_log, FIRMA_log_diff, FIRMA_diff_lin)


######################################################################
## MADS
# (1)
gene_xpr_treat = split(probe_int[  , chip_treat ], f = probe_map[  , "gene_id" ])
gene_xpr_treat = gene_xpr_treat[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_treat)) ]
gene_xpr_treat = sapply(gene_xpr_treat, median)
gene_xpr_treat = gene_xpr_treat[ match(x = probe_map[  , "gene_id" ], table = names(gene_xpr_treat)) ]
gene_xpr_ctrl = split(probe_int[  , chip_ctrl ], f = probe_map[  , "gene_id" ])
gene_xpr_ctrl = gene_xpr_ctrl[ match(x = unique(probe_map[  , "gene_id" ]), table = names(gene_xpr_ctrl)) ]
gene_xpr_ctrl = sapply(gene_xpr_ctrl, median)
gene_xpr_ctrl = gene_xpr_ctrl[ match(x = probe_map[  , "gene_id" ], table = names(gene_xpr_ctrl)) ]
Probe_int = probe_int
Probe_int[  , chip_treat ] = apply(X = Probe_int[  , chip_treat ], MARGIN = 2, FUN = function(X) return(X/gene_xpr_treat))
Probe_int[  , chip_ctrl ] = apply(X = Probe_int[  , chip_ctrl ], MARGIN = 2, FUN = function(X) return(X/gene_xpr_ctrl))
rm(gene_xpr_treat, gene_xpr_ctrl)

# (2)
ttest_g = apply(X = Probe_int, MARGIN = 1, FUN = function(X) {
	if(var(X) < 0.00001) return(1)
	return(t.test(x = as.numeric(X[ chip_treat ]), y = as.numeric(X[ chip_ctrl ]), var.equal = TRUE, alternative = "greater")[[ "p.value" ]])
})
ttest_l = apply(X = Probe_int, MARGIN = 1, FUN = function(X) {
	if(var(X) < 0.00001) return(1)
	return(t.test(x = as.numeric(X[ chip_treat ]), y = as.numeric(X[ chip_ctrl ]), var.equal = TRUE, alternative = "less")[[ "p.value" ]])
})
rm(Probe_int)

# (3)
overallpvalue<-function(p) {
	#overall p-value calculated from chi-square distribution from a list of p-values
	n=length(p)
	sum=0
	for (i in 1:n) {
		sum=sum-2*log(p[i])
		p_all=pchisq(sum,lower.tail=FALSE, df=n*2)
	}
	return(p_all)
}
ttest_g = split(ttest_g, f = probe_map[  , "exon_id" ])
ttest_g = ttest_g[ match(x = unique(probe_map[  , "exon_id" ]), table = names(ttest_g)) ]
ttest_g = sapply(ttest_g, overallpvalue)
names(ttest_g) = unique(probe_map[  , "exon_id" ])
ttest_l = split(ttest_l, f = probe_map[  , "exon_id" ])
ttest_l = ttest_l[ match(x = unique(probe_map[  , "exon_id" ]), table = names(ttest_l)) ]
ttest_l = sapply(ttest_l, overallpvalue)
names(ttest_l) = unique(probe_map[  , "exon_id" ])
rm(overallpvalue)

# (4)
MADS = pmin(ttest_g, ttest_l)
rm(ttest_g, ttest_l)

Methods = c(Methods, MADS = list(MADS))
rm(MADS)


######################################################################
## correlation
Correlation_treat = split(x = probe_int[  , chip_treat ], f = probe_map[  , "exon_id" ])
Correlation_treat = Correlation_treat[ match(x = unique(probe_map[  , "exon_id" ]), table = names(Correlation_treat)) ]
Correlation_treat = sapply(Correlation_treat, median, na.rm = TRUE)
Correlation_treat = split(Correlation_treat, f = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "gene_id" ])
Correlation_treat = Correlation_treat[ match(x = unique(probe_map[  , "gene_id" ]), table = names(Correlation_treat)) ]
Correlation_ctrl = split(x = probe_int[  , chip_ctrl ], f = probe_map[  , "exon_id" ])
Correlation_ctrl = Correlation_ctrl[ match(x = unique(probe_map[  , "exon_id" ]), table = names(Correlation_ctrl)) ]
Correlation_ctrl = sapply(Correlation_ctrl, median, na.rm = TRUE)
Correlation_ctrl = split(Correlation_ctrl, f = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "gene_id" ])
Correlation_ctrl = Correlation_ctrl[ match(x = unique(probe_map[  , "gene_id" ]), table = names(Correlation_ctrl)) ]
Correlation = lapply(X = 1:length(Correlation_treat), function(X) return(list(Correlation_treat[[ X ]], Correlation_ctrl[[ X ]])))
names(Correlation) = names(Correlation_treat)
rm(Correlation_treat, Correlation_ctrl)
Correlation = 1.01 + sapply(Correlation, function(X) {
	if(length(X[[ 1 ]]) < 3) return(NA)
	return(cor(X[[ 1 ]], X[[ 2 ]]))
})

Methods = c(Methods, Correlation = list(Correlation))
rm(Correlation)


######################################################################
## gene length
GES = probe_map[ !duplicated(probe_map[  , "exon_id" ]) , "gene_id" ]
GES = table(GES)
GES = GES[ match(x = unique(probe_map[  , "gene_id" ]), table = names(GES)) ]
GES_ = as.numeric(GES)
names(GES_) = unique(probe_map[  , "gene_id" ])
GES = GES_
rm(GES_)

NoExonDecrease = 1/GES
Methods = c(Methods, GeneLength = list(NoExonDecrease))
rm(GES, NoExonDecrease)



######################################################################
### Evaluation
######################################################################

######################################################################
## Methods
Lost = sapply(Methods, function(X) { !all(is.na(X[ 1:100 ])) })
Methods = Methods[ Lost ]
if(sum(!Lost) >= 1) { print(paste(sum(Lost), " predictions lost because only NA: ", paste(names(Methods)[ !Lost ], collapse = ", "), sep = "")) }

if(any(sapply(Methods, function(X) min(X, na.rm = TRUE)) < 0)) { stop("only predictions x > 0 allowed") }


######################################################################
## AS true positives
TP_gene = c(
	"ENSG00000100429", "ENSG00000105325", "ENSG00000005471", "ENSG00000010932",
	"ENSG00000131979", "ENSG00000135447", "ENSG00000197965", "ENSG00000143257",
	"ENSG00000170632", "ENSG00000171105", "ENSG00000163606", "ENSG00000082701",
	"ENSG00000183337", "ENSG00000101076", "ENSG00000015475", "ENSG00000142192",
	"ENSG00000148584", "ENSG00000106633"
)
TP_exon = c(
	"ENSE00000657605", "ENSE00000664434", "ENSE00000664435", "ENSE00000703613",
	"ENSE00000789744", "ENSE00000867228", "ENSE00000920061", "ENSE00000920062",
	"ENSE00000958356", "ENSE00001002783", "ENSE00001041739", "ENSE00001132747",
	"ENSE00001157509", "ENSE00001199958", "ENSE00001236373", "ENSE00001245499",
	"ENSE00001309242", "ENSE00001364818", "ENSE00001374791", "ENSE00001377541",
	"ENSE00001435515", "ENSE00001435651", "ENSE00001446058", "ENSE00001462105",
	"ENSE00001521926", "ENSE00001521954", "ENSE00001522887"
)


######################################################################
## prepare
Scores = lapply(X = Methods, FUN = Score_prepare)

Scores_pred = lapply(X = Scores, FUN = function(X) { 
	if(any(c(TP_exon, TP_gene) %in% names(X))) {
		if(all(names(X)[ 1:5 ] %in% probe_map[  , "exon_id" ])) {
			return(prediction(-X, as.numeric(names(X) %in% TP_exon)))
		} else if(all(names(X)[ 1:5 ] %in% probe_map[  , "gene_id" ])) {
			return(prediction(-X, as.numeric(names(X) %in% TP_gene)))
		} else { stop("predictions for exons or genes!") }
	} else {
		return(prediction(-X, c(rep(0, length(X)-(sum(is.infinite(X))+1)), rep(1, sum(is.infinite(X))+1))))
	}
})

Color = c("black", rainbow(length(Scores)-1, start = 2/3))
names(Color) = names(Scores)
LineType = c("solid", "42", "1242", "121242", "12121242", "124242", "12124242", "82", "1282", "11")
names(LineType) = names(Scores)


######################################################################
## AUC
O = lapply(Scores_pred, function(X) { return(performance(X, "auc")) })
AUC = sapply(O, function(X) { return(X@y.values[[ 1 ]]) })
write.table(AUC, file = "AUC.txt", quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
rm(O, AUC)


######################################################################
## ROC, tpr-fpr
O = lapply(Scores_pred, function(X) { return(performance(X, "tpr", "fpr")) })
pdf("TprFpr.pdf")
	plot(O[[ 1 ]], 
		col = Color[ 1 ], 
		lty = LineType[ 1 ], 
		main = "human tissue xprs, liver vs. pancreas, ROC/AEdb"
	)
	for(s in 2:length(O)) { 
		plot(O[[ s ]], 
			col = Color[ s ], 
			lty = LineType[ s ], 
			add = TRUE
		) 
	}
	legend("bottomright", legend = names(Scores), title = "legend", lwd = 1, col = Color, text.col = Color, lty = LineType)
dev.off()
rm(O)


print("DONE")
