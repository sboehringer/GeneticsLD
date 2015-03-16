#
#	visualization.R
#Fri Mar 21 16:44:54 2014

library('ggplot2');

# <!> copy from Rstatistics.R
vector.std = function(v, C = 1)(C * v / sum(v));

selectSNPsPhysical = function(snp = 'rs9269794', legend, marginNeg = 1e5, marginPos = 1e5) {
	snp = as.character(snp);
	i = which(legend$rs == snp);
	pos = legend$position[i];
	snps = which(legend$position >= pos - marginNeg & legend$position <= pos + marginPos);
	r = legend$rs[snps];
	r
}
selectSNPsOrder = function(snp = 'rs9269794', legend, marginNeg = 5, marginPos = 5) {
	i = which(legend$rs == snp);
	snps = max(1, i - marginNeg):min(nrow(legend), i + marginPos);
	r = legend$rs[snps];
	r
}

snpsAfs = function(snps, hts, legend) {
	r = nlapply(as.character(snps), function(snp) {
		i = which(legend$rs == snp);
		af = table.n(hts[, i], n = 1, min = 0);
		vector.std(af)[1]
	});
	r	
}

snpsMaf = function(snps, hts, legend, maf = 0.05) {
	snpsAfs = snpsAfs(snps, hts, legend);
	snps[snpsAfs > maf & snpsAfs < (1 - maf)];
}

dataSel = function(data, sel) {
	d1 = list(hts = data$hts[, sel, drop = F], legend = droplevels(data$legend[sel, , drop = F]));
	if (!is.null(data$afs)) d1 = c(d1, list(afs = data$afs[sel]));
	d1
}

SnpAf = function(snp, data) {
	af = table.n(data$hts[, which(data$legend$rs == as.character(snp))], n = 1, min = 0);
	vector.std(af)[1]
}
SnpsAfs = function(data)nlapply(data$legend$rs, SnpAf, data = data);

SnpsMaf = function(data, maf = 0.05, retain = rep(T, nrow(data$legend))) {
	snpsAfs = SnpsAfs(data);
	sel = (snpsAfs >= maf & snpsAfs <= (1 - maf)) | retain;
	data1 = dataSel(data, sel);
	data1 = c(data1, list(afs = unlist(snpsAfs[sel])));
	data1
}

SelectSNPsPhysical = function(snp = 'rs9269794', data, marginNeg = 1e5, marginPos = 1e5, maf = .05, buffer = 2) {
	i = which(legend$rs == snp);
	pos = legend$position[i];
	snps = which(legend$position >= pos - marginNeg & legend$position <= pos + marginPos);
	r = legend$rs[snps];
	r
}

#' Read phased hapmap data
#'
#' @param path path prefix
#' @param postfix
#' @export
readPhasedData = function(path, postfixHts = '.phase.gz', postfixLegend = '_legend.txt.gz') {
	hts = read.table(sprintf('%s%s', path, postfixHts));
	legend = read.table(sprintf('%s%s', path, postfixLegend), header = T, stringsAsFactors = T);
	data = list(hts = hts, legend = legend);
	data
}

#' Fetch data matrix of genotypes for SNPs sourrounding a given SNP
#'
#' @param snp name of SNP (rs-number), will be included in the result
#' @param data set as read by snpMatrix
#' @export
SelectSNPsOrder = function(snp = 'rs9269794', data, marginNeg = 5, marginPos = 5, maf = .05, buffer = 2) {
	snp = as.character(snp);	# remove factor status
	# <p> buffered selection
	i = which(data$legend$rs == snp);
	snps = max(1, i - marginNeg * buffer):min(nrow(data$legend), i + marginPos * buffer);
	#print(as.character(data$legend$rs[snps]));
	d1 = dataSel(data, snps);
	# <p> maf filtering
	d2 = SnpsMaf(d1, maf = maf, retain = (1:nrow(d1$legend)) == which(d1$legend$rs == snp) );
	# <p> prune
	i = which(d2$legend$rs == snp);
	snps = max(1, i - marginNeg):min(nrow(d2$legend), i + marginPos);
	d3 = dataSel(d2, snps);
	d3
}

snpCombinations = function(snp, data, N = 3, ..., selectionBy = selectSNPsPhysical) {
	snpI = which(as.character(snp) == data$legend$rs);
	snps = selectionBy(snp, data, ...)$legend$rs;
	# combinations w/o seed SNP
	cbs = combn(setdiff(which.indeces(snps, data$legend$rs), snpI), N - 1);
	cbs = cbind(t(cbs), snpI);	# correct SNP indeces
	r = list(snps = as.character(snps), combinations = cbs);
	r
}

#' Compute haplotype frequencies from matrix of alleles with rows representing haplotypes
#'
#' @param hts 0-1 valued matrix with rows representing haplotypes
#' @export
haplotypeFrequencies = function(hts) {
	N = ncol(hts);
	htsO = apply(hts, 1, bin2ord);
	htCounts = table.n(htsO, n = 2^N - 1, min = 0);
	hfs = vector.std(htCounts);
	hfs
}

haplotypeFrequenciesForSelection = function(sel, data) {
	r = lapply(1:nrow(sel$combinations), function(i) {
		hfs = haplotypeFrequencies(hts[, sel$combinations[i, ]]);
		hfs
	});
	r
}

hfsReparametrize = function(hfs) {
	r = sapply(hfs, function(hf) {
		p1s = par$multinomial2p1sStd(hf);
		p1sM = cbind(matrix(par$p1sMinMax(p1s), ncol = 2), p1s);
		dimnames(p1sM)[[2]][1:2] = paste('p1s', c('min', 'max'), sep = '_');

		pCumu = par$multinomial2cumuStd(hf);
		pCumuM = cbind(matrix(par$cumuMinMax(pCumu), ncol = 2), pCumu);
		dimnames(pCumuM)[[2]][1:2] = paste('pCumu', c('min', 'max'), sep = '_');
		r = cbind(p1sM, pCumuM);
		r[nrow(r), ]
	});
}

ldPlot = function() {
	tr = treeRects(st, heightScale = F);
	eps = 2e-1;
	trPlot = tr[tr$xmax - tr$xmin > eps, , drop = F];
	p = ggplot(trPlot, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label)) + 
		geom_rect() + 
		geom_text(aes(x = xlabel, y = ylabel, label = label), colour = "black", size = 3) +
		coord_polar(start = pi/2) + 
		scale_y_continuous(limits = c(0, max(trPlot$ymax)));
}

parValues = function(hfs) {
	r0 = hfsReparametrize(hfs);
	pCumu = r0['pCumu', ];
	# destandardize to compute d'
	cumuMin = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, min);
	cumuMax = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, max);
	cumuRange = cumuMax - cumuMin;
	cumuN = pCumu * cumuRange + cumuMin;
	dp = ifelse(cumuN > 0, cumuN / cumuMax, cumuN / cumuMin);
	
	r = list(
		minCumu = minimax(min(dp), 0, 1),
		maxCumu = minimax(max(dp), 0, 1)
	);
	r
}

parValuesQ = function(hfs) {
	r0 = hfsReparametrize(hfs);
	maxQts = quantile(r0['pCumu', ], c(.01, .99));
	pCumu = r0['pCumu', ];
	# destandardize to compute d'
	cumuMin = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, min);
	cumuMax = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, max);
	cumuRange = cumuMax - cumuMin;
	cumuN = pCumu * cumuRange + cumuMin;
	dp = ifelse(cumuN > 0, cumuN / cumuMax, cumuN / cumuMin);
	#if (all(dp == 1)) return(list(minCumu = 1, maxCumu = 1));
	if (all(dp == 1)) return(list(minCumu = NA, maxCumu = NA));
	dp = dp[dp > 0 & dp < 1];
	if (!length(dp)) return(list(minCumu = NA, maxCumu = NA));
	r = list(
		minCumu = minimax(min(dp), 0, 1),
		maxCumu = minimax(max(dp), 0, 1)
	);
	r
}

parValuesMean = function(hfs, threshold = .95) {
	r0 = hfsReparametrize(hfs);
	pCumu = r0['pCumu', ];
	# destandardize to compute d'
	cumuMin = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, min);
	cumuMax = apply(r0[c('pCumu_max', 'pCumu_min'), , drop = F], 2, max);
	cumuRange = cumuMax - cumuMin;
	cumuN = pCumu * cumuRange + cumuMin;
	dp = minimax(ifelse(cumuN > 0, cumuN / cumuMax, cumuN / cumuMin), 0, 1);
	dp = dp[dp < 1];
	if (length(dp) == 0) return(list(mean = NA));
	r = list(
		mean = mean(dp > threshold)
	);
	r
}
visualizeSnp = function(snp, data, models, parFunction = parValues, selectionBy = SelectSNPsOrder) {
	print(as.character(snp));
	r0 = iterateModels(models, function(level, window) {
		Log(Sprintf("Level: %{level}d, window: %{window}02d"));
		cbs = snpCombinations(snp, data, N = level,
			marginNeg = window, marginPos = window, selectionBy = selectionBy);
		hfs = haplotypeFrequenciesForSelection(cbs, data);
		parFunction(hfs);
	});
	ns = unlist(iterateModels(models, function(level, window)Sprintf("N:%{level}d,W:%{window}02d"))$results);
	ns1 = names(r0$results[[1]]);
	r1 = unlist.n(r0$results, 1);
	names(r1) = Sprintf("%{name}s:%{par}s", name = ns, par = ns1, sprintf_cartesian = T);
	r1 = c(list(af = unlist(SnpAf(snp, data))), r1);
	r1
}

visualizeSnpDiff = function(snp, data, models, maf = 0.05, selectionBy = SelectSNPsOrder,
	parFunction = parValues) {
	#r1 = visualizeWindow(snp, window, data, models, maf = maf, selectionBy = selectionBy, parFunction = parFunction);
	# <p> compute parameters
	r1 = visualizeSnp(snp, data, models, selectionBy = selectionBy, parFunction = parFunction);
	# <p> compute differences
	mds = iterateModels(models)$models_symbolic;
	lvls = models$level;
	r0 = r1[-1];
	r2 = lapply(2:length(lvls), function(i) {
		lapply(models$window, function(w) {
			j = which(apply(mds, 1, function(r)all(r == c(lvls[i], w))));
			jP = which(apply(mds, 1, function(r)all(r == c(lvls[i - 1], w))));
			r0[[2*j - 1]] - r0[[2*jP - 1]]
		});
	});
	r2 = unlist.n(r2, 1);
	names(r2) = Sprintf("N:%{level}d,W:%{window}02d", level = lvls[-1], window = models$window);
	r2
}

visualizeWindow = function(centerSnp, window = 3e4, data, models, maf = 0.05,
	selectionBy = SelectSNPsOrder, parFunction = parValues, visualizer = visualizeSnp) {
	d0 = selectionBy(centerSnp, data, marginNeg = window, marginPos = window, maf = maf);
	snps = as.character(d0$legend$rs);
	buffer = window + max(models$window);
	d1 = selectionBy(centerSnp, data, marginNeg = buffer, marginPos = buffer, maf = maf);
#snps = snps[1:3];
	#snps = "rs7749092";
	print(snps);
	r = nlapply(snps, visualizer, models = models, data = d1,
		selectionBy = selectionBy, parFunction = parFunction);
}

visualizeToPath = function(centerSnp, modelList, data, output, parFunction = parValues, window = 30, maf = 0.05, visualizer = visualizeSnp) {
	print(iterateModels(modelList)$models);
	v = visualizeWindow(centerSnp, window = window, data, modelList, maf = maf,
		parFunction = parFunction,
		selectionBy = SelectSNPsOrder,
		visualizer = visualizer);
	vM = sapply(v, identity);

	vmD = reshape.long(data.frame(vM), factorColumn = 'x', rowNamesAs = 'y');
	p = ggplot(vmD, aes(x, y, fill = value)) + geom_raster() +
		scale_fill_gradient2(low = "blue", mid = "white", high = "red",
			midpoint = 0.5, space = "rgb", na.value = "grey50", guide = "colourbar") +
		opts(axis.text.x = theme_text(angle=-90)) +
		opts(axis.text.x=theme_text(angle=-90)) + xlab(NULL) + ylab(NULL);
	ggsave(p, file = output, width = 14, height = 8);
}

#' Visualize a matrix of values
#'
#' @param mat matrix with values
#' @export
visualizeMatrix = function(mat) {
	vmD = reshape.long(data.frame(mat), factorColumn = 'x', rowNamesAs = 'y');
	p = ggplot(vmD, aes(x, y, fill = value)) + geom_raster() +
		scale_fill_gradient2(low = "blue", mid = "white", high = "red",
			midpoint = 0.5, space = "rgb", na.value = "grey50", guide = "colourbar") +
		opts(axis.text.x = theme_text(angle=-90)) +
		opts(axis.text.x=theme_text(angle=-90)) + xlab(NULL) + ylab(NULL);
	p
}
