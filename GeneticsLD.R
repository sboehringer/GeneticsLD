#
#	GeneticsLD.R
#Fri Feb  6 17:31:47 CET 2015

library('devtools');
library('roxygen2');
#source('GeneticsLD/R/Rdata.R')
#source('GeneticsLD/R/visualization.R')
library('Rcpp');

if (F) {
	#system('rm GeneticsLD/src/*.o GeneticsLD/src/*.so');

	roxygenise('GeneticsLD', roclets = 'rd');
	roxygenize('.', roclets = 'collate');
	roxygenise('GeneticsLD', roclets = 'namespace');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsLD', threads = 6);
}

if (F) {
	library('GeneticsLD');
	#M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
}


if (F) {
	source('RgenericAll.R');
	library('sets');
	library('haplo.stats');
	parallelize_setEnable(F);
}

if (0) {
	# count of 1's in standardized parametrization
	# one parameter vector per row
	hfss = list(
		c( .1, 0, .3, .05, .5, .1, 0, 0),
		c( 0, .1, .3, .05, .5, .1, 0, 0),
		c( 0, .1, 0, .05, .5, .1, 0, 0),
		c( 0, .1, .05, 0, .5, .1, 0, 0),
		c( 0, .1, 0, 0, .5, .1, .05, 0),
		c( 0, .1, 0, 0, .5, .1, 0, .05),
		c( 0, .1, 0, 0, 0, .1, 0, .05),
		c( 0, 0, 0, 0, 0, .1, .1, .05),
		c( 0, .1, .05, 0, .5, 0, 0, .1),
		c( .05, .1, .05, 0, .5, 0, 0, .1)
	);
	M = log2(length(hfss[[1]]));
	ns = apply(ord2bin(0:(length(hfss[[1]]) - 1), digits = M), 1, function(e)join(rev(e), ''));
	singleLoci = sapply(1:M, function(i)2^(i-1));
	eps = 1e-8;
	r = sapply(hfss, function(hfs) {
		p1s = par$multinomial2p1s(vector.std(hfs));
		std = par$multinomial2p1sStd(vector.std(hfs));
		p1sM = cbind(matrix(par$p1sMinMax(p1s), ncol = 2), p1s);
print(p1sM);
		stdP = std[-1]
		stdML = stdP[-(singleLoci)];	# standardized 1s minus single loci

		countH = sum(hfs < eps);
		countE = sum(abs(stdP - 1) < eps | abs(stdP) < eps);
		countEML = sum(abs(stdML -1) < eps | abs(stdML) < eps);
		c(countH, countE, countEML, std, p1s)
	});
	rDf = Df(t(r), names = c('#h', '#', '#S', ns, ns));
	print(round(rDf, 2));
}

std = function(v)(v/sum(v))
hfs2cumu = function(hfs) {
	hfs = std(hfs);
	P = new(parametrizer);
	cumu = P$multinomial2cumu(hfs);
	cumuStd = P$multinomial2cumuStd(hfs);
	rbind(hfs, cumu, cumuStd)
}

if (1) {
	require('GeneticsLD');
	print(hfs2cumu(1:4));
	print(hfs2cumu(c(0, 1, 1, 1)));
	print(hfs2cumu(c(0, 0, 1, 1)));
}
# data from
# http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2007-08_rel22/phased/
#
if (0) {
	d = readPhasedData('data/genotypes_chr6_CEU_r22_nr.b36_fwd');
	snps = SelectSNPsOrder('rs9269794', d, marginNeg = 1, marginPos = 1, maf = .05);
	htfs = haplotypeFrequencies(snps$hts);
}

if (1) {
	chr21 = readPhasedData('/Users/Lena/Desktop/internship/Data/genotypeschr21/genotypes_chr21_CEU_r22_nr.b36_fwd')
	s <- SelectSNPsOrder("rs8130901", chr21, marginNeg = 1, marginPos = 1, maf = .05)
	hs <- haplotypeFrequencies(s$hts)
	snp.list<-chr21$legend[100:199,1]  
}
