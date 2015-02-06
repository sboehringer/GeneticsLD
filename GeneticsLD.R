#
#	GeneticsLD.R
#Fri Feb  6 17:31:47 CET 2015

library('devtools');

if (T) {
	#system('rm GeneticsLD/src/*.o GeneticsLD/src/*.so');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsLD', threads = 6);
}

library('Rcpp');
#source('GeneticsLD/R/mcmc.R');

if (F) {
	library('GeneticsHaplotype');
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
