#
#	GeneticsLD.R
#Fri Feb  6 17:31:47 CET 2015

library('devtools');

if (T) {
	system('rm GeneticsLD/src/*.o GeneticsLD/src/*.so');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsLD', threads = 6);
}

library('Rcpp');
#source('GeneticsLD/R/mcmc.R');

if (F) {
	library('GeneticsHaplotype');
	#M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
}

