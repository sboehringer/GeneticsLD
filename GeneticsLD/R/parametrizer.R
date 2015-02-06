#
#	parametrizer.R
#Mon Sep 30 11:39:27 CEST 2013
#
#	build a module using inline
#

system("cd ~/src/Rprivate ; ./exportR.sh"); source('RgenericAll.R');

if (0) {
	library('Rcpp');
	library('inp1sline');
	dir = 'build';
	lib = 'haplotypes_rcpp';
	.libPaths(c(.libPaths(), dir));
	#if (!length(fetchRegexpr(lib, Sys.getenv('PKG_LIBS'))))
	Sys.setenv(`PKG_LIBS` = sprintf('%s -L"%s" -l%s', Sys.getenv('PKG_LIBS'), splitPath(dir)$absolute, lib));
	Sys.setenv(`PKG_CXXFLAGS` = sprintf('%s %s', Sys.getenv('PKG_LIBS'), stdOutFromCall(Rcpp:::CxxFlags())));


	dyn.load('build/libhaplotypes_rcpp.so');
	moduleRegex = '(?s:(?<=// -- begin inline Rcpp\n)(.*?)(?=// -- end inline Rcpp))';
	inc = join(sapply(c('parameterhelper.h', 'parametrizer.h', 'parametrizerrcpp.h', 'parametrizerrcpp_module.h'), function(f) fetchRegexpr(moduleRegex, readFile(f))), "\n");

	rcpp = cxxfunction( signature(), '' , includes = inc, plugin = 'Rcpp', verbose = T, settings = list(LinkingTo = 'haplotypes_rcpp') );
	mod = Module( "module_parametrizer", getDynLib(rcpp));

	par = new(mod$parametrizer);
	fqs = vector.std(1:8);
	f1 = par$multinomial2cumu(fqs);
	f0 = par$cumu2multinomial(f1);
}

if (1) {
	mod = createModule('module_parametrizer', 'build/libhaplotypes_rcpp.so',
		headers = c('parameterhelper.h', 'parametrizer.h', 'parametrizerrcpp.h', 'parametrizerrcpp_module.h'),
		output = 'mod'
	);
	System('tar czf mod.tgz mod');
}
if (0) {
	mod = createModule('module_parametrizer', 'libhaplotypes_rcpp.so',
		headers = c('parameterhelper.h', 'parametrizer.h', 'parametrizerrcpp.h', 'parametrizerrcpp_module.h'),
		output = 'mod'
	);
	System('tar czf mod.tgz mod');
}

if (0) {
	mod = activateModule('mod');

	par = new(mod$parametrizer);
	for (j in 1:7) {
		print(sprintf("#loci: %d", j));
		t = system.time({
		for (i in 1:1e5) {
			fqs = vector.std(runif(2^j));
			f1 = par$multinomial2cumu(fqs);
			f0 = par$cumu2multinomial(f1);
			if (any(abs(f0 - fqs) > 1e-10)) stop('mapping not bijective');
			#if (any(abs(f0 - fqs) > 1e-13)) cat('mapping inaccurate at 1e-13 level\n');
		}
		});
		print(t);
		print('so far so good');
	}
}

if (0) {
	parNames = function(i)sapply(0:(2^i - 1), function(j)paste(rev(ord2bin(j, i)), collapse = ''));
}
