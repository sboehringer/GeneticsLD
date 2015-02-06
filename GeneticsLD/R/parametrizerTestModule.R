#
#	parametrizerTestModule.R
#Tue Mar 25 18:34:37 2014
#
system("cd ~/src/Rprivate ; ./exportR.sh"); source('RgenericAll.R');

if (T) {
	mod = activateModule('mod');
	par = new(mod$parametrizer);
}

if (0) {

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

if (1) {
	j = 4;
	fqs = vector.std(runif(2^j));
	f1 = par$multinomial2p1s(fqs);
	f2 = par$multinomial2p1sStd(fqs);
	f1mm = cbind(matrix(par$p1sMinMax(f1), ncol = 2), f1, f2);

	f3 = par$multinomial2cumu(fqs);
	f4 = par$multinomial2cumuStd(fqs);
	f3mm = cbind(matrix(par$cumuMinMax(f3), ncol = 2), f3, f4);
	
	print(f1mm);
	print(f3mm);
}
