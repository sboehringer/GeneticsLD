#include <iostream>
#include <random>
#include <math.h>
#include "parameter.h"
#include "parametrizer.h"
using namespace std;

extern "C" {
#	include <sys/time.h>
}
vector<double>	runif(int N, default_random_engine &generator) {
	uniform_real_distribution<double>	distribution(0.0,1.0);
	vector<double>	r(N);
	for (int i = 0; i < N; i++) r[i] = distribution(generator);
	return r;
}

void	vector_print_idx(const Parameter<parameter_t> &p) {
	vector_print_idx_bin<parameter_t, 3>((vector<parameter_t>)p);
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;

#if 0
	vector< vector<uint> >	bits(4);
	for ( uint i = 0; i < 20; i++) {

		cout << "N: " << i << " # bits: " << countBitsSet(i) << "\n";
	}

	uint					Nloci = 3;
	vector< vector<uint> >	bitOrder(Nloci + 1);

	for ( uint j = 0; j < ((uint)1U << Nloci); j++) {
		bitOrder[countBitsSet(j)].push_back(j);
	}

	for ( uint j = 0; j < bitOrder.size(); j++) {
		cout << "Order " << j << ":";
		for (uint i = 0; i < bitOrder[j].size(); i ++) cout << " " << bitOrder[j][i];
		cout << "\n";
	}

#endif
#if 0
	for (uint i = 1; i < 6; i++) {
		PartitionCache	c(i);
		cout << "Partitions of: {0, .., " << i - 1 << "}\n";
		print_setPartitionBf(*c.partitions(i - 1));
		cout << "\n";
	}
#endif
#if 0
	ParameterHelper	ph(5);

	const parameter_t	htfs[] = {.1, .2, .3, .4};
	ph[2].printBitOrder();
	ParameterMultinomialL	p(2, htfs, ph);
	ph[2].printBitOrder();
	p.print();
	Parameter1sL	p1(p);
	p1.print();

	ParameterMultinomialL	p2(p1);
	p2.print();
#endif

#if 0
	ParameterHelper	ph(5);
	const parameter_t	htfs[] = {0.07677562, 0.20916280, 0.20237876, 0.09058649, 0.16948999, 0.02501366,
		0.16271583,0.06387685};
	ParameterMultinomialL	p(3, htfs, ph);
	ph[2].printBitOrder();
	p.print();
	Parameter1sL	p1(p);
	p1.print();

	ParameterMultinomialL	p2(p1);
	p2.print();

	cout << "p == p2: " << ( p2 == p ) << "\n";
	ParameterCumuL	p3(p2);
	p3.print();

	Parameter1sL	p4(p3);
	p4.print();

	ParameterMultinomialL	p5(p4);
	p5.print();

#endif

#define	benchmarkStart()\
	struct timeval t0, t1; \
	gettimeofday(&t0, 0)

#define	benchmarkStop(cycles)\
	gettimeofday(&t1, 0);\
	cout << "time taken per round (us): " \
		<< ((((double)t1.tv_sec - t0.tv_sec)*1e6 + ((double)t1.tv_usec - t0.tv_usec)) / (cycles)) << "\n"

#if 0
	ParameterHelper	ph(5);
	const parameter_t	htfs[] = {0.07677562, 0.20916280, 0.20237876, 0.09058649, 0.16948999, 0.02501366,
		0.16271583,0.06387685};

	ParameterMultinomialL	p(3, htfs, ph);

	p.print();
	struct timeval t0, t1;
	gettimeofday(&t0, 0);

	int	max = 1e7;
	for (int i = 0; i < max; i++) {
		Parameter1sL	p1(p);
		ParameterCumuL	p2(p1);
		Parameter1sL	p3(p2);
		ParameterMultinomialL	p4(p3);

		if (i == max - 1) p4.print();
	}
	gettimeofday(&t1, 0);

	cout << "time taken per round (us): "
		<< ((((double)t1.tv_sec - t0.tv_sec)*1e6 + ((double)t1.tv_usec - t0.tv_usec)) / max) << "\n";
#endif
#if 0
	for (int i = 1; i <= 4; i++) {
		PartitionCache	c(i);
		c.printBitOrder();
	}
#endif
#define	countOf(a)	(sizeof(a)/sizeof(a[0]))
#if 0
	Parametrizer	partrz;
	double			htfsRaw[] = {
		0.07677562, 0.20916280, 0.20237876, 0.09058649, 0.16948999, 0.02501366, 0.16271583,0.06387685
	};
	vector<double>	htfs(htfsRaw, htfsRaw + countOf(htfsRaw));
	vector_print<double>(htfs);

	vector<double>	&cumu = *partrz.multinomial2cumu(htfs);
	vector_print<double>(cumu);

	benchmarkStart();
	int	max = 1e7;
	for (int i = 0; i < max; i++) {
		vector<double>	&m = *partrz.cumu2multinomial(cumu);
		delete &m;
	}
	benchmarkStop(max);
#endif
#if 1
	typedef	vector<double>	vec_t;
	Parametrizer	partrz;
	// hf = runif(16); cat(paste(hf/sum(hf), collapse = ", "))
#	if 0
	double			htfsRaw[] = {
		0.0896666374000689, 0.108660904671737, 0.0816669844186323, 0.0274091349780836,
		0.0215752723966624, 0.0369767780105436, 0.077898739329432, 0.157790798703374,
		0.0321247895230421, 0.0682787993403603, 0.0925898640044041, 0.0701397006701806,
		0.00847218221107304, 0.0210534141261667, 0.00317682032732407, 0.102519179888916
	};
#	endif
#	if 0
	double			htfsRaw[] = {
		0.0415506891576117, 0.0410196622135817, 0.00161803517254192, 0.0235395605388719,
		0.0408408757794127, 0.0197700130856936, 0.0404997400938796, 0.0377055316487111,
		0.0422272578946939, 0.0487954821149679, 0.0223406819489681, 0.0271495083704544,
		0.0469946533374358, 0.0445242188712996, 0.0472065747390216, 0.0329942163867904,
		0.0415551828213593, 0.00896621033366813, 0.00535626886470554, 0.0499713483760888,
		0.000625652432077811, 0.00885615813621108, 0.0472792711977232, 0.0186285483397111,
		0.0188551258975844, 0.0356145255177474, 0.0244786982674912, 0.0505026483091139,
		0.0509329284055216, 0.0177078985496228, 0.0129840559892064, 0.0489087772082314
	};
#	endif
#	if 0
	double			htfsRaw[] = {
		0.0087848328141055, 0.0367626341660437, 0.0244560291917515, 0.0478164938335441,
		0.0094027197342296, 0.0632396485189338, 0.0347061696655587, 0.065933493912382,
		0.0891749449989024, 0.0778220832998305, 0.115625308083933, 0.101974726359158,
		0.0678744407732572, 0.0394853870720869, 0.111595945901245, 0.105345141675038
	};
#	endif

#	if 0
	cout.precision(4);
	vector<double>	htfs(htfsRaw, htfsRaw + countOf(htfsRaw));

	vector<double>	&p1s0 = *partrz.multinomial2p1s(htfs);

	vector<double>	&cumu = *partrz.multinomial2cumu(htfs);
	vector<double>	&cumuStd = *partrz.multinomial2cumuStd(htfs);
	cout << "cumu "; vector_print_idx_bin<double>(cumu);
	cout << "cumuStd "; vector_print_idx_bin<double>(cumuStd);

	vector<double>	&p1s = *partrz.cumu2p1s(cumu);
	vector<double>	&p1s1 = *partrz.cumuStd2p1s(cumuStd);
	cout << "p1s  "; vector_print_idx_bin<double>(p1s0);
	cout << "p1s  "; vector_print_idx_bin<double>(p1s);
	cout << "p1s  "; vector_print_idx_bin<double>(p1s1);

	vector<double>	&m = *partrz.cumu2multinomial(cumu);
	vector<double>	&m1 = *partrz.cumuStd2multinomial(cumuStd);
	cout << "htfs "; vector_print_idx_bin<double>(htfs);
	cout << "htfs "; vector_print_idx_bin<double>(m);
	cout << "htfs "; vector_print_idx_bin<double>(m1);
#	endif


#endif
#if 0
	typedef	vector<double>	vec_t;
	Parametrizer	partrz;
	// hf = runif(16); cat(paste(hf/sum(hf), collapse = ", "))
	// hf = runif(32); cat(paste(hf/sum(hf), collapse = ", "))
	double			htfsRaw[] = {
		0.0896666374000689, 0.108660904671737, 0.0816669844186323, 0.0274091349780836,
		0.0215752723966624, 0.0369767780105436, 0.077898739329432, 0.157790798703374,
		0.0321247895230421, 0.0682787993403603, 0.0925898640044041, 0.0701397006701806,
		0.00847218221107304, 0.0210534141261667, 0.00317682032732407, 0.102519179888916
	};
	vector<double>	htfs(htfsRaw, htfsRaw + countOf(htfsRaw));
	vector_print<double>(htfs);

	vector<double>	&cumu = *partrz.multinomial2p1s(htfs);
	vector_print<double>(cumu);

	vector<double>	&m = *partrz.p1s2multinomial(cumu);
	vector_print<double>(m);
	

#endif
#	if 0
	vector<double>	r = runif(16);
	cout << "Runif: ";
	vector_print(r);
#	endif

#	if 0
#	define	N	8
	//Parametrizer	partrz;
	default_random_engine	generator;
	cout.precision(2);
	for (int j = 0; j < (int)1e1; j++) {
		vector<double>	r = runif(N, generator);
		r[0] = 0;
		vector<double>	&m1 = *partrz.cumuStd2multinomial(r);
		vector<double>	&m2 = *partrz.cumuStd2p1s(r);
		vector<double>	&m3 = *partrz.p1s2multinomial(m2);
		for (int i = 0; i < m1.size(); i++) {
			if (m1[i] < 0 || m1[i] > 1) {
				cout << j << "\n";
				cout << "cumuStd: "; vector_print_idx_bin<double>(r);
				cout << "1s:      "; vector_print_idx_bin<double>(m2);
				cout << "multin:  "; vector_print_idx_bin<double>(m1);
				break;
			}
		}
		delete &m1;
	}

#	endif

#	if 0
	ParameterHelper			h(5);
	double					htfsRaw[] = { 0.125, 0, .375, .5 };
	vector<double>			htfs(htfsRaw, htfsRaw + countOf(htfsRaw));
	ParameterMultinomialL	pm(htfs, h);
	Parameter1sL			p1s(pm);
	Parameter1sMinMaxL		p1sMnMx = p1s.minMax();
	Parameter1sStdL			p1std(p1s);
	
	pm.print();
	p1s.print();
	p1std.print();
	p1sMnMx.min.printP("min");
	p1sMnMx.max.printP("max");
#	endif
#	if 0
	ParameterHelper			h(5);
	double					htfsRaw[] = { 0.125, 0, .375, .5, 0.125, 0, .375, .5 };
	vector<double>			htfs(htfsRaw, htfsRaw + countOf(htfsRaw));
	ParameterMultinomialL	pm(htfs, h);
	Parameter1sL			p1s(pm);
	Parameter1sMinMaxL		p1sMnMx = p1s.minMax();
	Parameter1sStdL			p1std(p1s);
	
	pm.print();
	p1s.print();
	p1std.print();
	p1sMnMx.min.printP("min");
	p1sMnMx.max.printP("max");
#	endif

#	if 0
#	define	N	8
	default_random_engine	generator;
	const ParameterHelper			h(5);
	for (int j = 0; j < (int)1e5; j++) {
		// Parameter1sStdL
		//double					rRaw[] = { 1, .46, .22, .68, .93, .52, .035, .53 };
		//double					rRaw[] = { 1, .067, .69, .93, .53, .65, .7, .76 };
		//double					rRaw[] = { 1, .46, .22, .68, .93, .52, .035, .53 };
		//double					rRaw[] = { 1, .56, .49, .96, .2, .63, .65, .8 };
		//double					rRaw[] = { 1, .28, .60, .063, .95, .28, .43, .199 };
		//double					rRaw[] = { 1, .46, .45, .93, .22, .91, .86, .51 };
		//vector<double>			r(rRaw, rRaw + countOf(rRaw));
		const vector<double>	r(runif(N, generator));
		Parameter1sStdL			p1std(r, h);
		Parameter1sL			p1s(p1std);
		ParameterMultinomialL	pm(p1s);
		Parameter1sL			p1sB(pm);
		Parameter1sStdL			p1stdB(p1sB);
		Parameter1sMinMaxL		p1sMnMx = p1s.minMax();

		if (!(j % (int)1e4)) cout << "Iteration: " << j << "\n";
		for (int i = 0; i < pm.size(); i++) {
			//if (abs(p1stdB[i] - p1std[i]) > 1e-4 || pm[i] < 0 || pm[i] > 1) {
			if (abs(p1stdB[i] - p1std[i]) > 1e-4 || pm[i] > 1) {
				cout.precision(6);
				cout << "1sStd    "; vector_print_idx(p1std);
				cout << "1s:      "; vector_print_idx_bin<parameter_t, 3>(p1s);
				cout << "multin:  "; vector_print_idx_bin<parameter_t, 3>(pm);
				cout << "1sB:     "; vector_print_idx_bin<parameter_t, 3>(p1sB);
				cout << "1sStdB:  "; vector_print_idx_bin<parameter_t, 3>(p1stdB);
				cout.precision(2);
				p1sMnMx.min.printP("min");
				p1sMnMx.max.printP("max");
				break;
			}
		}
	}
#	endif

#	if 0
#	define	N	8
	default_random_engine	generator;
	const ParameterHelper	h(5);
	for (int j = 0; j < (int)1e5; j++) {
		// Parameter1sStdL
		//double					rRaw[] = { 1, .46, .22, .68, .93, .52, .035, .53 };
		//double					rRaw[] = { 1, .067, .69, .93, .53, .65, .7, .76 };
		//double					rRaw[] = { 1, .46, .22, .68, .93, .52, .035, .53 };
		//double					rRaw[] = { 1, .56, .49, .96, .2, .63, .65, .8 };
		//double					rRaw[] = { 1, .28, .60, .063, .95, .28, .43, .199 };
		//double					rRaw[] = { 1, .46, .45, .93, .22, .91, .86, .51 };
		//vector<double>			r(rRaw, rRaw + countOf(rRaw));
		const vector<double>	r(runif(N, generator));
		ParameterMultinomialL	pm(r, h);
		Parameter1sL			p1s(pm);
		Parameter1sStdL			p1std(p1s);
		Parameter1sL			p1sB(p1std);
		ParameterMultinomialL	pmB(p1sB);

		if (!(j % (int)1e4)) cout << "Iteration: " << j << "\n";
		for (int i = 0; i < pm.size(); i++) {
			//if (abs(p1stdB[i] - p1std[i]) > 1e-4 || pm[i] < 0 || pm[i] > 1) {
			if (abs(pm[i] - pmB[i]) > 1e-4 || pm[i] > 1) {
				cout.precision(3);
				cout << "multin:  "; vector_print_idx_bin<parameter_t, 3>(pm);
				cout << "1s:      "; vector_print_idx_bin<parameter_t, 3>(p1s);
				cout << "1sStd    "; vector_print_idx_bin<parameter_t, 3>(p1std);
				cout << "1sB:     "; vector_print_idx_bin<parameter_t, 3>(p1sB);
				cout << "multinB: "; vector_print_idx_bin<parameter_t, 3>(pmB);
				cout.precision(2);
				break;
			}
		}
	}

#	endif

#	if 0
#	define	N	8
	default_random_engine	generator;
	const ParameterHelper	h(5);
	for (int j = 0; j < (int)1e0; j++) {
		// Std = -1
		//double					rRaw[] = { 0, 0, .05, 0.1, 0, 0, .3, .5 };
		// |std| > 3
		//double					rRaw[] = { 0, 0, .05, 0, .1, 0, .3, .5 };
		//double					rRaw[] = { .1, 0, .05, 0, .1, .2, .3, .5 };
		double					rRaw[] = { .476, .0952, 0, 0, .0952, 0, .286, .0476 };	// <!> 111: 1.5
		vector<double>			r(rRaw, rRaw + countOf(rRaw));
		//const vector<double>	r(runif(N, generator));
		ParameterMultinomialL	pm(r, h);
		Parameter1sL			p1s(pm);
		Parameter1sStdL			p1std(p1s);
		Parameter1sL			p1sB(p1std);
		ParameterMultinomialL	pmB(p1sB);

		if (!(j % (int)1e4)) cout << "Iteration: " << j << "\n";
		for (int i = 0; i < pm.size(); i++) {
			//if (abs(p1stdB[i] - p1std[i]) > 1e-4 || pm[i] < 0 || pm[i] > 1) {
			if (1 || abs(pm[i] - pmB[i]) > 1e-4 || pm[i] > 1) {
				cout.precision(3);
				cout << "multin:  "; vector_print_idx_bin<parameter_t, 3>(pm);
				cout << "1s:      "; vector_print_idx_bin<parameter_t, 3>(p1s);
				cout << "1sStd    "; vector_print_idx_bin<parameter_t, 3>(p1std);
				cout << "1sB:     "; vector_print_idx_bin<parameter_t, 3>(p1sB);
				cout << "multinB: "; vector_print_idx_bin<parameter_t, 3>(pmB);
				cout.precision(2);
				break;
			}
		}
	}

#	endif

#	if 1
#	define	Nloci	3
#	define	N		(1 << Nloci)
	default_random_engine	generator;
	const ParameterHelper	h(5);
	// Std = -1
	//double					rRaw[] = { 0, 0, .05, 0.1, 0, 0, .3, .5 };
	// |std| > 3
	//double					rRaw[] = { 0, 0, .05, 0, .1, 0, .3, .5 };
	//double					rRaw[] = { .1, 0, .05, 0, .1, .2, .3, .5 };
	//double					rRaw[] = { .1, 0, .05, .5, .1, .2, .3, 0 };
	//double					rRaw[] = { .1, .3, .05, .5, .1, .2, 0, 0 };
	double					rRaw[] = { .1, 0, .3, .05, .5, .1, 0, 0 };	// <!> 111: 1.5
	//double					rRaw[] = { .1, 0, 0, .05, .5, .3, 0, 0 };
	//double					rRaw[] = { 0, .3, .05, .5, .1, .2, .2, .1 };
	vector<double>			r(rRaw, rRaw + countOf(rRaw));

	ParameterMultinomialL	pm(r, h);
	Parameter1sL			p1s(pm);

	for (int j = 0; j < (int)N; j++) {
		
		Parameter1sL			p1sf(p1s, (haplotype_t)j);
		Parameter1sStdL			p1std(p1sf);
		Parameter1sL			p1sB(p1std);
		ParameterMultinomialL	pmB(p1sB);

		cout.precision(3);
		cout << "Flipmask:" << j << endl;
		cout << "multin:  "; vector_print_idx_bin<parameter_t, 3>(pmB);
		cout << "1s:      "; vector_print_idx_bin<parameter_t, 3>(p1sf);
		cout << "1sStd    "; vector_print_idx_bin<parameter_t, 3>(p1std);
		cout.precision(2);
	}

#	endif

	cout.flush();
	return 0;
}