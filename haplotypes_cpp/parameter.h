/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2013  Stefan Boehringer <email>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef PARAMETER_H
#define PARAMETER_H

#include "parameterhelper.h"
#include <iomanip>
#include <algorithm>
#include <memory>
#include <vector>

using namespace std;

#define	Nhts(loci)	((1 << (loci)))
inline int	log2ceil(unsigned int i) {
	if (!i) return -1;
	int	lg2 = sizeof(unsigned int)*8 - __builtin_clz (i) - 1;
	if (!lg2) return 0;
	lg2 += !!(i & ((1 << lg2) - 1));	// check for remainder
	return lg2;
}

template<class T> void vector_print(const vector<T> &v) {
	for (int i = 0; i < v.size(); i++) {
		cout << (i? ", ": "") << v[i];
	}
	cout << "\n";
}

template<class T> void vector_print_idx_bin(const vector<T> &v) {
	for (int i = 0; i < v.size(); i++) {
		cout << (i? ", ": "") << bitset<5>(i) << ':' << v[i];
	}
	cout << "\n";
}
template<class T, int C> void vector_print_idx_bin(const vector<T> &v) {
	for (int i = 0; i < v.size(); i++) {
		cout << (i? ", ": "") << bitset<C>(i) << ':' << v[i];
	}
	cout << "\n";
}

template<class TF, class TT> TT	identity(TF from) { return (TT)from; }
// convert vector type from class TF to class TT
template<class TF, class TT> vector<TT>	*vectorConvert(const vector<TF> &vf) {
	vector<TT>	&vt = * new vector<TT>();
	vt.resize( vf.size() );
	transform( vf.begin(), vf.end(), vt.begin(), identity<TF, TT> );
	return &vt;
}

typedef long double		parameter_t;
typedef unsigned int	haplotype_t;
template <class T> class Parameter1s;
template <class T> class Parameter1sStd;
template <class T> class ParameterCumu;
template <class T> class ParameterCumuStd;
/*
 * multinomial parameter vector
 * overparametrized to include all probabilities
 * vector is indexed by haplotype number, which is lexicographic order of vector of alleles
 */
template <class T> class Parameter : public vector<T>
{
protected:
	int 					Nloci;
	const PartitionCache	&c;

public:
    Parameter(const vector<double> &values, const ParameterHelper &ph)
    : vector<T>(*auto_ptr< vector<T> >(vectorConvert<double, T>(values)))
	, Nloci(log2ceil(values.size()))
	, c(ph[Nloci])
	{};
// 	: Nloci(log2ceil(values.size()))
// 	, c(ph[Nloci])
// 	{
// 		std::auto_ptr< vector<T> >	v(vectorConvert<double, T>(values));
// 		(vector<T>)(*this) = *v;
// 	};
    Parameter(int nloci_, ParameterHelper &ph)
	: vector<T>(Nhts(nloci_), (T)0.0), Nloci(nloci_), c(ph[nloci_])
	{};
	Parameter(int nloci_, const T *values, ParameterHelper &ph)
	: vector<T>(values, values + Nhts(nloci_)), Nloci(nloci_), c(ph[nloci_])
	{};
    Parameter(const Parameter<T>& other)
	: vector<T>((vector<T>)other), Nloci(other.Nloci), c(other.c)
	{};
    Parameter(const PartitionCache &c_)
    : vector<T>(Nhts(c_.Nelements), (T)0.0), Nloci(c_.Nelements), c(c_)
	{};
    virtual ~Parameter();
    virtual Parameter<T>& operator=(const Parameter<T>& other);
    virtual bool operator==(const Parameter<T>& other) const;
	void	print(int precision = 3, const char *prefix = "") const;
	void	printP(const char *prefix = "") const { this->print(3, prefix); };
	virtual	const char	*className(void) const { return "Generic"; }
	inline int	nloci(void) { return Nloci; }
	inline const PartitionCache	partitionCache(void) { return c; }
	inline const vector<T> *as_vector(void) { return this; }
};
typedef Parameter<parameter_t>	ParameterL;

/*
 * implementation
 */

template <class T> Parameter<T>::~Parameter()
{

}

template <class T> Parameter<T>&  Parameter<T>::operator=(const Parameter<T>& other)
{
    return *this;
}

template <class T> bool Parameter<T>::operator==(const Parameter<T>& other) const
{
	return ((vector<T>) (*this) == (vector<T>) other);
}

template <class T> void Parameter<T>::print(int precision, const char * prefix) const {
	using namespace	std;
	cout << setprecision(5);

	cout << "[" << className() << "]" << prefix << ":";
	for (int i = 0; i < this->size(); i++) {
		cout << " " << bitset<4>(i) << ':' << (*this)[i];
	}
	cout << "\n";
}

/*
 * haplotype frequencies
 */
template <class T> class ParameterMultinomial : public Parameter<T>
{
public:
    ParameterMultinomial(const vector<double> &values, const ParameterHelper &ph) : Parameter<T>(values, ph) {
		T	sum = 0;
		for (int i = 0; i < this->size(); i++) sum += (*this)[i];
		for (int i = 0; i < this->size(); i++) (*this)[i] /= sum;
	};
	ParameterMultinomial(int nloci_, const T *values, const ParameterHelper &ph) : Parameter<T>(nloci_, values, ph) {};
    ParameterMultinomial(const Parameter1s<T> &p) : Parameter<T>((Parameter<T>)p)  {
		from1s(p);
	};
	void from1s(const Parameter1s<T>& other);
	virtual	const char	*className(void) const { return "Multinomial"; }
};

template <class T> void ParameterMultinomial<T>::from1s(const Parameter1s<T>& other) {
	// iterate from long haplotypes to small ones
	//cout << "\nmultinomial from1s\n";
	parameter_t	sum = 0;
	for (int i = this->Nloci; i > 0; i--) {
		for (int j = 0; j < this->c.bitOrder[i].size(); j++) {
			haplotype_t	h = this->c.bitOrder[i][j];
			T	f = other[h];	// raw haplotype frequency
			bitidx_t	zeros = bitCountUnset(h, this->Nloci);	// count number or zeros
			//cout << "FrequencyS: " << f << "\n";
			// iterate over all values with less zeros for which ones are fixed
			for (int k = 1; k < (1 << zeros); k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(h, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, this->Nloci);	// count number or zeros
				//cout << "Haplotypes: " << h << ", " << h1 << " k: " << k << " zeros: " << zeros << ", " << zeros1 << '\n';
				// even or odd number of less zeros, assume h1-parameter to already be compouted
				if ((zeros1 - zeros) & 1)
					f -= other[h1];
				else
					f += other[h1];
			}
			(*this)[h] = f;
			sum += f;
		}
	}
	(*this)[0] = 1 - sum;
}
typedef ParameterMultinomial<parameter_t>	ParameterMultinomialL;

/*
 * probability vector in terms of all marginal only-1-allele haplotypes
 * vector is 
 */
template<class T> struct Parameter1sMinMax {
	const Parameter1s<T>	min, max;
};
typedef Parameter1sMinMax<parameter_t>	Parameter1sMinMaxL;

template <class T> class Parameter1s : public Parameter<T>
{
public:
    Parameter1s(const vector<double> &values, const ParameterHelper &ph) : Parameter<T>(values, ph) {};
	Parameter1s(int nloci_, const ParameterHelper &ph) : Parameter<T>(nloci_, ph) {};
    Parameter1s(const ParameterMultinomial<T> &p) : Parameter<T>((Parameter<T>)p)  {
		fromMultinomial(p);
	};
    Parameter1s(const Parameter1sStd<T> &p) : Parameter<T>((Parameter<T>)p)  {
		from1sStd(p);
	};
    Parameter1s(const ParameterCumu<T> &p) : Parameter<T>((Parameter<T>)p)  {
		fromCumu(p);
	};
    Parameter1s(const ParameterCumuStd<T> &p) : Parameter<T>((Parameter<T>)p)  {
		fromCumuStd(p);
	};
    Parameter1s(const PartitionCache &c_) : Parameter<T>(c_)
	{};
    Parameter1s(const Parameter1s<T> &other, haplotype_t flipMask) : Parameter<T>(other.c) {
		flipped(other, flipMask);
	};

	void	fromMultinomial(const ParameterMultinomial<T>& other);
	void	from1sStd(const Parameter1sStd<T>& other);
	void	fromCumu(const ParameterCumu<T>& other);
	void	fromCumuStd(const ParameterCumuStd<T>& other);
	void	flipped(const Parameter1s<T>& other, haplotype_t flipMask);
	virtual	const char	*className(void) const { return "1s"; }
	Parameter1sMinMax<T>	minMax(void) const;
};
typedef Parameter1s<parameter_t>	Parameter1sL;

template <class T> void Parameter1s<T>::fromMultinomial(const ParameterMultinomial<T>& other) {
	// iterate from long haplotypes to small ones in terms of contained 1-alleles
	//cout << "\n1s fromMultinomial\n";
	(*this)[0] = 1;	// "all dot" haplotype, marginal of everything
					// formally, a redundant parameter
	for (int i = this->Nloci; i > 0; i--) {
		for (int j = 0; j < this->c.bitOrder[i].size(); j++) {
			haplotype_t	h = this->c.bitOrder[i][j];
			T			f = other[h];	// raw haplotype frequency
			bitidx_t	zeros = bitCountUnset(h, this->Nloci);	// count number or zeros
			//cout << "Haplotype: " << h << " zeros: " << zeros << '\n';
			// iterate over all values with less zeros for which ones are fixed
			for (int k = 1; k < (1 << zeros); k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(h, k, zeros);
				//cout << "Haplotype1: " << h1 << " k: " << k << '\n';
				bitidx_t	zeros1 = bitCountUnset(h1, this->Nloci);	// count number or zeros
				// even or odd number of less zeros, assume h1-parameter to already be compouted
				if ((zeros1 - zeros) % 2)
					f += (*this)[h1];
				else
					f -= (*this)[h1];
			}
			(*this)[h] = f;
		}
	}
}

template <class T> void Parameter1s<T>::flipped(const Parameter1s<T>& other, haplotype_t flipMask) {
	// see fromMulitnomial for comments
	(*this)[0] = 1;
	// iterate haplotypes, order does not matter
	for (int ones = this->Nloci; ones >= 1; ones--) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		haplotype_t	hflip = h & flipMask;	// which alleles to flip in this haplotype?
		bitidx_t	onesF = bitCountSet(hflip, this->Nloci);	// number of loci to flip
		T			f = 0;	// unflipped haplotype frequency
		//cout << "Haplotype: " << h << " loci to flip: " << onesF << " flipMask: " << flipMask << '\n';
		// iterate over all values with less ones for which outer ones are fixed
		for (int k = 0; k < (1 << onesF); k++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(hflip, k, onesF);
			haplotype_t	h2 = (h & ~flipMask) | h1;
			bitidx_t	ones1 = bitCountSet(~k, onesF);	// count number or zeros
			// even or odd number of less zeros, assume h1-parameter to already be compouted
			T	inc = ((onesF - ones1) & 1)? -(other)[h2]: (other)[h2];
			f += inc;
			//cout << "Haplotype1: " << h1 << " h2:" << h2 << " k: " << k << " inc:" << inc << " f:" << f << '\n';
		}
		(*this)[h] = f;
	}
}

template <class T> struct Parameter1sMinMax<T> Parameter1s<T>::minMax(void) const {
	// <p> return vector
	Parameter1s<T>	maxFs(this->c);
	Parameter1s<T>	minFs(this->c);

	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			minFs[h] = 0;
			maxFs[h] = 1;
			continue;
		}
		// <p> minimum/maximum computation
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
// 		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, (*this)[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -(*this)[h2]: (*this)[h2];
			}
			//min1s = max(min1s, (zeros & 1)? r: -r);
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
		}
		minFs[h] = min1s;
		maxFs[h] = max1s;
	}
	minFs[0] = maxFs[0] = 1;
	return (Parameter1sMinMax<T>){ minFs, maxFs };
}

template <class T> class Parameter1sStd : public Parameter<T>
{
public:
    Parameter1sStd(const vector<double> &values, const ParameterHelper &ph) : Parameter<T>(values, ph) {
		(*this)[0] = 1;
	};
	Parameter1sStd(int nloci_, const ParameterHelper &ph) : Parameter<T>(nloci_, ph) {};
    Parameter1sStd(const Parameter1s<T> &p) : Parameter<T>((Parameter<T>)p)  {
		from1s(p);
	};
    Parameter1sStd(const PartitionCache &c_) : Parameter<T>(c_)
	{};

	void	from1s(const Parameter1s<T>& other);
	virtual	const char	*className(void) const { return "1sStd"; }
};
typedef Parameter1sStd<parameter_t>	Parameter1sStdL;

/*
 * Standardization of 1-allele parameterization
 * Needs to combine functions flipped, both loops from minMax above into one loop
 * This is to maintain symmetry in the computation to be able to reverse the computation
 */

template <class T> void Parameter1sStd<T>::from1s(const Parameter1s<T>& other) {
	// iterate haplotypes, from low to high order
	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			(*this)[h] = other[h];
			continue;
		}
		// <p> minimum/maximum computation
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
// 		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, other[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -other[h2]: other[h2];
			}
			//min1s = max(min1s, (zeros & 1)? r: -r);
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
		}
		(*this)[h] = !(max1s - min1s)? 1: ((other[h] - min1s) / (max1s - min1s));
	}
}


template <class T> void Parameter1s<T>::from1sStd(const Parameter1sStd<T>& other) {
	// iterate haplotypes, from low to high order
	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			(*this)[h] = other[h];
			continue;
		}
		// <p> minimum/maximum computation
		// maximum depends on minimum of lower order; start at pairs
		//	we iterate the partitions, it is only necessary to look at haplotypes of length -1
		//	as such data structures are not available, iterate partitions which is a superset of needed haplotypes
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, (*this)[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -(*this)[h2]: (*this)[h2];
				//cout << "   mult-resid-l[" << h2 << "]: " << (*this)[h2] << " r:" << r << endl;
			}
			//min1s = max(min1s, (zeros & 1)? r: -r);
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
			//cout << "  mult-resid[" << hm << "]: " << r << " (min1s, max1s) = (" << min1s << ", " << max1s << ")" << endl;
		}
		(*this)[h] = other[h] * (max1s - min1s) + min1s;
		//cout << "halpotype:" << h << "  min1s:" << min1s << " max1s:" << max1s << endl;
		//cout << "--" << endl;
	}
}

/*
 * vector of cumulants of haplotype frequencies
 */

template<class T> struct ParameterCumuMinMax {
	const ParameterCumu<T>	min, max;
};
typedef ParameterCumuMinMax<parameter_t>	ParameterCumuMinMaxL;

template <class T> class ParameterCumu : public Parameter<T>
{
public:
    ParameterCumu(vector<double> &values, ParameterHelper &ph) : Parameter<T>(values, ph) {};
	ParameterCumu(int nloci_, ParameterHelper &ph) : Parameter<T>(nloci_, ph) {};
    ParameterCumu(const Parameter1s<T> &p) : Parameter<T>((Parameter<T>)p)  {
		from1s(p);
	};


	void	from1s(const Parameter1s<T>& other);
	virtual	const char	*className(void) const { return "Cumu"; }
	struct ParameterCumuMinMax<T>	minMax(void) const;
};
typedef ParameterCumu<parameter_t>	ParameterCumuL;

#include <cmath>

template <class T> void ParameterCumu<T>::from1s(const Parameter1s<T>& other) {
	// iterate from cumulants of low order to cumulants of high order
	//cout << "\ncumu from1s\n";
	(*this)[0] = 0;	// undefined, there is no haplotype with 0 loci
	for (int i = 1; i <= this->Nloci; i++) {
		for (int j = 0; j < this->c.bitOrder[i].size(); j++) {
			// binary represenation: 1 represents retained locus, 0: marginalization
			haplotype_t	h = this->c.bitOrder[i][j];
			T			cum = 0;
			bitidx_t	ones = bitCountSet(h, this->Nloci);				// count number or ones
			const SetPartitionBf	&parts = *(this->c.partitions(ones - 1));	// all paritions with 'ones' time components
			// iterate over all values with more zeros for which zeros are fixed at positions from h
			// skip first element to make implementations symmetric
			for (int k = 1; k < parts.size(); k++) {
				int	partSize = parts[k].size();
				T	t = (partSize % 2)? tgamma(partSize): -tgamma(partSize);
				for (int l = 0; l < partSize; l++) {
					haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, parts[k][l], ones);
					t *= other[h1];
				}
				cum += t;
			}
			(*this)[h] = other[h] + cum;
		}
	}
}

template <class T> void Parameter1s<T>::fromCumu(const ParameterCumu<T>& other) {
	// iterate from cumulants of low order to cumulants of high order
	//cout << "\n1s from cumu\n";
	(*this)[0] = 0;	// undefined, there is no haplotype with 0 loci
	for (int i = 1; i <= this->Nloci; i++) {
		for (int j = 0; j < this->c.bitOrder[i].size(); j++) {
			// binary represenation: 1 represents retained locus, 0: marginalization
			haplotype_t	h = this->c.bitOrder[i][j];
			T			cum = 0;
			bitidx_t	ones = bitCountSet(h, this->Nloci);				// count number or ones
			const SetPartitionBf	&parts = *(this->c.partitions(ones - 1));	// all paritions with 'ones' time components
			// iterate over all values with less zeros for which ones are fixed, skip first element
			for (int k = 1; k < parts.size(); k++) {
				int	partSize = parts[k].size();
				T	t = (partSize % 2)? tgamma(partSize): -tgamma(partSize);
				for (int l = 0; l < partSize; l++) {
					haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, parts[k][l], ones);
					t *= (*this)[h1];
				}
				cum += t;
			}
			(*this)[h] = other[h] - cum;
		}
	}
}
template <class T> struct ParameterCumuMinMax<T> ParameterCumu<T>::minMax(void) const {
	// <p> return vector
	ParameterCumu<T>	maxFs(this->c);
	ParameterCumu<T>	minFs(this->c);
	Parameter1s<T>		other(this->c);
	other[0] = 1;

	// iterate haplotypes, from low to high order
	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			other[h] = (*this)[h];
			minFs[h] = 0;
			maxFs[h] = 1;
			continue;
		}
		// <p> minimum/maximum computation
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
// 		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, other[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -other[h2]: other[h2];
			}
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
		}

		// cumulant part
		T			cum = 0;
		// all paritions with 'ones' time components
		const SetPartitionBf	&parts = *(this->c.partitions(ones - 1));
		// iterate over all values with more zeros for which zeros are fixed at positions from h
		// skip first element to make implementations symmetric
		for (int k = 1; k < parts.size(); k++) {
			int	partSize = parts[k].size();
			T	t = (partSize % 2)? tgamma(partSize): -tgamma(partSize);
			for (int l = 0; l < partSize; l++) {
				haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, parts[k][l], ones);
				t *= other[h1];
			}
			cum += t;
		}
		minFs[h] = min1s + cum;
		maxFs[h] = max1s + cum;
		other[h] = (*this)[h] - cum;
	}
	minFs[0] = maxFs[0] = 1;
	return (ParameterCumuMinMax<T>){ minFs, maxFs };
}

template <class T> class ParameterCumuStd : public Parameter<T>
{
public:
    ParameterCumuStd(vector<double> &values, ParameterHelper &ph) : Parameter<T>(values, ph) {};
	ParameterCumuStd(int nloci_, ParameterHelper &ph) : Parameter<T>(nloci_, ph) {};
    ParameterCumuStd(const Parameter1s<T> &p) : Parameter<T>((Parameter<T>)p)  {
		from1s(p);
	};


	void	from1s(const Parameter1s<T>& other);
	virtual	const char	*className(void) const { return "CumuStd"; }
};
typedef ParameterCumuStd<parameter_t>	ParameterCumuStdL;

template <class T> void ParameterCumuStd<T>::from1s(const Parameter1s<T>& other) {
	(*this)[0] = 0;	// undefined, there is no haplotype with 0 loci

	// iterate haplotypes, from low to high order
	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			(*this)[h] = other[h];
			continue;
		}
		// <p> minimum/maximum computation
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
// 		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, other[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -other[h2]: other[h2];
			}
			//min1s = max(min1s, (zeros & 1)? r: -r);
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
		}
		// cumulant part
		T			cum = 0;
		// all paritions with 'ones' time components
		const SetPartitionBf	&parts = *(this->c.partitions(ones - 1));
		// iterate over all values with more zeros for which zeros are fixed at positions from h
		// skip first element to make implementations symmetric
		for (int k = 1; k < parts.size(); k++) {
			int	partSize = parts[k].size();
			T	t = (partSize % 2)? tgamma(partSize): -tgamma(partSize);
			for (int l = 0; l < partSize; l++) {
				haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, parts[k][l], ones);
				t *= other[h1];
			}
			cum += t;
		}
		T	cumMin = min1s + cum, cumMax = max1s + cum;
		cum += other[h];
		(*this)[h] = !(cumMax - cumMin)? 1: ((cum - cumMin) / (cumMax - cumMin));
	}

}

template <class T> void Parameter1s<T>::fromCumuStd(const ParameterCumuStd<T>& other) {
	(*this)[0] = 0;	// undefined, there is no haplotype with 0 loci
	// iterate haplotypes, from low to high order
	for (int ones = 1; ones <= this->Nloci; ones++) for (int j = 0; j < this->c.bitOrder[ones].size(); j++) {
		haplotype_t	h = this->c.bitOrder[ones][j];
		// <p> single loci remain unchanged
		if (ones == 1) {
			(*this)[h] = other[h];
			continue;
		}
		// <p> minimum/maximum computation
		// maximum depends on minimum of lower order; start at pairs
		//	we iterate the partitions, it is only necessary to look at haplotypes of length -1
		//	as such data structures are not available, iterate partitions which is a superset of needed haplotypes
		T	min1s = 0, max1s = 1;
		// we only need to iterate partitions with two components as higher components
		//	are implicitely handled
		//	by induction, we iterate over un-standardized values
		// iterate partitions of two, interpret index as binary representation of parition
		for (int l = 0; l < ones; l++) {
			haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, ~(1 << l), ones);
			max1s = min(max1s, (*this)[h1]);
		}
		// <p> determine min by computing multinomial residuals
		for (haplotype_t hm = 0; hm < (1 << ones) -1; hm++) {
			T			r = 0;	// residual for haplotype h
			bitidx_t	zeros = bitCountUnset(hm, ones);	// count number or zeros
			for (int k = 0; k < (1 << zeros) - 1; k++) {
				haplotype_t	h1 = bitEmbedValueInZeros<haplotype_t>(hm, k, zeros);
				bitidx_t	zeros1 = bitCountUnset(h1, ones);	// count number or zeros
				haplotype_t	h2 = bitEmbedValueInOnes<haplotype_t>(h, h1, ones);
				r += ((zeros1 - zeros) & 1)? -(*this)[h2]: (*this)[h2];
				//cout << "   mult-resid-l[" << h2 << "]: " << (*this)[h2] << " r:" << r << endl;
			}
			//min1s = max(min1s, (zeros & 1)? r: -r);
			if (!(zeros & 1)) {
				min1s = max(min1s, max((T)0, -r));
				max1s = min(max1s, 1 - r);
			} else {
				max1s = min(max1s, r);
				min1s = max(min1s, r - 1);
			}
			//cout << "  mult-resid[" << hm << "]: " << r << " (min1s, max1s) = (" << min1s << ", " << max1s << ")" << endl;
		}

		// cumulant part
		T			cum = 0;
		const SetPartitionBf	&parts = *(this->c.partitions(ones - 1));	// all paritions with 'ones' time components
		// iterate over all values with less zeros for which ones are fixed, skip first element
		for (int k = 1; k < parts.size(); k++) {
			int	partSize = parts[k].size();
			T	t = (partSize % 2)? tgamma(partSize): -tgamma(partSize);
			for (int l = 0; l < partSize; l++) {
				haplotype_t	h1 = bitEmbedValueInOnes<haplotype_t>(h, parts[k][l], ones);
				t *= (*this)[h1];
			}
			cum += t;
		}
		T	cumMin = min1s - cum, cumMax = max1s - cum;
		(*this)[h] = other[h] * (cumMax - cumMin) + cumMin;
	}
}

#endif // PARAMETER_H
