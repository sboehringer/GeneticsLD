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


#ifndef PARAMETERHELPER_H
#define PARAMETERHELPER_H
// -- begin inline Rcpp
#include <iostream>
#include <vector>
#include <bitset>

using namespace std;
typedef unsigned int				uint;
typedef unsigned int				element_t;
typedef vector< vector<uint> >		bitOrder_t;
// partion, blocks, elements
typedef vector<uint>				indexSet_t;
typedef vector<indexSet_t>			partitionBlocks_t;
typedef vector<partitionBlocks_t>	SetPartition;


inline void print_setPartition(const SetPartition &p) {
	for (uint i = 0; i < p.size(); i++) {
		cout << "Partition: " << i << ": ";
		for (uint j = 0; j < p[i].size(); j++) {
			if (j) cout << ", ";
			cout << "{";
			for (uint k = 0; k < p[i][j].size(); k++) {
				if (k) cout << ", ";
				cout << p[i][j][k];
			}
			cout << "}";
		}
		cout << "\n";
	}
}
// partitions are interpreted in terms of N + 1 elements, thus partitions[0] partitions a single element
// bitfield based version partion, blocks, elements as bitfield
typedef vector<uint>				partitionBlocksBf_t;
// partition with only one element garuanteed to be first element
typedef vector<partitionBlocksBf_t>	SetPartitionBf;

inline void print_setPartitionBf(const SetPartitionBf &p) {
	for (uint i = 0; i < p.size(); i++) {
		cout << "Partition: " << i << ": ";
		for (uint j = 0; j < p[i].size(); j++) {
			if (j) cout << ", ";
			cout << "{" << bitset<8>(p[i][j]) << "}";
		}
		cout << "\n";
	}
}

typedef	int	bitidx_t;

// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
inline int countBitsSet(unsigned int i) {
    i = i - ((i >> 1) & 0x55555555U);
    i = (i & 0x33333333U) + ((i >> 2) & 0x33333333U);
    return (((i + (i >> 4)) & 0x0F0F0F0FU) * 0x01010101U) >> 24;
}
inline int	bitCountSet(unsigned int i, bitidx_t maxBits) {
	return countBitsSet(i & ( (1 << maxBits) - 1));
}
inline int bitCountUnset(unsigned int i, bitidx_t maxBits) {
	bitidx_t		bitsInT = sizeof(int) * 8;
	unsigned int	mask = ~(unsigned int)0 - ( (1 << maxBits) - 1);
	return bitsInT - countBitsSet(i | mask);
}

template<class T> inline T	bitSet(T v, bitidx_t i) {
	return v | (((T)1) << i);
}
template<class T> inline T	bitClear(T v, bitidx_t i) {
	return v & (~ (((T)1) << i));
}
template<class T> inline T	bitAt(T v, bitidx_t i) {
	return !!(v & (((T)1) << i));
}
template<class T> inline T	bitEmbedValueInZeros(T d, T v, bitidx_t maxBits) {
	//bitidx_t	bitsToSet = bitCountUnset(d, maxBits);

	for (bitidx_t i = 0, j = 0; j < maxBits; i++) {
		if (!bitAt<T>(d, i)) {
			if (bitAt<T>(v, j)) d = bitSet<T>(d, i);
			j++;
		}
	}
	return d;
}
template<class T> inline T	bitEmbedValueInOnes(T d, T v, bitidx_t maxBits) {
	//bitidx_t	bitsToSet = bitCountSet(d, maxBits);

	for (bitidx_t i = 0, j = 0; j < maxBits; i++) {
		if (bitAt<T>(d, i)) {
			if (!bitAt<T>(v, j)) d = bitClear<T>(d, i);
			j++;
		}
	}
	return d;
}

inline void printBo(bitOrder_t &bo) {
	int	N = bo.size() - 1;
	cout << "bitOrder [" << N << "]: \n";
	for (int j = 0; j <= N; j++) {
		cout << "\t#bits [" << j << "] #els: " << bo[j].size() << ": ";
		for (int k = 0; k < bo[j].size(); k++) {
			if (k) cout << ", ";
			cout << bitset<5>(bo[j][k]);
		}
		cout << '\n';
	}
}
inline void computeBo(bitOrder_t &bo) {
	for ( uint j = 0; j < ((uint)1U << (bo.size() - 1)); j++) bo[countBitsSet(j)].push_back(j);
}



class PartitionCache {
public:
	int						Nelements;
	bitOrder_t				bitOrder;
	vector<SetPartitionBf*>	Partitions;

public:
	PartitionCache(int Nelements_)
	: Nelements(Nelements_), bitOrder(Nelements + 1, vector<uint>(0, 0) ), Partitions(Nelements) {
		computeBitorder();
		computePartitions();
	}

	const SetPartitionBf	*partitions(int i) const { return Partitions[i]; }
	void printBitOrder(void) const;
private:
	void computeBitorder(void);
	void computePartitions(void);
	
	
};

typedef  vector<PartitionCache *>	CacheVector;

class ParameterHelper : CacheVector
{
public:
	ParameterHelper(int maxElements_) : CacheVector(maxElements_, (PartitionCache *)0) {}
	const PartitionCache& operator[](int i) const {
		const PartitionCache *	el = (*(CacheVector *)this)[i - 1];
		if (!el) {
			el = (*(CacheVector *)this)[i - 1] = new PartitionCache(i);
		}
		return *el;
	}
	
};
// -- end inline Rcpp

#endif // PARAMETERHELPER_H
