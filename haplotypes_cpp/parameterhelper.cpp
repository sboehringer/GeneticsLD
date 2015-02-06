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


#include "parameterhelper.h"

void PartitionCache::printBitOrder(void) const {
	cout << "bitOrder [" << Nelements << "]: \n";
	for (int j = 0; j <= Nelements; j++) {
		cout << "\t#bits [" << j << "] #els: " << bitOrder[j].size() << ": ";
		for (int k = 0; k < bitOrder[j].size(); k++) {
			if (k) cout << ", ";
			cout << bitOrder[j][k];
		}
		cout << '\n';
	}
}

void PartitionCache::computeBitorder(void) {

	for ( uint j = 0; j < ((uint)1U << Nelements); j++) {
		//cout << "Pushing " << j << " to bucket " << countBitsSet(j) << " [W: " << Nelements << "]\n";
		bitOrder[countBitsSet(j)].push_back(j);
	}
}

// all partitions of set { 0, ..., N }
SetPartition	*partitionSet(uint N) {

	if (N == 0) {
		return new SetPartition(1, partitionBlocks_t(1, indexSet_t(1, (uint)N)));
	} else {
		SetPartition	&p = *partitionSet(N - 1);
		SetPartition	&n = *new SetPartition();	// heuristic capacity
		// copy-step
		for (uint i = 0; i < p.size(); i++) {
			partitionBlocks_t	&b = *new partitionBlocks_t(p[i]);
			b.push_back( indexSet_t(1, (uint)N) );
			n.push_back(b);
		}
		// add-step
		for (uint i = 0; i < p.size(); i++) {
			for (uint j = 0; j < p[i].size(); j++) {
				partitionBlocks_t	&b = *new partitionBlocks_t(p[i]);
				b[j].push_back( N );
				n.push_back(b);
			}
		}
		delete &p;
		return &n;
	}
}
#if 0
// recursive algorithm:
//	new partition is c( partitions one element less + new block with new element, partitions one element less add new element to each block )
void PartitionCache::computePartitions(void) {
	if (Nelements) {
		Partitions = *partitionSet(Nelements - 1);	
	} else {
		// throw
	}
}
#endif

SetPartitionBf	*partitionSetBf(uint N) {
	if (N == 0) {
		return new SetPartitionBf(1, partitionBlocksBf_t(1, bitSet<uint>(0, 0)));
	} else {
		SetPartitionBf	&p = *partitionSetBf(N - 1);
		SetPartitionBf	&n = *new SetPartitionBf();	// heuristic capacity
		// add-step
		for (uint i = 0; i < p.size(); i++) {
			for (uint j = 0; j < p[i].size(); j++) {
				partitionBlocksBf_t	b(p[i]);
				b[j] = bitSet<uint>(b[j], N);
				n.push_back(b);
			}
		}
		// copy-step
		for (uint i = 0; i < p.size(); i++) {
			partitionBlocksBf_t	b(p[i]);
			b.push_back(bitSet<uint>(0, N));
			n.push_back(b);
		}
		delete &p;
		return &n;
	}
}

void PartitionCache::computePartitions(void) {
	if (Nelements) {
		for (int i = 0; i < Nelements; i++) {
			Partitions[i] = partitionSetBf(i);
		}
	} else {
		// throw
	}
}

