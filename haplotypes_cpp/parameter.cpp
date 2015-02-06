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


#include "parameter.h"
#include "parameterhelper.h"


/*
 * 
 */

#if 0
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


