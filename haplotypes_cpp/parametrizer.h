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


#ifndef PARAMETRIZER_H
#define PARAMETRIZER_H

#include "parameter.h"

// -- begin inline Rcpp

#define	MAX_LOCI	32
class Parametrizer
{
private:
	ParameterHelper	h;
public:
	Parametrizer() : h(MAX_LOCI) {}
	vector<double>	*cumu2multinomial(vector<double> &i);
	vector<double>	*multinomial2cumu(vector<double> &i);
	vector<double>	*cumuMinMax(vector<double> &i);

	vector<double>	*multinomial2p1s(vector<double> &i);
	vector<double>	*multinomial2p1sStd(vector<double> &i);
	vector<double>	*p1s2multinomial(vector<double> &i);
	vector<double>	*p1sMinMax(vector<double> &i);
	
	vector<double>	*cumu2p1s(vector<double> &i);
	vector<double>	*p1s2cumu(vector<double> &i);

	vector<double>	*multinomial2cumuStd(vector<double> &i);
	vector<double>	*cumuStd2multinomial(vector<double> &i);

	vector<double>	*p1s2cumuStd(vector<double> &i);
	vector<double>	*cumuStd2p1s(vector<double> &i);
};
// -- end inline Rcpp

#endif // PARAMETRIZER_H
