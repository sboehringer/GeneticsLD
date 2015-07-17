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


#include "parametrizer.h"

vector<double>	*Parametrizer::multinomial2cumu(vector<double> &i) {
	const ParameterMultinomialL	pMultinom(i, h);
	const Parameter1sL			p1s(pMultinom);
	const ParameterCumuL		pCumu(p1s);
	return vectorConvert<parameter_t, double>(pCumu);
}

vector<double>	*Parametrizer::multinomial2p1sStd(vector<double> &i) {
	const ParameterMultinomialL	pMultinom(i, h);
	const Parameter1sL			p1s(pMultinom);
	const Parameter1sStdL		p1sStd(p1s);
	return vectorConvert<parameter_t, double>(p1sStd);
}


vector<double>	*Parametrizer::multinomial2cumuStd(vector<double> &i) {
	const ParameterMultinomialL	pMultinom(i, h);
	const Parameter1sL			p1s(pMultinom);
	const ParameterCumuStdL		pCumuStd(p1s);
	return vectorConvert<parameter_t, double>(pCumuStd);
}
vector<double>	*Parametrizer::multinomial2cumuStd1(vector<double> &i) {
	const ParameterMultinomialL	pMultinom(i, h);
	const Parameter1sL			p1s(pMultinom);
	const ParameterCumuStd1L	pCumuStd(p1s);
	return vectorConvert<parameter_t, double>(pCumuStd);
}
vector<double>	*Parametrizer::cumuMinMax(vector<double> &i) {
	const ParameterCumuL	pCumu(i, h);
	ParameterCumuMinMaxL 	l = pCumu.minMax();	//limits
	vector<double>			&mn = *vectorConvert<parameter_t, double>(l.min);
	vector<double>			&mx = *vectorConvert<parameter_t, double>(l.max);
	mn.insert(mn.end(), mx.begin(), mx.end());
	delete &mx;
	return &mn;
}


vector<double>	*Parametrizer::cumu2multinomial(vector<double> &i) {
	const ParameterCumuL		pCumu(i, h);
	const Parameter1sL			p1s(pCumu);
	const ParameterMultinomialL	pMultinom(p1s);
	return vectorConvert<parameter_t, double>(pMultinom);
}

vector<double>	*Parametrizer::cumuStd2multinomial(vector<double> &i) {
	const ParameterCumuStdL		pCumuStd(i, h);
	const Parameter1sL			p1s(pCumuStd);
	const ParameterMultinomialL	pMultinom(p1s);
	return vectorConvert<parameter_t, double>(pMultinom);
}

vector<double>	*Parametrizer::cumuStd12multinomial(vector<double> &i) {
	const ParameterCumuStd1L	pCumuStd(i, h);
	const Parameter1sStd1L		p1s(pCumuStd);
	const ParameterMultinomialL	pMultinom(p1s);
	return vectorConvert<parameter_t, double>(pMultinom);
}

vector<double>	*Parametrizer::multinomial2p1s(vector<double> &i) {
	const ParameterMultinomialL	pMultinom(i, h);
	const Parameter1sL			p1s(pMultinom);
	return vectorConvert<parameter_t, double>(p1s);
}
vector<double>	*Parametrizer::p1sMinMax(vector<double> &i) {
	const Parameter1sL	p1s(i, h);
	Parameter1sMinMaxL 	l = p1s.minMax();	//limits
	vector<double>		&mn = *vectorConvert<parameter_t, double>(l.min);
	vector<double>		&mx = *vectorConvert<parameter_t, double>(l.max);
	mn.insert(mn.end(), mx.begin(), mx.end());
	delete &mx;
	return &mn;
}

vector<double>	*Parametrizer::p1s2multinomial(vector<double> &i) {
	const Parameter1sL			p1s(i, h);
	const ParameterMultinomialL	pMultinom(p1s);
	return vectorConvert<parameter_t, double>(pMultinom);
}

vector<double>	*Parametrizer::p1s2cumu(vector<double> &i) {
	const Parameter1sL			p1s(i, h);
	const ParameterCumuL		pCumu(p1s);
	return vectorConvert<parameter_t, double>(pCumu);
}
vector<double>	*Parametrizer::p1s2cumuStd(vector<double> &i) {
	const Parameter1sL			p1s(i, h);
	const ParameterCumuStdL		pCumuStd(p1s);
	return vectorConvert<parameter_t, double>(pCumuStd);
}

vector<double>	*Parametrizer::cumu2p1s(vector<double> &i) {
	const ParameterCumuL		pCumu(i, h);
	const Parameter1sL			p1s(pCumu);
	return vectorConvert<parameter_t, double>(p1s);
}

vector<double>	*Parametrizer::cumuStd2p1s(vector<double> &i) {
	const ParameterCumuStdL		pCumuStd(i, h);
	const Parameter1sL			p1s(pCumuStd);
	return vectorConvert<parameter_t, double>(p1s);
}
