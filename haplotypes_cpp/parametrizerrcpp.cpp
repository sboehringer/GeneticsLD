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


#include "parametrizerrcpp.h"

SEXP	ParametrizerRcpp::cumu2multinomial(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::cumu2multinomial(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::cumuStd2multinomial(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::cumuStd2multinomial(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::cumuStd12multinomial(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::cumuStd12multinomial(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}

SEXP	ParametrizerRcpp::multinomial2cumu(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::multinomial2cumu(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::multinomial2cumuStd(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::multinomial2cumuStd(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::multinomial2cumuStd1(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::multinomial2cumuStd1(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::cumuMinMax(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::cumuMinMax(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}

SEXP	ParametrizerRcpp::p1s2multinomial(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::p1s2multinomial(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}

SEXP	ParametrizerRcpp::multinomial2p1s(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::multinomial2p1s(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}

SEXP	ParametrizerRcpp::multinomial2p1sStd(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::multinomial2p1sStd(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
SEXP	ParametrizerRcpp::p1sMinMax(SEXP p) {
	vector<double>	p_ = as< vector<double> >(p);
	vector<double>	*r = Parametrizer::p1sMinMax(p_);
	SEXP	v = wrap(*r);
	delete r;
	return v;
}
