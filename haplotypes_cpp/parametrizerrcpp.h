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


#ifndef PARAMETRIZERRCPP_H
#define PARAMETRIZERRCPP_H

#include "parametrizer.h"
// -- begin inline Rcpp
#include <Rcpp.h>
using namespace Rcpp;

class ParametrizerRcpp : public Parametrizer
{

public:
	ParametrizerRcpp() {}
	SEXP	cumu2multinomial(SEXP p);
	SEXP	multinomial2cumu(SEXP p);
	SEXP	multinomial2cumuStd(SEXP p);
	SEXP	cumuMinMax(SEXP p);

	SEXP	p1s2multinomial(SEXP p);
	SEXP	multinomial2p1s(SEXP p);
	SEXP	multinomial2p1sStd(SEXP p);
	SEXP	p1sMinMax(SEXP p);

};
// -- end inline Rcpp

#endif // PARAMETRIZERRCPP_H
