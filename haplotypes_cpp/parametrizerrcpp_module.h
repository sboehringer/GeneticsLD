/*
 * parametrizerrcpp_module.h
 * Mon Sep 30 12:28:59 2013
 */

// -- begin inline Rcpp
RCPP_MODULE(module_parametrizer) {
	using namespace Rcpp;
	class_<ParametrizerRcpp>("parametrizer")
	.constructor()
	.method("cumu2multinomial", &ParametrizerRcpp::cumu2multinomial)
	.method("multinomial2cumu", &ParametrizerRcpp::multinomial2cumu)
	.method("multinomial2cumuStd", &ParametrizerRcpp::multinomial2cumuStd)
	.method("cumuMinMax", &ParametrizerRcpp::cumuMinMax)
	.method("p1s2multinomial", &ParametrizerRcpp::p1s2multinomial)
	.method("multinomial2p1s", &ParametrizerRcpp::multinomial2p1s)
	.method("multinomial2p1sStd", &ParametrizerRcpp::multinomial2p1sStd)
	.method("p1sMinMax", &ParametrizerRcpp::p1sMinMax)
	;
}

// -- end inline Rcpp
