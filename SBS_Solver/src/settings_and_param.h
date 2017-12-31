#pragma once

#include <string>
#include <map>
#include <complex>
#include "json11.hpp"

namespace Constants {
	constexpr double pi = 3.1415926535897;
	constexpr double c = 3e8;
}

class Config {
	std::map< std::string, double > values;
public:
	Config();
	void Init(const json11::Json &);

	const double operator[] (const std::string &) const;
	double & operator[] (const std::string &);
};

class Param {
	std::map< std::string, double > values;
	std::map< std::string, std::vector< std::complex<double> > > func;
public:
	Param();
	void Init(const json11::Json &);

	const double operator[] (const std::string &) const;
	double & operator[] (const std::string &);

	const std::complex<double> Rho_Initial_Condition(const double x) const;
	const std::complex<double> Es_Boundary_Value(const double t) const;
	const std::complex<double> Ep_Boundary_Value(const double t) const;
};
