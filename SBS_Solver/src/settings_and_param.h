#pragma once

#include <string>
#include <map>
#include <complex>
#include <limits>
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
protected:
	const std::complex<double> Boundary_Initial_Value(const std::string &name, const unsigned long ind) const;
public:
	Param();
	void Init(const json11::Json &);

	const double operator[] (const std::string &) const;
	double & operator[] (const std::string &);

	const std::complex<double> Rho_Initial_Condition(const unsigned long ind, const double x = std::numeric_limits<double>::quiet_NaN()) const;
	const std::complex<double> Es_Boundary_Value(const unsigned long ind, const double t = std::numeric_limits<double>::quiet_NaN()) const;
	const std::complex<double> Ep_Boundary_Value(const unsigned long ind, const double t = std::numeric_limits<double>::quiet_NaN()) const;
};
