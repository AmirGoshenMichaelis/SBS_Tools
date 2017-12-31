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
public:
	Param();
	Init(const json11::Json &);

	const std::complex<double> Es_boundary_value(const double t) const;
	const std::complex<double> Ep_boundary_value(const double t) const;
};

/*
std::map_values<std::string, double>

*/