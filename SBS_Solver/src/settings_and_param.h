#pragma once

#include <string>
#include <map>
#include <complex>
#include "json11.hpp"

class Config {
	std::map<std::string, double> map;

public:
	double & operator [](const std::string &);
	const double operator[](const std::string &) const;

	Config();
	Config(const json11::Json &);
};

class Param {
	std::map<std::string, double> map;

public:
	double & operator [](const std::string &);
	const double operator[](const std::string &) const;

	static const double pi;
	static const double c;
	
	Param();
	Param(const json11::Json &);

	const std::complex<double> Es_boundary_value(const double t) const;
	const std::complex<double> Ep_boundary_value(const double t) const;

};
const double Param::pi = 3.1415926535897;
const double Param::c = 3e8;
