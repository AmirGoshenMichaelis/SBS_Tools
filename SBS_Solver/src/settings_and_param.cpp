
#include <map>
#include <vector>
#include <string>
#include <complex>
#include <random>
#include <limits>
#include <cmath>
#include "json11.hpp"
#include "settings_and_param.h"

Config::Config()
{
	json11::Json j = json11::Json::object{ { "L", 50.0 },{ "N", 1e3 },{ "MaxErr", 1e-5 },{ "MaxIteration", 10 }, { "InerGridPointIteration", 3 }, { "TimeLength", 2} };
	Init(j);
}

void Config::Init(const json11::Json & j)
{
	for (auto it = j.object_items().begin(); it != j.object_items().end(); ++it)
		if (it->second.is_number())
			values[it->first] = it->second.number_value();
}

const double Config::operator[](const std::string & key) const
{
	return (values.at(key));
}

double & Config::operator[] (const std::string & key)
{
	return (values[key]);
}

Param::Param()
{
	json11::Json j = json11::Json::object{
		{ "n_fg", 1.45 },
		{ "Gamma_B", 2 * Constants::pi * 40e6 },
		{ "g0", 5e-11 },
		{ "Delta_omega", 0.0 },
		{ "alpha", 2.3e-5 },
		{ "rho", json11::Json::object{ { "imag" , 0.0 }, {"real", 0.0 } } },
		{ "Ep",  json11::Json::object{ { "imag" , 0.0 }, {"real", 0.0 } } },
		{ "Es",  json11::Json::object{ { "imag" , 0.0 }, {"real", 0.0 } } },
	};
	Init(j);
}

void Param::Init(const json11::Json & j)
{
	for (auto it = j.object_items().begin(); it != j.object_items().end(); ++it)
		if (it->second.is_number())
			values[it->first] = it->second.number_value();
	std::vector< std::string > func_names{ "rho", "Ep", "Es" };
	for (auto name : func_names) {
		//func[name] = std::vector< std::complex<double> >{ std::complex<double>(0.0,0.0) };
		if (j[name].is_object()) {
			size_t size = 0;
			std::map< std::string, std::vector<double> > map_vec{ {"imag",std::vector<double>() }, {"real",std::vector<double>() } };
			for (auto it = map_vec.begin(); it != map_vec.end(); ++it) {
				if (j[name][it->first].is_number())
					it->second.push_back(j[name][it->first].number_value());
				else if (j[name][it->first].is_array())
					for (auto array_valye : j[name][it->first].array_items())
						it->second.push_back(array_valye.is_number() ? array_valye.number_value() : 0.0);
				if(it->second.size() > size)
					size = it->second.size();
			}
			func[name].clear();
			for (unsigned long ele = 0; ele < size; ++ele) {
				double b = (ele < map_vec["imag"].size()) ? map_vec["imag"][ele] : 0.0;
				double a = (ele < map_vec["real"].size()) ? map_vec["real"][ele] : 0.0;
				func[name].push_back(std::complex<double>(a, b));
			}
		}
	}
}

const double Param::operator[](const std::string & key) const
{
	return (values.at(key));
}

double & Param::operator[](const std::string & key)
{
	return (values[key]);
}

const std::complex<double> Param::Boundary_Initial_Value(const std::string &name, const unsigned long ind) const
{
	std::complex<double> value;

	if ((func.at(name).size() == 0) || (func.at(name).size() >= ind))
		return (value);
	if (func.at(name).size() == 1)
		return func.at(name)[0];
	value = func.at(name)[ind];
	return (value);
}

const std::complex<double> Param::Rho_Initial_Condition(const unsigned long ind, const double x) const
{
	//static std::random_device rd;
	//static std::mt19937 gen(rd());
	//static std::uniform_real_distribution<double> dis(-1.0, 1.0);
	//return func.at("rho")[0]*std::conj(func.at("rho")[0])* std::complex<double>(dis(gen), dis(gen));

	return Boundary_Initial_Value(std::string("rho"), ind);
}

const std::complex<double> Param::Es_Boundary_Value(const unsigned long ind, const double t) const
{
	double FWHM = 120e-09 * Constants::c;
	double sigma = FWHM / (2.0*std::sqrt(2.0*std::log(2.0)) );
	double mu = 3.0 * sigma;
	double lambda = func.at("Es")[0].real(); // 1/std::sqrt( 2.0 * Constants::pi * sigma*sigma )
	double value = lambda*std::exp(- (t-mu)*(t - mu) / (2.0*sigma*sigma));
	return std::complex<double>(value, 0.0);
	return Boundary_Initial_Value(std::string("Es"), ind);
}

const std::complex<double> Param::Ep_Boundary_Value(const unsigned long ind, const double t) const
{
	return Boundary_Initial_Value(std::string("Ep"), ind);
}

