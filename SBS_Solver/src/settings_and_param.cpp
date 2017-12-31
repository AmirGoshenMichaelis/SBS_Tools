
#include <map>
#include <vector>
#include <string>
#include "settings_and_param.h"

Config::Config()
{
	json11::Json j = json11::Json::object{ { "L", 50.0 },{ "N", 1e3 },{ "MaxErr", 1e-5 },{ "MaxIteration", 10 } };
	Init(j);
}

void Config::Init(const json11::Json & j)
{
	for (auto it = j.object_items().begin(); it != j.object_items().end(); ++it) {
		if (it->second.is_number())
			values[it->first] = it->second.number_value();
	}
}

const double Config::operator[](const std::string & key) const
{
	return (values.at(key));
}

double & Config::operator[] (const std::string & key)
{
	return (values[key]);
}


//Param::Param(double In_fg = 1.45, double IGamma_B = 2 * pi*40e6, double Ig0 = 5e-11, double IDelta_omega = 0.0, double iAlpha = 0.0) // 2.3e-5
//	: n_fg(In_fg), Gamma_B(IGamma_B), g0(Ig0), Delta_omega(IDelta_omega), rho0(0.0), alpha(iAlpha) {
//	;
//};
//const double n_fg, Gamma_B, g0, Delta_omega;
//const double rho0;
//const double alpha;

