
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <limits>
#include <vector>
#include <chrono>
#include <algorithm>
#include "json11.hpp"
#include "sbs_solver.h"


// Stimulated Brillouin Scattering

const double SBS_Solver::SolveEq::a[3][3] = { { 0.0, 0.0, 0.0 },{ 5.0 / 24.0, 1.0 / 3.0, -1.0 / 24.0 },{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } };
const double SBS_Solver::SolveEq::a_hat[3][3] = { { 0.0, 0.0, 0.0 },{ 3.0 / 8.0, 0.0, 9.0 / 8.0 },{ 4.0 / 3.0, -8.0 / 3.0, 10.0 / 3.0 } };

SBS_Solver::SBS_Solver(const std::string & file_name)
{
	std::ifstream json_file(file_name);
	if (json_file.good()) {
		std::stringstream settings_json;
		settings_json << json_file.rdbuf();
		std::string err;
		json11::Json j = json11::Json::parse(settings_json.str(), err);
		if (!err.empty())
			throw err;
		if (j["configuration"].is_object())
			config.Init(j["configuration"]);
		if (j["parameters"].is_object())
			param.Init(j["parameters"]);
	};
}

void SBS_Solver::Solve(void)
{
	Init_Variables();
	iteration_info.duration_elapsed_time_seconds.clear();
	iteration_info.error.clear();

	double pre_err = 0;
	for (unsigned long iteration = 0; iteration < config["MaxIteration"]; iteration++) {
		auto start = std::chrono::system_clock::now();
		double err = Run_Iteration_Solver();
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		iteration_info.duration_elapsed_time_seconds.push_back(elapsed_seconds.count());
		iteration_info.error.push_back(err);
		std::cout << "Iteration #" << iteration << " elapsed time: " << elapsed_seconds.count() << "s (total error=" << err << ")\n";

		if (
			(err < config["MaxErr"]) || // absolute error
			(std::fabs(pre_err - err) < config["MaxErr"]) || // relative absolute error
			(std::fabs(pre_err - err) / err < config["MaxErr"]) // relative error
			)
			break;
		pre_err = err;
	}
}

void SBS_Solver::Stringify(const std::string & file_name)
{
	std::ofstream data(file_name);

	data << "{ " << std::endl;

	std::map< std::string, Vec * > map_vec{ { "t", &t }, { "z", &z }};
	bool add_comma = false;
	for (auto it : map_vec) {
		if (add_comma)
			data << "," << std::endl;
		else 
			add_comma = true;
		data << "\"" << it.first << "\" : [";
		for (unsigned long i = 0; i < it.second->size(); ++i)
			data << (*it.second)[i] << ((i<it.second->size()-1)? "," : "]");
	}

	data << "," << std::endl;

	std::map< std::string, Matrix * > map_matrix{ { "Ep", &Ep }, { "Es", &Es }, { "Rho", &Rho } };
	add_comma = false;
	for (auto it : map_matrix) {
		if (add_comma)
			data << "," << std::endl;
		else
			add_comma = true;
		data << "\"" << it.first << "\" : {" << std::endl;
		data << "\t\"real\" : \n\t[" << std::endl;
		for (unsigned long i = 0; i < it.second->size(); ++i) {
			data << "\t\t[";
			for (unsigned long j = 0; j < (*it.second)[i].size(); ++j)
				data << (*it.second)[i][j].real() << ((j < (*it.second)[i].size() - 1) ? "," : "");
			data << ((i < it.second->size() - 1) ? "]," : "]") << std::endl;
		}
		data << "\t], " << std::endl << "\t\"imag\" : [" << std::endl;
		for (unsigned long i = 0; i < it.second->size(); ++i) {
			data << "\t\t[";
			for (unsigned long j = 0; j < (*it.second)[i].size(); ++j)
				data << (*it.second)[i][j].imag() << ((j < (*it.second)[i].size() - 1) ? "," : "");
			data << ((i < it.second->size() - 1) ? "]," : "]") << std::endl;
		}
		data << "\t]" << std::endl;
		data << "\t}";
	}
	data << "," << std::endl;
	data << "\"log\" : {" << std::endl;
	map_vec.clear();
	map_vec["duration_elapsed_time_seconds"] = &iteration_info.duration_elapsed_time_seconds;
	map_vec["error"] = &iteration_info.error;
	add_comma = false;
	for (auto it : map_vec) {
		if (add_comma)
			data << "," << std::endl;
		else
			add_comma = true;
		data << "\t\"" << it.first << "\" : [";
		for (unsigned long i = 0; i < it.second->size(); ++i)
			data << (*it.second)[i] << ((i<it.second->size() - 1) ? "," : "]");
	}
	data << "\n\t}" << std::endl;

	data << "}" << std::endl;
}

void SBS_Solver::Dump_H5(const std::string &)
{
}

void SBS_Solver::Init_Variables(void)
{
	config["h"] = config["L"] / config["N"];
	sim_points.mid_pints.spatial = config["N"];
	sim_points.mid_pints.temporal = 2 * config["N"];
	sim_points.points.spatial = config["N"] + 1;
	sim_points.points.temporal = 2 * config["N"] + 1;

	std::vector<Matrix *> mat_ptr = { &Ep, &Es, &Rho };
	for (auto it : mat_ptr)
		(*it) = Matrix(sim_points.points.temporal, std::vector< std::complex<double> >(sim_points.points.spatial, std::complex<double>()));
	std::vector<Matrix *> mat_mid_pint_ptr = { &Ep_mid, &Es_mid, &Rho_mid };
	for (auto it : mat_mid_pint_ptr)
		(*it) = Matrix(sim_points.mid_pints.temporal, std::vector< std::complex<double> >(sim_points.mid_pints.spatial, std::complex<double>()));

	t = Vec(sim_points.points.temporal, 0.0);
	z = Vec(sim_points.points.spatial, 0.0);

	// Boundary condition
	for (unsigned long i = 0; i < sim_points.points.temporal; i++) {
		t[i] = i*config["h"];
		Ep[i][sim_points.points.spatial - 1] = param.Ep_Boundary_Value(t[i]);
		Es[i][0] = param.Es_Boundary_Value(t[i]);
	}
	// Initial condition
	for (unsigned long k = 0; k < sim_points.points.spatial; k++) {
		z[k] = k*config["h"];
		Ep[0][k] = 0.0;
		Es[0][k] = 0.0;
		Rho[0][k] = param.Rho_Initial_Condition(z[k]);
	}

}

double SBS_Solver::Run_Iteration_Solver(void)
{
	double represented_err = 0.0;

	for (unsigned long i = 0; i < sim_points.points.temporal; i++) {
		for (unsigned long k = 0; k < sim_points.points.spatial; k++) {
			std::vector<double> err;
			err.push_back(SolveEs(i, k));
			err.push_back(SolveEp(i, k));
			err.push_back(SolveRho(i, k));
			represented_err += *std::max_element(err.begin(), err.end()) / static_cast<double>(sim_points.points.temporal*sim_points.points.spatial);
		}
	}

	return (represented_err);
}

std::vector< std::complex<double> > SBS_Solver::Get_Grid_Values(const unsigned long i, const unsigned long k)
{
	return (std::vector< std::complex<double> >{Es[i][k], Ep[i][k], Rho[i][k]});
}

std::vector< std::complex<double> > SBS_Solver::Get_Grid_Mid_Values(const unsigned long i, const unsigned long k)
{
	return (std::vector< std::complex<double> >{Es_mid[i][k], Ep_mid[i][k], Rho_mid[i][k]});
}

double SBS_Solver::SolveEs(const unsigned long i, const unsigned long k)
{
	if (k > sim_points.points.spatial - 2 || i > sim_points.points.temporal - 2)
		return (0.0 /* std::numeric_limits<double>::quiet_NaN() */);
	// NumericError err;
	EsEquation Es_eq(param, Matrix{ { Get_Grid_Values(i,k) }, { Get_Grid_Mid_Values(i,k) }, { Get_Grid_Values(i + 1,k + 1) } });
	Es_mid[i][k] = Es[i][k] + config["h"] * Es_eq.Mid();
	Es[i + 1][k + 1] = Es[i][k] + config["h"] * Es_eq.Point();

	if (k <= sim_points.points.spatial - 3 && i <= sim_points.points.temporal - 3) {
		Es_mid[i + 1][k + 1] = Es[i][k] + config["h"] * Es_eq.Extrapolate_Mid();
		Es[i + 2][k + 2] = Es[i][k] + config["h"] * Es_eq.Extrapolate_Point();
	}


	return 0.0;
}

double SBS_Solver::SolveEp(const unsigned long i, const unsigned long k)
{

	if (k < 1 || i > sim_points.points.temporal - 2)
		return (0.0 /* std::numeric_limits<double>::quiet_NaN() */);
	// NumericError err;
	EpEquation Ep_eq(param, Matrix{ { Get_Grid_Values(i,k) },{ Get_Grid_Mid_Values(i,k - 1) },{ Get_Grid_Values(i + 1,k - 1) } });

	Ep_mid[i][k - 1] = Ep[i][k] + config["h"] * Ep_eq.Mid();
	Ep[i + 1][k - 1] = Ep[i][k] + config["h"] * Ep_eq.Point();

	if (k >= 2 && i <= sim_points.points.temporal - 3) {
		Ep_mid[i + 1][k - 2] = Ep[i][k] + config["h"] * Ep_eq.Extrapolate_Mid();
		Ep[i + 2][k - 2] = Ep[i][k] + config["h"] * Ep_eq.Extrapolate_Point();
	}

	return 0.0;
}

double SBS_Solver::SolveRho(const unsigned long i, const unsigned long k)
{
	// NumericError err;
	if (i < sim_points.points.temporal - 1) {
		if (i >= 1) {
			RhoEquation rho_eq(param, Matrix{ { Get_Grid_Values(i,k) },{ Get_Grid_Values(i - 1,k) } });
			Rho[i + 1][k] = Rho[i][k] + config["h"] / 2.0 * rho_eq.Predictor();
		}
		RhoEquation rho_eq(param, Matrix{ { Get_Grid_Values(i + 1,k) },{ Get_Grid_Values(i,k) } });
		Rho[i + 1][k] = Rho[i][k] + config["h"] / 2.0 * rho_eq.Corrector();
	}

	if (i < sim_points.mid_pints.temporal-1 && k < sim_points.mid_pints.spatial) {
		if (i >= 1) {
			RhoEquation rho_eq(param, Matrix{ { Get_Grid_Mid_Values(i,k) },{ Get_Grid_Mid_Values(i - 1,k) } });
			Rho_mid[i + 1][k] = Rho_mid[i][k] + config["h"] / 2.0 * rho_eq.Predictor();
		}
		RhoEquation rho_eq(param, Matrix{ { Get_Grid_Mid_Values(i + 1,k) },{ Get_Grid_Mid_Values(i,k) } });
		Rho_mid[i + 1][k] = Rho_mid[i][k] + config["h"] / 2.0 * rho_eq.Corrector();;
	}

	return (0.0 /* std::numeric_limits<double>::quiet_NaN() */);
}

SBS_Solver::SolveEq::SolveEq(const Matrix & mat) : values(mat)
{
}

std::complex<double> SBS_Solver::SolveEq::Mid(void)
{
	std::complex<double> result;
	result =
		a[1][0] * Eq(values[0][0], values[0][1], values[0][2]) +
		a[1][1] * Eq(values[1][0], values[1][1], values[1][2]) +
		a[1][2] * Eq(values[2][0], values[2][1], values[2][2]);
	return (result);
}

std::complex<double> SBS_Solver::SolveEq::Point(void)
{
	std::complex<double> result;
	result =
		a[2][0] * Eq(values[0][0], values[0][1], values[0][2]) +
		a[2][1] * Eq(values[1][0], values[1][1], values[1][2]) +
		a[2][2] * Eq(values[2][0], values[2][1], values[2][2]);
	return (result);
}

std::complex<double> SBS_Solver::SolveEq::Extrapolate_Mid(void)
{
	std::complex<double> result;
	result =
		a_hat[1][0] * Eq(values[0][0], values[0][1], values[0][2]) +
		a_hat[1][1] * Eq(values[1][0], values[1][1], values[1][2]) +
		a_hat[1][2] * Eq(values[2][0], values[2][1], values[2][2]);
	return (result);
}

std::complex<double> SBS_Solver::SolveEq::Extrapolate_Point(void)
{
	std::complex<double> result;
	result =
		a_hat[2][0] * Eq(values[0][0], values[0][1], values[0][2]) +
		a_hat[2][1] * Eq(values[1][0], values[1][1], values[1][2]) +
		a_hat[2][2] * Eq(values[2][0], values[2][1], values[2][2]);
	return (result);
}

std::complex<double> SBS_Solver::EsEquation::Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho)
{
	std::complex<double> value;
	value = -0.25*param["alpha"] * Es + std::complex<double>(0.0, 1) * Ep * std::conj(Rho);
	return (value);
}

SBS_Solver::EsEquation::EsEquation(const Param & parameters, const Matrix & mat) : SolveEq(mat), param(parameters)
{

}

std::complex<double> SBS_Solver::EpEquation::Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho)
{
	std::complex<double> value;
	value = +0.25*param["alpha"] * Ep + std::complex<double>(0.0, 1) * Es * Rho;
	return (value);
}

SBS_Solver::EpEquation::EpEquation(const Param & parameters, const Matrix & mat) : SolveEq(mat), param(parameters)
{

}

std::complex<double> SBS_Solver::RhoEquation::Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho)
{
	std::complex<double> value;
	value = Rho * std::complex<double>(-(param["n_fg"] / Constants::c)*(0.5*param["Gamma_B"]), -(param["n_fg"] / Constants::c)*(-param["Delta_omega"])) +
		std::complex<double>(-(param["n_fg"] / Constants::c)*(param["Gamma_B"] * 0.25*param["g0"]), 0.0) * Ep * std::conj(Es);
	return (value);
}

SBS_Solver::RhoEquation::RhoEquation(const Param & parameters, const Matrix & mat) : values(mat), param(parameters)
{

}

std::complex<double> SBS_Solver::RhoEquation::Corrector(void)
{
	std::complex<double> result;
	result = Eq(values[0][0], values[0][1], values[0][2]) +
		Eq(values[1][0], values[1][1], values[1][2]);
	return result;
}

std::complex<double> SBS_Solver::RhoEquation::Predictor(void)
{
	std::complex<double> result;
	result = 3.0 * Eq(values[0][0], values[0][1], values[0][2]) -
		Eq(values[1][0], values[1][1], values[1][2]);
	return result;
}

/*
double SBS_Solver::SolveEp(unsigned long i, unsigned long j)
{
// \partial_\tau E_p = +0.25*alpha*E_p - i*0.5*E_s*rho
// Ep(t, 0) = A(t) Ep(0, z) = 0
std::vector<double> err(2, std::numeric_limits<double>::max());
std::complex<double> ypre, ycurr;
if (j < 1 || i > 2 * config->N - 1)
return 0.0;
ypre = Ep_mid[i][j - 1];
Ep_mid[i][j - 1] = Ep[i][j] + config->h *(
a[1][0] * (+0.25*param->alpha*Ep[i][j] - std::complex<double>(0.0, 1)*Es[i][j] * Rho[i][j]) +
a[1][1] * (+0.25*param->alpha*Ep_mid[i][j - 1] - std::complex<double>(0.0, 1)*Es_mid[i][j - 1] * Rho_mid[i][j - 1]) +
a[1][2] * (+0.25*param->alpha*Ep[i + 1][j - 1] - std::complex<double>(0.0, 1)*Es[i + 1][j - 1] * Rho[i + 1][j - 1])
);
ycurr = Ep_mid[i][j - 1];
if (std::abs(ycurr) == 0.0)
err[0] = std::abs(ycurr - ypre);
else
err[0] = std::abs(ycurr - ypre) / std::abs(ycurr);

ypre = Ep[i + 1][j - 1];
Ep[i + 1][j - 1] = Ep[i][j] + config->h *(
a[2][0] * (+0.25*param->alpha*Ep[i][j] - std::complex<double>(0.0, 1)*Es[i][j] * Rho[i][j]) +
a[2][1] * (+0.25*param->alpha*Ep_mid[i][j - 1] - std::complex<double>(0.0, 1)*Es_mid[i][j - 1] * Rho_mid[i][j - 1]) +
a[2][2] * (+0.25*param->alpha*Ep[i + 1][j - 1] - std::complex<double>(0.0, 1)*Es[i + 1][j - 1] * Rho[i + 1][j - 1])
);
ycurr = Ep[i + 1][j - 1];
if (std::abs(ycurr) == 0.0)
err[1] = std::abs(ycurr - ypre);
else
err[1] = std::abs(ycurr - ypre) / std::abs(ycurr);

if (j >= 2 && i <= 2 * config->N - 2) {
Ep_mid[i + 1][j - 2] = Ep[i][j] + config->h *(
a_hat[1][0] * (+0.25*param->alpha*Ep[i][j] - std::complex<double>(0.0, 1)*Es[i][j] * Rho[i][j]) +
a_hat[1][1] * (+0.25*param->alpha*Ep_mid[i][j - 1] - std::complex<double>(0.0, 1)*Es_mid[i][j - 1] * Rho_mid[i][j - 1]) +
a_hat[1][2] * (+0.25*param->alpha*Ep[i + 1][j - 1] - std::complex<double>(0.0, 1)*Es[i + 1][j - 1] * Rho[i + 1][j - 1])
);

Ep[i + 2][j - 2] = Ep[i][j] + config->h *(
a_hat[2][0] * (+0.25*param->alpha*Ep[i][j] - std::complex<double>(0.0, 1)*Es[i][j] * Rho[i][j]) +
a_hat[2][1] * (+0.25*param->alpha*Ep_mid[i][j - 1] - std::complex<double>(0.0, 1)*Es_mid[i][j - 1] * Rho_mid[i][j - 1]) +
a_hat[2][2] * (+0.25*param->alpha*Ep[i + 1][j - 1] - std::complex<double>(0.0, 1)*Es[i + 1][j - 1] * Rho[i + 1][j - 1])
);
}

return (*std::max_element(err.begin(), err.end()));
}

*/
/*
double SBS_Solver::SolveEs(unsigned long i, unsigned long j)
{
// \partial_\zeta E_s = -0.25*alpha*E_s + i*g2*E_p*rho^*
// Es(t, L) = 0    Es(0, z) = 0

std::vector<double> err(2, std::numeric_limits<double>::max());
std::complex<double> ypre, ycurr;
if (j > config->N - 1 || i > 2 * config->N - 1)
return 0.0;
ypre = Es_mid[i][j];
Es_mid[i][j] = Es[i][j] + config->h *(
a[1][0] * (-0.25*param->alpha*Es[i][j] + std::complex<double>(0.0, 1)*Ep[i][j] * Rho[i][j]) +
a[1][1] * (-0.25*param->alpha*Es_mid[i][j] + std::complex<double>(0.0, 1)*Ep_mid[i][j] * Rho_mid[i][j]) +
a[1][2] * (-0.25*param->alpha*Es[i + 1][j + 1] + std::complex<double>(0.0, 1)*Ep[i + 1][j + 1] * Rho[i + 1][j + 1])
);
ycurr = Ep_mid[i][j];
if (std::abs(ycurr) == 0.0)
err[0] = std::abs(ycurr - ypre);
else
err[0] = std::abs(ycurr - ypre) / std::abs(ycurr);

ypre = Es[i + 1][j + 1];
Es[i + 1][j + 1] = Es[i][j] + config->h *(
a[2][0] * (-0.25*param->alpha*Es[i][j] + std::complex<double>(0.0, 1)*Ep[i][j] * Rho[i][j]) +
a[2][1] * (-0.25*param->alpha*Es_mid[i][j] + std::complex<double>(0.0, 1)*Ep_mid[i][j] * Rho_mid[i][j]) +
a[2][2] * (-0.25*param->alpha*Es[i + 1][j + 1] + std::complex<double>(0.0, 1)*Ep[i + 1][j + 1] * Rho[i + 1][j + 1])
);
ycurr = Ep[i + 1][j + 1];
if (std::abs(ycurr) == 0.0)
err[1] = std::abs(ycurr - ypre);
else
err[1] = std::abs(ycurr - ypre) / std::abs(ycurr);

if (j <= config->N - 2 && i <= 2 * config->N - 2) {
Es_mid[i + 1][j + 1] = Es[i][j] + config->h *(
a_hat[1][0] * (-0.25*param->alpha*Es[i][j] + std::complex<double>(0.0, 1)*Ep[i][j] * Rho[i][j]) +
a_hat[1][1] * (-0.25*param->alpha*Es_mid[i][j] + std::complex<double>(0.0, 1)*Ep_mid[i][j] * Rho_mid[i][j]) +
a_hat[1][2] * (-0.25*param->alpha*Es[i + 1][j + 1] + std::complex<double>(0.0, 1)*Ep[i + 1][j + 1] * Rho[i + 1][j + 1])
);

Es[i + 2][j + 2] = Es[i][j] + config->h *(
a_hat[2][0] * (-0.25*param->alpha*Es[i][j] + std::complex<double>(0.0, 1)*Ep[i][j] * std::conj(Rho[i][j])) +
a_hat[2][1] * (-0.25*param->alpha*Es_mid[i][j] + std::complex<double>(0.0, 1)*Ep_mid[i][j] * std::conj(Rho_mid[i][j])) +
a_hat[2][2] * (-0.25*param->alpha*Es[i + 1][j + 1] + std::complex<double>(0.0, 1)*Ep[i + 1][j + 1] * std::conj(Rho[i + 1][j + 1]))
);
}

return (*std::max_element(err.begin(), err.end()));
}

*/
/*
double SBS_Solver::SolveRho(unsigned long i, unsigned long j)
{
// \partial_\t rho = -n_fg*c^-1 * (0.5*Gamma_B - i*delta_w)rho + n_fg*c^-1 * Gamma_B*0.25*g_0*E_p*E_s^*
// Rho(0, z) = rho0
coefficient{ std::complex<double>(-(iParam->n_fg / iParam->c)*(0.5*iParam->Gamma_B), -(iParam->n_fg / iParam->c)*(-iParam->Delta_omega)),
std::complex<double>(-(iParam->n_fg / iParam->c)*(iParam->Gamma_B*0.25*iParam->g0), 0.0) },

/ *
predictor (i,j) (i-1,j) --> (i+1,j)
corrector (i,j) (i+1,j) --> (i+1,j)
*********
predictor(mid) (i,j) (i-1,j) --> (i+1,j)
corrector(mid) (i,j) (i+1,j) --> (i+1,j)
* /

std::vector<double> err(2, std::numeric_limits<double>::max());
std::complex<double> ypre, ycurr;

if (i < 2 * config->N) {
ypre = Rho[i + 1][j];
if (i >= 1)
Rho[i + 1][j] = Rho[i][j] + config->h / 2 * (
3.0*(coefficient[0] * Rho[i][j] + coefficient[1] * Ep[i][j] * std::conj(Es[i][j]))
- (coefficient[0] * Rho[i - 1][j] + coefficient[1] * Ep[i - 1][j] * std::conj(Es[i - 1][j])));

Rho[i + 1][j] = Rho[i][j] + config->h / 2 * (
(coefficient[0] * Rho[i + 1][j] + coefficient[1] * Ep[i + 1][j] * std::conj(Es[i + 1][j]))
+ (coefficient[0] * Rho[i][j] + coefficient[1] * Ep[i][j] * std::conj(Es[i][j])));
ycurr = Rho[i + 1][j];
if (std::abs(ycurr) == 0.0)
err[0] = std::abs(ycurr - ypre);
else
err[0] = std::abs(ycurr - ypre) / std::abs(ycurr);
}
else
err[0] = 0;

if (i < 2 * config->N - 1 && j < config->N) {
ypre = Rho_mid[i + 1][j];
if (i >= 1)
Rho_mid[i + 1][j] = Rho_mid[i][j] + config->h / 2 * (
3.0*(coefficient[0] * Rho_mid[i][j] + coefficient[1] * Ep_mid[i][j] * std::conj(Es_mid[i][j]))
- (coefficient[0] * Rho_mid[i - 1][j] + coefficient[1] * Ep_mid[i - 1][j] * std::conj(Es_mid[i - 1][j])));

Rho_mid[i + 1][j] = Rho_mid[i][j] + config->h / 2 * (
(coefficient[0] * Rho_mid[i + 1][j] + coefficient[1] * Ep_mid[i + 1][j] * std::conj(Es_mid[i + 1][j]))
+ (coefficient[0] * Rho_mid[i][j] + coefficient[1] * Ep_mid[i][j] * std::conj(Es_mid[i][j])));
ycurr = Rho_mid[i + 1][j];
if (std::abs(ycurr) == 0.0)
err[1] = std::abs(ycurr - ypre);
else
err[1] = std::abs(ycurr - ypre) / std::abs(ycurr);
}
else
err[1] = 0;

return (*std::max_element(err.begin(), err.end()));
}

*/
