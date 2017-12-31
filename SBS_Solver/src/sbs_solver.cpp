#include "sbs_solver.h"
/*
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

*/
//stringify



//class SBS_Solver { // Stimulated Brillouin Scattering
//				   //   Probe E_p;   // probe wave eq.
//				   //   Pump E_s;    // stock wave eq.
//				   //   Density Rho; // rho   wave eq.
//
//	std::vector< std::vector< std::complex<double> > > Ep_mid, Es_mid, Rho_mid;
//	const Param * param;
//	const Config * config;
//	const std::complex<double> coefficient[2];
//	const double a[3][3] = { { 0.0, 0.0, 0.0 },{ 5.0 / 24.0, 1.0 / 3.0, -1.0 / 24.0 },{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } };
//	const double a_hat[3][3] = { { 0.0, 0.0, 0.0 },{ 3.0 / 8.0, 0.0, 9.0 / 8.0 },{ 4.0 / 3.0, -8.0 / 3.0, 10.0 / 3.0 } };
//
//	double SolveEp(unsigned long i, unsigned long j);
//	double SolveEs(unsigned long i, unsigned long j);
//	double SolveRho(unsigned long i, unsigned long j);
//
//public:
//	void Solve(void);
//	std::vector< std::vector< std::complex<double> > > Ep, Es, Rho;
//	std::vector<double>  z, t;
//};
/*

/ *
\partial_\tau E_p = +0.25*alpha*E_p -i*0.5*E_s*rho
\partial_\zeta E_s = -0.25*alpha*E_s +i*g2*E_p*rho^*
\partial_\t rho = -n_fg*c^-1*(0.5*Gamma_B-i*delta_w)rho+n_fg*c^-1*Gamma_B*0.25*g_0*E_p*E_s^*
Ep(t,L) = A(t) Ep(0,z)=0
Es(t,0) = 0    Es(0,z)=0
Rho(0,z) = rho0
* /
class SBS_Solver { // Stimulated Brillouin Scattering
				   //   Probe E_p;   // probe wave eq.
				   //   Pump E_s;    // stock wave eq.
				   //   Density Rho; // rho   wave eq.

	std::vector< std::vector< std::complex<double> > > Ep_mid, Es_mid, Rho_mid;
	const Param * param;
	const Config * config;
	const std::complex<double> coefficient[2];
	const double a[3][3] = { { 0.0, 0.0, 0.0 },{ 5.0 / 24.0, 1.0 / 3.0, -1.0 / 24.0 },{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } };
	const double a_hat[3][3] = { { 0.0, 0.0, 0.0 },{ 3.0 / 8.0, 0.0, 9.0 / 8.0 },{ 4.0 / 3.0, -8.0 / 3.0, 10.0 / 3.0 } };

	double SolveEp(unsigned long i, unsigned long j);
	double SolveEs(unsigned long i, unsigned long j);
	double SolveRho(unsigned long i, unsigned long j);

public:
	SBS_Solver(const Param *, const Config *);
	void Solve(void);
	std::vector< std::vector< std::complex<double> > > Ep, Es, Rho;
	std::vector<double>  z, t;
};


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

double SBS_Solver::SolveRho(unsigned long i, unsigned long j)
{
	// \partial_\t rho = -n_fg*c^-1 * (0.5*Gamma_B - i*delta_w)rho + n_fg*c^-1 * Gamma_B*0.25*g_0*E_p*E_s^*
	// Rho(0, z) = rho0

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

SBS_Solver::SBS_Solver(const Param * iParam, const Config * iConfig) :
	coefficient{ std::complex<double>(-(iParam->n_fg / iParam->c)*(0.5*iParam->Gamma_B), -(iParam->n_fg / iParam->c)*(-iParam->Delta_omega)),
	std::complex<double>(-(iParam->n_fg / iParam->c)*(iParam->Gamma_B*0.25*iParam->g0), 0.0) },
	z(iConfig->N + 1, 0.0),
	t(2 * iConfig->N + 1, 0.0),
	Ep(2 * iConfig->N + 1, std::vector< std::complex<double> >(iConfig->N + 1, std::complex<double>(iParam->E_real, iParam->E_imag))),
	Es(2 * iConfig->N + 1, std::vector< std::complex<double> >(iConfig->N + 1, std::complex<double>(iParam->E_real, iParam->E_imag))),
	Rho(2 * iConfig->N + 1, std::vector< std::complex<double> >(iConfig->N + 1, iParam->rho0)),
	Ep_mid(2 * iConfig->N, std::vector< std::complex<double> >(iConfig->N, std::complex<double>(iParam->E_real, iParam->E_imag))),
	Es_mid(2 * iConfig->N, std::vector< std::complex<double> >(iConfig->N, std::complex<double>(iParam->E_real, iParam->E_imag))),
	Rho_mid(2 * iConfig->N, std::vector< std::complex<double> >(iConfig->N, iParam->rho0))
{
	param = iParam;
	config = iConfig;
	for (unsigned long i = 0; i < 2 * iConfig->N + 1; i++) {
		t[i] = i*iConfig->h;
		Ep[i][config->N] = iParam->A(t[i]);
		Es[i][0] = iParam->A(t[i]) / 1000.0;
	}
	for (unsigned long i = 0; i < iConfig->N + 1; i++) {
		z[i] = i*iConfig->h;
		Ep[0][i] = 0;
		Es[0][i] = 0;
		Rho[0][i] = 0;
	}
}

void SBS_Solver::Solve(void)
{
	std::vector<std::string> convergence_rate;
	std::vector<double> err(3, std::numeric_limits<double>::max());
	double pre_total_err = 0;
	for (unsigned long iteration = 0; iteration < config->max_iteration; iteration++) {
		auto start = std::chrono::system_clock::now();
		double total_err = 0; // std::numeric_limits<double>::max();
		for (unsigned long i = 0; i < 2 * config->N + 1; i++) {
			for (unsigned long j = 0; j < config->N + 1; j++) {
				err[0] = SolveEs(i, j); err[1] = SolveEp(i, j); err[2] = SolveRho(i, j);
				double local_err = *max_element(err.begin(), err.end());
				if (local_err > total_err)
					total_err = local_err;
				//for (unsigned long iteration = 0; max_element(err) < config->max_err && iteration < config->max_iteration; iteration++)
				//err[0] = SolveEp(i, j); err[1] = SolveEs(i, j); err[2] = SolveRho(i, j);
			}
		}
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Iteration #" << iteration << " elapsed time: " << elapsed_seconds.count() << "s (total error=" << total_err << ")\n";
		if (total_err < config->max_err)
			break;
		if (std::fabs(pre_total_err - total_err) < config->max_err)
			break;
		pre_total_err = total_err;
	}
}

bool Read_Json_Param(Config & conf)
{
	const std::string file_name("conf.json");

	std::ifstream json_file(file_name);
	if (!json_file.good())
		return false;

	std::stringstream conf_json;
	conf_json << json_file.rdbuf();
	std::string parse_err;
	json11::Json j = json11::Json::parse(conf_json.str(), parse_err);
	std::cerr << parse_err << std::endl;
	if (j.type() != j.OBJECT)
		return false;
	if (j["L"].is_number())
		conf.L = j["L"].number_value();
	if (j["N"].is_number())
		conf.N = j["N"].number_value();
	if (j["max_err"].is_number())
		conf.max_err = j["max_err"].number_value();
	if (j["max_iteration"].is_number())
		conf.max_iteration = j["max_iteration"].number_value();
	conf.h = conf.L / conf.N;
	return true;
}

int main(int argc, char *argv, char *env)
{
	Param param;
	Config conf;

	Read_Json_Param(conf);

	SBS_Solver sbs(&param, &conf);
	sbs.Solve();

	std::ofstream data("data.json");

	data << "{ ";
	data << R"("t":[)";
	for (unsigned long i = 0; i < 2 * conf.N; ++i)
		data << sbs.t[i] << ",";
	data << sbs.t[2 * conf.N] << "]";
	data << "," << std::endl;
	data << R"("z":[)";
	for (unsigned long i = 0; i < conf.N; ++i)
		data << sbs.z[i] << ",";
	data << sbs.z[conf.N] << "]";
	data << "," << std::endl;
	data << R"("Ep":{)";
	data << R"("real":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Ep[i][j].real() << ",";
		data << sbs.Ep[i][conf.N].real() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "," << std::endl;
	data << R"("imag":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Ep[i][j].imag() << ",";
		data << sbs.Ep[i][conf.N].imag() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "}" << std::endl;

	data << "," << std::endl;
	data << R"("Es":{)";
	data << R"("real":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Es[i][j].real() << ",";
		data << sbs.Es[i][conf.N].real() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "," << std::endl;
	data << R"("imag":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Es[i][j].imag() << ",";
		data << sbs.Es[i][conf.N].imag() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "}" << std::endl;

	data << "," << std::endl;
	data << R"("Rho":{)";
	data << R"("real":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Rho[i][j].real() << ",";
		data << sbs.Rho[i][conf.N].real() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "," << std::endl;
	data << R"("imag":[)";
	for (unsigned long i = 0; i < 2 * conf.N + 1; ++i) {
		data << "[";
		for (unsigned long j = 0; j < conf.N; ++j)
			data << sbs.Rho[i][j].imag() << ",";
		data << sbs.Rho[i][conf.N].imag() << (i == 2 * conf.N ? "]" : " ],") << std::endl;
	}
	data << " ]" << std::endl;
	data << "}" << std::endl;
	data << "}" << std::endl;
	return (0);
}

*/

SBS_Solver::SBS_Solver(const std::string &)
{
}

void SBS_Solver::Solve(void)
{
}

void SBS_Solver::Stringify(const std::string &)
{
}

void SBS_Solver::Dump_H5(const std::string &)
{
}
