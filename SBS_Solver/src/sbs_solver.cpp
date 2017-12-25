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
