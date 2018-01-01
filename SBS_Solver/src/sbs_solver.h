#pragma once

#include <string>
#include "settings_and_param.h"

class SBS_Solver {
	struct NumSimPoints {
		struct SpatialTemporalPoints {
			unsigned long spatial, temporal;
		};
		SpatialTemporalPoints mid_pints, points;
	};
	typedef std::vector< std::vector< std::complex<double> > > Matrix;
	typedef std::vector<double> Vec;
	struct IterationInfo {
		std::vector<double> duration_elapsed_time_seconds;
		std::vector<double> error;
	};

	class SolveEq {
	protected:
		const static double a[3][3];
		const static double a_hat[3][3];
		virtual std::complex<double> Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho) = 0;
		Matrix values;
	public:
		SolveEq(const Matrix &);
		std::complex<double> Mid(void);
		std::complex<double> Point(void);
		std::complex<double> Extrapolate_Mid(void);
		std::complex<double> Extrapolate_Point(void);

	};
	class EsEquation : public SolveEq {
		const Param & param;
		virtual std::complex<double> Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho);
	public:
		EsEquation(const Param &, const Matrix &);
	};
	class EpEquation : public SolveEq {
		const Param & param;
		virtual std::complex<double> Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho);
	public:
		EpEquation(const Param &, const Matrix &);
	};
	class RhoEquation {
		const Param & param;
		virtual std::complex<double> Eq(std::complex<double> Es, std::complex<double> Ep, std::complex<double> Rho);
		Matrix values;
	public:
		RhoEquation(const Param &, const Matrix &);
		std::complex<double> Corrector(void);
		std::complex<double> Predictor(void);
	};

	IterationInfo iteration_info;
	NumSimPoints sim_points;

	Matrix Ep_mid, Es_mid, Rho_mid;
	Matrix Ep, Es, Rho;
	Vec  z, t;

	void Init_Variables(void);
	double Run_Iteration_Solver(void);

	std::vector< std::complex<double> > Get_Grid_Mid_Values(const unsigned long, const unsigned long);
	std::vector< std::complex<double> > Get_Grid_Values(const unsigned long, const unsigned long);

	double SolveEs(const unsigned long, const unsigned long);
	double SolveEp(const unsigned long, const unsigned long);
	double SolveRho(const unsigned long, const unsigned long);
protected:
	Param  param;
	Config config;
public:
	SBS_Solver(const std::string &);
	void Solve(void);
	void Stringify(const std::string &);
	void Dump_H5(const std::string &);
};
