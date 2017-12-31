#pragma once

#include <string>
#include "settings_and_param.h"

class SBS_Solver {
	const static double a[3][3];
	const static double a_hat[3][3];
protected:
	Param  param;
	Config config;

public:
	SBS_Solver(const std::string &);
	void Solve(void);
	void Stringify(const std::string &);
	void Dump_H5(const std::string &);
};
