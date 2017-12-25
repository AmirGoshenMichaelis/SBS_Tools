#pragma once

#include <string>
#include "settings_and_param.h"

class SBS_Solver {
	Param  param;
	Config config;

public:
	SBS_Solver(const std::string &);
	void Solve(void);
	void Dump(const std::string &);
};
