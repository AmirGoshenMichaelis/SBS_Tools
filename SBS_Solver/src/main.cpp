
#include <iostream>
#include <fstream>
#include <string>
#include "sbs_solver.h"

int main(int argc, char *argv, char *env)
{
	SBS_Solver sbs("settings.json");
	
	sbs.Solve();
	sbs.Dump("data.json");

	return (0);
}