
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <map>
#include "sbs_solver.h"

int main(int argc, char **argv, char **env)
{
	std::map<std::string, std::string> opt{ {"Format", "json"}, {"Data", "data.json"}, {"Settings", "settings.json"} };

	for (int i_arg = 1; i_arg < argc; ++i_arg)
		for (auto it = opt.begin(); it != opt.end(); ++it)
			if (std::strcmp(argv[i_arg], ("-" + it->first).c_str())==0)
				if (i_arg + 1 < argc) {
					it->second = argv[++i_arg];
					break;
				}
				else {
					std::cerr << it->first << " is missing argument" << std::endl;
					return (-1);
				}

	SBS_Solver sbs(opt["settings"]);

	sbs.Solve();

	if (_strcmpi(opt["Format"].c_str(), "json")==0)
		sbs.Stringify(opt["Data"]);
	else if (_strcmpi(opt["Format"].c_str(), "hdf5")==0)
		sbs.Dump_H5(opt["Data"]);

	return (0);
}
