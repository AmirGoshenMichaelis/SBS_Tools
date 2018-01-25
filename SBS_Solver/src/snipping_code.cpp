
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

class Settings {
protected:
	json11::Json settings;
public:
	Settings(const json11::Json &);
	virtual Settings & Parse(const json11::Json &);
	virtual const double & operator[] (const std::string &) const;
	virtual double & operator[] (const std::string &);

	virtual const std::string & operator() (const std::string &) const;
	virtual std::string  & operator() (const std::string &);

	virtual const bool & is(const std::string &) const;
	virtual bool & be(const std::string &);

	json11::Json operator[] (const unsigned long i);
	json11::Json at(const std::string &);
};


	std::hash<std::string> hash_fn;
	
	if (hash_fn(key) == hash_fn("h"))
		values["h"] = values["L"] / values["N"];
