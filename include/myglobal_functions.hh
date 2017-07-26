/*
 * myglobal_functions.hh
 *
 *  Created on: Sep 5, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_MYGLOBAL_FUNCTIONS_HH_
#define INCLUDE_MYGLOBAL_FUNCTIONS_HH_

#include <string>     // std::string, std::to_string
#include <sstream>
#include <iostream>
#include <iomanip>

#define MyPI_ 3.141592653589793

namespace CGAL{

// Radian to degree angle conversion
template <typename NT>
NT Rad2Deg (NT & angle) {
	return angle*180.0/MyPI_;
}

// degree to radian angle conversion
template <typename NT>
NT Deg2Rad (NT & angle) {
	return angle*MyPI_/180.0;
}

std::string print_seconds_to_hours (double &secs){
	std::string out;
	std::ostringstream hours_ss, minutes_ss, seconds_ss;
	int hours, minutes, seconds;

	seconds = int(secs) % 60;
	minutes = int(secs / 60) % 60;
	hours 	= int(seconds / 3600);

	hours_ss << hours;
	minutes_ss << minutes;
	seconds_ss << seconds;

	out = hours_ss.str() + "h" + minutes_ss.str() + "m" + seconds_ss.str() + "s";

	return out;
}

std::string int2str (const int &a){
	std::string out;
	std::ostringstream out_ss;
	out_ss << a;

	out = out_ss.str();

	return out;
}

std::string int2str_setw (const int &a,const int &w){
	std::string out;
	std::ostringstream out_ss;
	out_ss << std::setfill('0') << std::setw(w) << a;
	out = out_ss.str();

	return out;
}

int count_digits (const int &the_integer){
	int nDigits = std::floor(std::log10(std::abs(the_integer))) + 1;

	return nDigits;
}

int getnumline (const std::string filename){

	int numline(0);
	std::string line;
	std::ifstream    infile(filename.c_str());

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file \""<< filename << "\"."<< std::endl;
		exit(-1);
	}

	while (!infile.eof()){
		std::getline(infile, line);
		numline++;
	}
	return (numline - 1);
}

double pround(double x, int precision)
{
    if (x == 0.)
        return x;
    int ex = floor(log10(abs(x))) - precision + 1;
    double div = pow(10, ex);
    return floor(x / div + 0.5) * div;
}

} // end of namespace


#endif /* INCLUDE_MYGLOBAL_FUNCTIONS_HH_ */
