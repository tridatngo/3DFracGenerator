/*
 * getNumLine.hh
 *
 *  Created on: Jan 17, 2017
 *      Author: ngotr
 */

#ifndef INCLUDE_GETNUMLINES_HH_
#define INCLUDE_GETNUMLINES_HH_


#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include <cstring>
#include <iomanip> // setprecision
#include <list>

namespace CGAL{
int getNumLines(const std::string &filename) {

	int numline(0);
	std::string line;
	std::ifstream    infile(filename.c_str());

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file \""<< filename << "\"."<< std::endl;
		//exit(-1);
		return(-1);
	}

	while (!infile.eof()){
		std::getline(infile, line);
		//numline++;
		if (line.size() > 0)
			numline++;
	}

	//std::cout << "numline = " << (numline - 1) << std::endl;

	return numline;
}
} // end of namespace

#endif /* INCLUDE_GETNUMLINES_HH_ */
