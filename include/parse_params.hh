/*
 * parse_params.hh
 *
 *  Created on: Sep 15, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_PARSE_PARAMS_HH_
#define INCLUDE_PARSE_PARAMS_HH_

#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<cstring>

namespace CGAL{

int parse_params(int argc, char* argv[]){
	const char* inputFile = (argc > 1) ? argv[1] : "CL_dfnGen.input";
	std::ifstream    infile(inputFile);
	int n4;

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		exit(-1);
	}

	std::string line;
	std::istringstream ist(line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> n4;
	}

	return n4;
}

} // end of namespace
#endif /* INCLUDE_PARSE_PARAMS_HH_ */
