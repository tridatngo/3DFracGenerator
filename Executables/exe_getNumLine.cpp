/*
 * exe_getNumLine.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: ngotr
 */

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


int main(int argc, char** argv) {

	std::string filename = "";
	if( argc == 2 ) {
		filename = argv[1];
	}
	else {
		std::cout << "Usage: ./exefile filename \n";
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	int numline(0);
	std::string line;
	std::ifstream    infile(filename.c_str());

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file \""<< filename << "\"."<< std::endl;
		exit(-1);
	}

	while (!infile.eof()){
		std::getline(infile, line);
		//numline++;
		if (line.size() > 0)
			numline++;
	}

	std::cout << "numline = " << (numline - 1) << std::endl;

	return 1;
}
