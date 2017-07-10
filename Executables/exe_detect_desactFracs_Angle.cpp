/*
 * exe_create_parameter_mlgi_files.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: ngotr
 */

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iomanip> // setprecision
#include <list>
#include <unistd.h>
#include <ctime>
#include <boost/utility.hpp>
#include <math.h>
#include <algorithm>
#include <iterator>     /* std::next & std::prev */
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>

#define VERY_VERBOSE 0

int count_digits (const int &the_integer){
	int nDigits = std::floor(std::log10(std::abs(the_integer))) + 1;

	return nDigits;
}


std::string int2str_setw (const int &a,const int &w){
	std::string out;
	std::ostringstream out_ss;
	out_ss << std::setfill('0') << std::setw(w) << a;
	out = out_ss.str();

	return out;
}
//int main(int argc, char** argv) {

int main() {

	int num_frac;
	std::string filename1("../output/params.txt");
	std::string filename2("../output/poly_reorderedparentName.dat");
	std::string filename3("../output/insidebbox_frac.txt");
	std::string filename4("../desactFracsFile.dat");

	std::ifstream    infile1(filename1.c_str());
	std::string line, linenospace;
	if (! infile1.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename1 << "\"."<< std::endl;
		exit(-1);
	}

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile1, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_frac;
	}

	std::getline(infile1, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile1, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	std::list<int> fracNameList, fracParentNameList, newfracParentNameList, oldfracParentNameList, finalfracParentNameList;

	int num_count_frac = 0;
	while(num_count_frac < num_frac && !infile1.eof()){
		//std::cout << "Line : " << line << std::endl;
		int dummy_int;
		double dummy_double;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> dummy_int >> dummy_double;

		//if (dummy_double <= -89.9 && dummy_double >= -90.1){
		if (dummy_double <= -89 && dummy_double >= -91){
			fracNameList.push_back(dummy_int);
		}
		num_count_frac++;
		std::getline(infile1, line);
	}

	fracNameList.sort();
	fracNameList.unique();

	infile1.close();

	std::ifstream    infile2(filename2.c_str());
	if (! infile2.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename2 << "\"."<< std::endl;
		exit(-1);
	}
	while(!infile2.eof()){
		std::getline(infile2, line);
		int dummy_int;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> dummy_int;
		fracParentNameList.push_back(dummy_int);
	}

	for (std::list<int>::iterator it = fracNameList.begin(); it != fracNameList.end(); it++){
		std::cout << "FracName = " << *it << std::endl;
		std::list<int>::iterator fracParent_iter(fracParentNameList.begin());
		std::advance(fracParent_iter, *it-1);
		newfracParentNameList.push_back(*fracParent_iter);
	}
	newfracParentNameList.sort();
	newfracParentNameList.unique();

	infile2.close();

	std::ifstream    infile3(filename3.c_str());
	if (! infile3.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename3 << "\"."<< std::endl;
		exit(-1);
	}
	while(!infile3.eof()){
		std::getline(infile3, line);
		int dummy_int;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> dummy_int;
		oldfracParentNameList.push_back(dummy_int);
	}

	for (std::list<int>::iterator it = newfracParentNameList.begin(); it != newfracParentNameList.end(); it++){
		std::cout << "FracParentName = " << *it << std::endl;
		std::list<int>::iterator oldfracParent_iter(oldfracParentNameList.begin());
		std::advance(oldfracParent_iter, *it-1);
		finalfracParentNameList.push_back(*oldfracParent_iter);
	}
	finalfracParentNameList.sort();
	finalfracParentNameList.unique();

	infile3.close();

	std::ofstream    outfile(filename4.c_str());
	for (std::list<int>::iterator it = finalfracParentNameList.begin(); it != finalfracParentNameList.end(); it++){
		std::cout << "FinalFracParentName = " << *it << std::endl;
		outfile << *it << std::endl;
	}
	return 1;
}



