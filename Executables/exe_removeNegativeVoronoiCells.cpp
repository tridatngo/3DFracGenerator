/*
 * exe_removeNegativeVoronoiCells.cpp
 *
 *  Created on: Jan 5, 2017
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
	std::string filename1("output/params.txt");
	std::string filename2("lagrit_logs/rmNegativeVoronoiCells.dat");


	double MinVoronoiArea;

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
	infile1.close();

	std::list<int> fracNameList;
	int num_count_frac = 0;


	std::cout << "Writing : merge_poly_new.lgi \n";
	std::string outfile_("merge_poly_new.lgi");
	std::ofstream    outfile;
	//outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);
	outfile.open(outfile_.c_str());

	int digits(count_digits(num_frac));

	for (int iter = 1; iter!=num_frac+1; iter++){
		std::string files("lagrit_logs/log_lagrit_" + int2str_setw(iter,digits));
		std::cout << "file : " <<  files << std::endl;

		std::ifstream    infile(files.c_str());
		if (! infile.is_open()) {
			std::cerr << "Failed to open the file "<< files << std::endl;
			exit(-1);
		}

		while (!infile.eof())
		{
			std::getline(infile, line);
			//std::cout << "line : " <<  line << std::endl;

			if (line.find("Minimum Voronoi area        =") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> MinVoronoiArea;
				std::cout << "MinVoronoiArea = " << std::scientific << std::setprecision(6) << MinVoronoiArea << std::endl;
				if (MinVoronoiArea < 0) {
					fracNameList.push_back(MinVoronoiArea);
				}
				else{

					outfile << "read / mesh_" << int2str_setw(iter,digits) << ".inp / mo_" << int2str_setw(iter,digits) << " \n";
					outfile << "define / MO_NAME_FRAC / mo_" << int2str_setw(iter,digits) << " \n";

					outfile << "cmo / addatt / MO_NAME_FRAC / volume / evol_onen \n";
					outfile << "math / sum / MO_NAME_FRAC / evol_sum / 1 0 0 / MO_NAME_FRAC / evol_one \n";

					outfile << "addmesh / merge / cmo_tmp / cmo_tmp / mo_" << int2str_setw(iter,digits) << " \n";
					outfile << "cmo / delete / mo_" << int2str_setw(iter,digits) << " \n";
					outfile << "\n";
				}
			}
		}
		infile.close();
	}


	outfile << " # Writing out merged fractures\n";

	outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
	outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
	outfile << "cmo select cmo_tmp\n";

	outfile << "dump lagrit part_1.lg cmo_tmp\n";

	outfile << "finish \n";

	std::cout << "fracNameList.size() = " << fracNameList.size() << std::endl;

	/*
	fracNameList.sort();
	fracNameList.unique();

	std::ofstream    outfile(filename2.c_str());
	for (std::list<int>::iterator it = fracNameList.begin(); it != fracNameList.end(); it++){
		std::cout << "fracNameList = " << *it << std::endl;
		outfile << *it << std::endl;
	}
	*/

	// Path to lagrit executable
	std::string lagrit_path( std::string("/work/irlin104_1/ngotr/Outils/lagrit/lagrit_ulin3.2") );

	std::string cmd;
	cmd = lagrit_path + std::string(" < merge_poly_new.lgi > lagrit_logs/log_merge_poly_new");
	std::system(cmd.c_str());

	return 1;
}


