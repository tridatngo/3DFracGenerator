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

	std::cout<< "\nCreating parameter_i.mlgi files" << std::endl;
	std::system("rm -rf output/parameters");
	std::system("mkdir output/parameters");
	std::string filename_("output/params.txt");

	int E_name, num_ell, digits;
	double theta, dL0, x1, y1, z1, x2, y2, z2;

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	std::istringstream ist(line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell;
	}
	digits = count_digits(num_ell);

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	for (int it =0; it!=num_ell; it++){
		E_name = it;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> E_name >> theta >> dL0 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;

		std::string  	outfile_("output/parameters/parameters_" + int2str_setw(E_name,digits) + ".mlgi");
		std::ofstream   outfile;
		outfile.open(outfile_.c_str());

		// header
		outfile << "define / ID / " << E_name << "\n";
		outfile << "define / OUTFILE_GMV / mesh_" << int2str_setw(E_name,digits) <<  ".gmv\n";
		outfile << "define / OUTFILE_AVS / mesh_" <<  int2str_setw(E_name,digits) << ".inp\n";

		// NTD 03/11/2016 [IMP --> change fracture connectivity when using std::setprecision(12), to be improve]
		outfile << "define / THETA  / " << theta << "\n";
		//outfile << "define / THETA  / " << std::setprecision(12) << theta << "\n";

		outfile << "define / dL0 / " << std::setprecision(12) << dL0 << "\n";
		outfile << "define / X1 / " << std::setprecision(12) << x1 << "\n";
		outfile << "define / Y1 / " << std::setprecision(12) << y1 << "\n";
		outfile << "define / Z1 / " << std::setprecision(12) << z1 << "\n";
		outfile << "define / X2 / " << std::setprecision(12) << x2 << "\n";
		outfile << "define / Y2 / " << std::setprecision(12) << y2 << "\n";
		outfile << "define / Z2 / " << std::setprecision(12) << z2 << "\n";
		outfile << "finish \n";

		outfile.close();

		std::getline(infile, line);
	}

	infile.close();

	std::cout << "Creating parameter_i.mlgi files: Complete." << std::endl;

	return 1;
}



