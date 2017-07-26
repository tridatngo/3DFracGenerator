/*
 * check_polygon_triangulation.hh
 *
 *  Created on: Oct 7, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CHECK_POLYGON_TRIANGULATION_HH_
#define INCLUDE_CHECK_POLYGON_TRIANGULATION_HH_
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
#include <iomanip> // setprecision
#include "myglobal_functions.hh"

namespace CGAL{

bool check_polygon_triangulation(const int &ell_num,const int &digits) {

	std::string fname = std::string("lagrit_logs/log_lagrit_" + CGAL::int2str_setw(ell_num,digits));
	std::ifstream    infile(fname.c_str());
	std::string line, linenospace;
	int num_psmooth;

	while (!infile.eof())
	{
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		//std::cout << "line : " << line << std::endl;

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("THEPSETpsmoothHAS") != std::string::npos){
				unsigned first = line.find("HAS")+2;
				unsigned last = line.size();
				if (linenospace.find("POINTS") != std::string::npos){
					last = line.find("POINTS")-1;
				}

				//std::cout << "linenospace = " << linenospace << std::endl;
				//std::cout << "first = " << first << std::endl;
				//std::cout << "last = " << last << std::endl;

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> num_psmooth;
				//std::cout << "num_psmooth = " << num_psmooth << std::endl;

				break;
			}
		}
	}

	if (num_psmooth > 0){
		return true;
	}
	else{
		return false;
	}
}
}

#endif /* INCLUDE_CHECK_POLYGON_TRIANGULATION_HH_ */
