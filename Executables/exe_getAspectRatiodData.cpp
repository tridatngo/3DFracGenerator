/*
 * exe_getAspectRatioData.cpp
 *
 *  Created on: Dec 7, 2016
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
#include<cstring>
#include <iomanip> // setprecision
#include <math.h>
#include <cmath>
#include "../include/myglobal_functions.hh"

int main() {

	std::string fname = std::string("lagrit_logs/log_merge_all");
	std::ifstream    infile(fname.c_str());
	std::string line;
	int num_psmooth, a1, a2, a3, a4, a5, a6, a7;

	while (!infile.eof())
	{
		std::getline(infile, line);
		//std::cout << line << std::endl;

		if (line.substr (0,1) != "#" ){
			if (line.find("elements with aspect ratio < .01:") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a1;
			}
			else if (line.find("elements with aspect ratio b/w .01 and .02:") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a2;
			}
			else if (line.find("elements with aspect ratio b/w .02 and .05:") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a3;
			}
			else if (line.find("elements with aspect ratio b/w .05 and .1 :") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a4;
			}
			else if (line.find("elements with aspect ratio b/w .1  and .2 :") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a5;
			}
			else if (line.find("elements with aspect ratio b/w .2  and .5 :") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a6;
			}
			else if (line.find("elements with aspect ratio b/w .5  and 1. :") != std::string::npos){
				unsigned first = line.find(":");
				unsigned last = line.size();

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> a7;
			}
		}
	}
	std::string ofname = std::string("output/MeshAspectRatio.txt");
	std::ofstream    outfile(ofname.c_str());
	//std::cout << a1 << ", " << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ", " << a7 << std::endl;
	outfile << "# a1<0.01 -- 0.01<a2<0.02 -- 0.02<a3<0.05 -- 0.05<a4<0.1 -- 0.1<a5<0.2 -- 0.2<a5<0.5 -- 0.5<a5<1.0 -- all" << std::endl;
	outfile << a1 << std::endl;
	outfile << a2 << std::endl;
	outfile << a3 << std::endl;
	outfile << a4 << std::endl;
	outfile << a5 << std::endl;
	outfile << a6 << std::endl;
	outfile << a7 << std::endl;
	outfile << a1+a2+a3+a4+a5+a6+a7 << std::endl;

	return 0;
}

