/*
 * check_quality_voronoi_mesh.hh
 *
 *  Created on: Oct 6, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CHECK_QUALITY_VORONOI_MESH_HH_
#define INCLUDE_CHECK_QUALITY_VORONOI_MESH_HH_

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

namespace CGAL{
bool check_quality_voronoi_mesh() {

	std::ifstream    infile("lagrit_logs/log_merge_all");
	if (! infile.is_open()) {
		std::cerr << "Failed to open the lagrit_logs/log_merge_all file." << std::endl;
		exit(-1);
	}
	std::string line, linenospace;
	double min_aspect_ratio, min_volume;

	while (!infile.eof())
	{
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		//std::cout << "line : " << line << std::endl;

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("minaspectratio") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("max") != std::string::npos){
					last = line.find("max")-1;
				}

				//std::cout << "linenospace = " << linenospace << std::endl;
				//std::cout << "first = " << first << std::endl;
				//std::cout << "last = " << last << std::endl;


				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> min_aspect_ratio;
				std::cout << "min_aspect_ratio = " << std::scientific << std::setprecision(6) << min_aspect_ratio << std::endl;
			}


			if (linenospace.find("minvolume") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("max") != std::string::npos){
					last = line.find("max")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> min_volume;
				std::cout << "min_volume = " << std::scientific << std::setprecision(6) << min_volume << std::endl;
			}
		}
	}

	if (min_volume > 0 && min_aspect_ratio > 0){
		return true;
	}
	else{
		return false;
	}
}

}
#endif /* INCLUDE_CHECK_QUALITY_VORONOI_MESH_HH_ */
