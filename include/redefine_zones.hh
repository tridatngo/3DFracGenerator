/*
 * redefine_zones.hh
 *
 *  Created on: Oct 25, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_REDEFINE_ZONES_HH_
#define INCLUDE_REDEFINE_ZONES_HH_

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
#include <sys/types.h>
#include <sys/stat.h>

namespace CGAL{

void redefine_zones(std::string &lagrit_path_, double &eps_bound){

	// Section 8 : redefine zones
	// Creates lagrit script to define domain size

	std::string outfile_("domainattr.lgi");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	outfile << "read / gmv / full_mesh.gmv / mo\n";
	outfile << "cmo / printatt / mo/ -xyz- / minmax\n";
	outfile << "finish\n";
	outfile.close();

	std::string cmd;
	std::string lagrit_path(lagrit_path_);

	cmd = lagrit_path + std::string(" < domainattr.lgi > printxyz.out");
	std::system(cmd.c_str());

	// Read data from lagrit output file
	std::string infile_("printxyz.out");
	std::ifstream    infile;
	infile.open(infile_.c_str());

	std::string line, linenospace;

	float x_min, x_max, y_min, y_max, z_min, z_max;

	//infile.open(inputFile);
	while (!infile.eof())
	{
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("xic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}
				/*
				std::cout << " first : " << first << std::endl;
				std::cout << " last : " << last << std::endl;
				 */

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> x_min >> x_max;
			}
		}

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("yic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> y_min >> y_max;
			}
		}


		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("zic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> z_min >> z_max;
			}
		}
	}
	infile.close();

	// lagrit scripts to create new zone files: boundary zones

	std::string outfile_1("bound_zones.lgi");
	std::ofstream    outfile1;
	outfile1.open(outfile_1.c_str());

	outfile1 << "read / gmv / full_mesh.gmv / mo\n";
	outfile1 << "define / XMAX / "<< std::setprecision(12) << x_max - eps_bound <<"\n";
	outfile1 << "define / XMIN / "<< std::setprecision(12) << x_min + eps_bound <<"\n";
	outfile1 << "define / YMAX / "<< std::setprecision(12) << y_max - eps_bound <<"\n";
	outfile1 << "define / YMIN / "<< std::setprecision(12) << y_min + eps_bound <<"\n";
	outfile1 << "define / ZMAX / "<< std::setprecision(12) << z_max - eps_bound <<"\n";
	outfile1 << "define / ZMIN / "<< std::setprecision(12) << z_min + eps_bound <<"\n";

	outfile1 << "\n";

	outfile1 << "pset / top/ attribute / zic / 1,0,0/ gt /ZMAX\n";
	outfile1 << "pset / bottom/ attribute/ zic/ 1,0,0/ lt/ZMIN\n";
	outfile1 << "pset / left_w / attribute/ xic/ 1,0,0 /lt / XMIN\n";

	outfile1 << "\n";

	outfile1 << "pset / front_s / attribute/ yic / 1,0,0 / gt/YMAX\n";
	outfile1 << "pset / right_e / attribute/ xic/1,0,0/ gt/XMAX\n";
	outfile1 << "pset / back_n / attribute/ yic/ 1,0,0 / lt/YMIN\n";
	outfile1 << "pset/-all-/ zone / boundary / ascii\n";
	outfile1 << "finish\n";
	outfile1.close();

	//f.write(lagrit_input%parameters)

	cmd = lagrit_path + " < bound_zones.lgi > boundary_output.txt ";
	std::system(cmd.c_str());

	std::system("cp boundary_bottom.zone pboundary_bottom.zone");
	std::system("cp boundary_left_w.zone pboundary_left_w.zone");
	std::system("cp boundary_front_s.zone pboundary_front_s.zone");
	std::system("cp boundary_right_e.zone pboundary_right_e.zone");
	std::system("cp boundary_back_n.zone pboundary_back_n.zone");
	std::system("cp boundary_top.zone pboundary_top.zone");

	/*
	 * sed '$d' mon_fichier.txt : Delete the last line of file
	 * sed '1d' mon_fichier.txt : Delete the first line of file
	 * sed '7,9d' mon_fichier.txt : Delete the lines between 7th and 9th lines
	 */

	for (int i = 0; i!=2; i++){
		std::system("sed -i '$d' boundary_top.zone ");
		std::system("sed -i '$d' boundary_bottom.zone ");
		std::system("sed -i '$d' boundary_left_w.zone ");
		std::system("sed -i '$d' boundary_front_s.zone ");
		std::system("sed -i '$d' boundary_right_e.zone ");
	}

	std::system("sed -i '1d' boundary_bottom.zone ");
	std::system("sed -i '1d' boundary_left_w.zone ");
	std::system("sed -i '1d' boundary_front_s.zone ");
	std::system("sed -i '1d' boundary_right_e.zone ");
	std::system("sed -i '1d' boundary_back_n.zone ");
	std::system("cat boundary_top.zone boundary_bottom.zone boundary_left_w.zone boundary_front_s.zone boundary_right_e.zone boundary_back_n.zone  > allboundaries.zone ");

}


void redefine_zones_dimScaling(const std::string &lagrit_path_, const double &eps_bound, const int &dimScaling){

	// Section 8 : redefine zones
	// Creates lagrit script to define domain size

	std::string outfile_("domainattr.lgi");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	outfile << "read / gmv / full_mesh.gmv / mo\n";
	outfile << "cmo / printatt / mo/ -xyz- / minmax\n";
	outfile << "finish\n";
	outfile.close();

	std::string cmd;
	std::string lagrit_path(lagrit_path_);

	cmd = lagrit_path + std::string(" < domainattr.lgi > printxyz.out");
	std::system(cmd.c_str());

	// Read data from lagrit output file
	std::string infile_("printxyz.out");
	std::ifstream    infile;
	infile.open(infile_.c_str());

	std::string line, linenospace;

	float x_min, x_max, y_min, y_max, z_min, z_max;

	//infile.open(inputFile);
	while (!infile.eof())
	{
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("xic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}
				/*
				std::cout << " first : " << first << std::endl;
				std::cout << " last : " << last << std::endl;
				 */

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> x_min >> x_max;
			}
		}

		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("yic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> y_min >> y_max;
			}
		}


		if (linenospace.substr (0,1) != "#" ){
			if (linenospace.find("zic") != std::string::npos){

				unsigned first = line.find("c");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> z_min >> z_max;
			}
		}
	}
	//std::cout <<"z_min = " << z_min << ", zmax = " << z_max<< ", dimScaling = "<< dimScaling << std::endl;

	infile.close();

	// lagrit scripts to create new zone files: boundary zones

	std::system("rm -rf tri_fracture.stor");

	std::string outfile_1("bound_zones.lgi");
	std::ofstream    outfile1;
	outfile1.open(outfile_1.c_str());

	//outfile1 << "read / gmv / full_mesh.gmv / mo\n";
	outfile1 << "read / avs / full_mesh_DFNWorks.inp / mo\n";
	outfile1 << "\n";
	outfile1 << "define / XMAX / "<< std::setprecision(12) << (x_max - eps_bound)*dimScaling <<"\n";
	outfile1 << "define / XMIN / "<< std::setprecision(12) << (x_min + eps_bound)*dimScaling <<"\n";
	outfile1 << "define / YMAX / "<< std::setprecision(12) << (y_max - eps_bound)*dimScaling <<"\n";
	outfile1 << "define / YMIN / "<< std::setprecision(12) << (y_min + eps_bound)*dimScaling <<"\n";
	outfile1 << "define / ZMAX / "<< std::setprecision(12) << (z_max - eps_bound)*dimScaling <<"\n";
	outfile1 << "define / ZMIN / "<< std::setprecision(12) << (z_min + eps_bound)*dimScaling <<"\n";

	outfile1 << "\n";

	outfile1 << "pset / top/ attribute / zic / 1,0,0/ gt /ZMAX\n";
	outfile1 << "pset / bottom/ attribute/ zic/ 1,0,0/ lt/ZMIN\n";
	outfile1 << "pset / left_w / attribute/ xic/ 1,0,0 /lt / XMIN\n";

	outfile1 << "\n";

	outfile1 << "pset / front_s / attribute/ yic / 1,0,0 / gt/YMAX\n";
	outfile1 << "pset / right_e / attribute/ xic/1,0,0/ gt/XMAX\n";
	outfile1 << "pset / back_n / attribute/ yic/ 1,0,0 / lt/YMIN\n";
	outfile1 << "pset/-all-/ zone / boundary / ascii\n";
	outfile1 << "\n";
	outfile1 << "finish\n";
	outfile1.close();

	//f.write(lagrit_input%parameters)

	cmd = lagrit_path + " < bound_zones.lgi > boundary_output.txt ";
	std::system(cmd.c_str());

	std::system("cp boundary_bottom.zone pboundary_bottom.zone");
	std::system("cp boundary_left_w.zone pboundary_left_w.zone");
	std::system("cp boundary_front_s.zone pboundary_front_s.zone");
	std::system("cp boundary_right_e.zone pboundary_right_e.zone");
	std::system("cp boundary_back_n.zone pboundary_back_n.zone");
	std::system("cp boundary_top.zone pboundary_top.zone");

	/*
	 * sed '$d' mon_fichier.txt : Delete the last line of file
	 * sed '1d' mon_fichier.txt : Delete the first line of file
	 * sed '7,9d' mon_fichier.txt : Delete the lines between 7th and 9th lines
	 */

	for (int i = 0; i!=2; i++){
		std::system("sed -i '$d' boundary_top.zone ");
		std::system("sed -i '$d' boundary_bottom.zone ");
		std::system("sed -i '$d' boundary_left_w.zone ");
		std::system("sed -i '$d' boundary_front_s.zone ");
		std::system("sed -i '$d' boundary_right_e.zone ");
	}

	std::system("sed -i '1d' boundary_bottom.zone ");
	std::system("sed -i '1d' boundary_left_w.zone ");
	std::system("sed -i '1d' boundary_front_s.zone ");
	std::system("sed -i '1d' boundary_right_e.zone ");
	std::system("sed -i '1d' boundary_back_n.zone ");
	std::system("cat boundary_top.zone boundary_bottom.zone boundary_left_w.zone boundary_front_s.zone boundary_right_e.zone boundary_back_n.zone  > allboundaries.zone ");

}

}

#endif /* INCLUDE_REDEFINE_ZONES_HH_ */
