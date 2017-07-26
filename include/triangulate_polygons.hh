/*
 * triangulate_polygons.hh
 *
 *  Created on: Sep 13, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_TRIANGULATE_POLYGONS_HH_
#define INCLUDE_TRIANGULATE_POLYGONS_HH_

#include <string.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdlib.h>		/* system, NULL, EXIT_FAILURE */
#include "myglobal_functions.hh"
#include "check_polygon_triangulation.hh"
#include <sys/types.h>
#include <sys/stat.h>

struct stat info_2;

namespace CGAL{

void modify_lgi_file(std::string &input){

	std::ifstream    infile(input.c_str());
	std::string line, linenospace;
	int num_line(0);

	while (!infile.eof())
	{
		std::getline(infile, line);
		num_line++;
	}
	infile.close();

	std::string 	output = input + "_mod";
	std::ofstream 	outfile(output.c_str());

	infile.open(input.c_str());
	for (int i = 0; i!=num_line; i++){
		std::getline(infile, line);
		if (i < 26  || i > 51){
			outfile << line << std::endl;
		}
	}
	outfile.close();
}

void triangulate_polygons (std::string &lagrit_path_, const int &nPoly, const bool correctMesh){
	/*
	Triangulate each polygon. This python script runs lagrit to mesh each individual fractures
	Line Format std::system( "lagrit < mesh_polyi.lgi ")
	 */

	std::cout << "Triangulating polygons..." << std::endl;
	int digits(CGAL::count_digits(nPoly));
	std::string cmd;

	// Path to lagrit executable
	std::string lagrit_path(lagrit_path_);

	for (int ell_num = 1; ell_num!=nPoly+1; ell_num++){

		cmd = std::string("ln -s output/parameters/parameters_") + CGAL::int2str_setw(ell_num,digits) + std::string(".mlgi")
		+ std::string(" parameters_") + CGAL::int2str_setw(ell_num,digits) + std::string(".mlgi");
		std::system (cmd.c_str());

		cmd = std::string("ln -s intersections/intersections_") + CGAL::int2str_setw(ell_num,digits) + std::string(".inp")
		+ std::string(" intersections_") + CGAL::int2str_setw(ell_num,digits) + std::string(".inp");
		std::system (cmd.c_str());

		cmd = std::string("ln -s output/lagrit/ellipse_cutoff_") +CGAL::int2str_setw(ell_num,digits) + std::string(".inp")
		+ std::string(" ellipse_cutoff_") + CGAL::int2str_setw(ell_num,digits)  + std::string(".inp");
		std::system (cmd.c_str());

		cmd = std::string("ln -s output/lagrit/mesh_poly_") +CGAL::int2str_setw(ell_num,digits) + std::string(".lgi")
		+ std::string(" mesh_poly_") + CGAL::int2str_setw(ell_num,digits)  + std::string(".lgi");

		std::system (cmd.c_str());

		cmd = lagrit_path + std::string(" < mesh_poly_" + CGAL::int2str_setw(ell_num,digits) + ".lgi" +
				" > lagrit_logs/log_lagrit_" + CGAL::int2str_setw(ell_num,digits) );
		std::system(cmd.c_str());

		std::string checkfile = std::string("mesh_" + CGAL::int2str_setw(ell_num,digits)) + ".inp";
		std::string inputFile = std::string("mesh_poly_" + CGAL::int2str_setw(ell_num,digits)) + ".lgi";
		
		//if (!boost::filesystem::exists(checkfile) && !CGAL::check_polygon_triangulation(ell_num,digits) ){

		/* Correct the mesh by adding tiny triangulated bits */
		if ( stat( checkfile.c_str(), &info_2) != 0 && !CGAL::check_polygon_triangulation(ell_num,digits) && correctMesh){
			//std::cout << "\n";
			std::cout << "==> Correcting the mesh of the close outline number " << CGAL::int2str_setw(ell_num,digits);
			CGAL::modify_lgi_file(inputFile);
			cmd = lagrit_path + std::string(" < " + inputFile + "_mod > lagrit_logs/log_lagrit_" + CGAL::int2str_setw(ell_num,digits) );
			std::system(cmd.c_str());
			std::cout << ". Done." << std::endl;
			std::cout << "\n";
		}
	}
}

}

#endif /* INCLUDE_TRIANGULATE_POLYGONS_HH_ */
