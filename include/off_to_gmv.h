//
//  off_to_gmv.h
//
//  Created by ngotr on 5/1/13.
//

#ifndef off_to_gmv_h
#define off_to_gmv_h

#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<cstring>
#include <math.h>

// Get current time and date
#include <ctime>
//#include <chrono>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

namespace CGAL {
//template<typename VtkMultiWriter>

bool off_to_gmv(const std::string &filename_)
{
	// Open the input file and read the two polygons from it.
	std::string off_file = filename_ + ".off";
	const char* filename = off_file.c_str();
	std::ifstream    infile(filename);
	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	int num_vertices, num_cells, dummy, vertices_per_cell;

	std::getline(infile, line);

	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.clear();

		//std::cout << line << std::endl;
		ist >> num_vertices >> num_cells >> dummy;
	}
	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));
	//std::vector< std::vector<int> > cells_(num_cells, std::vector<int>(3));


	std::getline(infile, line);
	if (line.size() == 0)
	{
		//std::cout << "Skip this line " << std::endl;
		//std::cout << "Starting to read vertice data" << std::endl;
		for (int ivert=0; ivert!=num_vertices; ivert++)
		{
			std::getline(infile, line);
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> vertices_[ivert][0] >> vertices_[ivert][1] >> vertices_[ivert][2];
		}
	}
	else{
		//std::cout << "Starting to read vertice data" << std::endl;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> vertices_[0][0] >> vertices_[0][1] >> vertices_[0][2];

		for (int ivert=1; ivert!=num_vertices; ivert++)
		{
			std::getline(infile, line);
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> vertices_[ivert][0] >> vertices_[ivert][1] >> vertices_[ivert][2];
		}
	}

	//std::cout << "Starting to read face data" << std::endl;

	std::getline(infile, line);
	std::istringstream ist(line);
	ist >> vertices_per_cell;
	//std::cout << "Number of vertices per cell, vertices_per_cell = " << vertices_per_cell << std::endl;

	std::vector< std::vector<int> > cells_(num_cells, std::vector<int>(vertices_per_cell));

	for (int i =0; i < vertices_per_cell; i++)
		ist >> cells_[0][i];

	for (int icells=1; icells < num_cells; icells++)
	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist >> dummy;
		for (int i =0; i < vertices_per_cell; i++)
			ist >> cells_[icells][i];
	}

	// Get current time and date
	time_t tt;
	//tt = std::chrono::system_clock::to_time_t ( today );

	infile.close();

	//std::string  outfile_(filename_.substr(0,filename_.length()-4)+".gmv");
	std::string  outfile_(filename_+".gmv");
	//std::ofstream    outfile;
	std::fstream    outfile;
	FILE * fp;

	// File pointer
	fp = fopen (outfile_.c_str(), "w");

	//outfile.open(outfile_.c_str());

	int indent_size = 2;
	std::string indent_unit(indent_size, ' ');
	std::string indent = indent_unit;

	// header
	//outfile << "simdate " << ctime(&tt) << std::endl;
	fprintf (fp,"gmvinput ascii \n");
	fprintf (fp, "codename CGAL \n");
	fprintf (fp, "nodes %s %i \n", indent.c_str() , num_vertices);
	//fprintf (fp, "\n");

	// Write vertices
	outfile.setf(std::ios::scientific);
	int num_output_col(10);
	for (int i=0; i != 3; ++i){
		for (int it=0; it != num_vertices; ++it){
			if ( it % num_output_col == 0 && it >0) {
				fprintf (fp, "\n");
				char buffer [13];
				sprintf (buffer, "%11.6E", vertices_[it][i]);
				fprintf (fp, " %13.13s", buffer);
			}
			else{
				char buffer [13];
				sprintf (buffer, "%11.6E", vertices_[it][i]);
				fprintf (fp, " %13.13s", buffer);

			}
			if ( it == num_vertices-1){
				fprintf (fp, "\n");
			}
		}
	}



	/*
	for (int it=0; it != num_vertices; ++it)
	{
		outfile << indent;
		outfile << std::scientific << vertices_[it][0] << " " << vertices_[it][1]  << " " << vertices_[it][2]  << std::endl;
	}
	 */

	// Write cells
	fprintf (fp, "cells %s %i \n", indent.c_str() , num_cells);

	for (int it=0; it != num_cells; ++it)
	{
		fprintf (fp, "tri %s %i", indent.c_str(), 3);
		for (int i=0; i < vertices_per_cell; i++)
			// An iterator starts by 0 in C++ but by 1 in Fortran
			// Must adapt cell numbering in gmv output
			fprintf (fp, " %i %s", cells_[it][i] + 1, " ");

		fprintf (fp, "\n");
	}

	fprintf (fp, "endgmv \n");

	fclose(fp);
	return true;
}
} // end namespace CGAL

#endif

