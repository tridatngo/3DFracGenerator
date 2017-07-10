/*
 * exe_read_inp_write_vtk_legacy_file.cpp
 *
 *  Created on: Sep 14, 2016
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

int main(int argc, char** argv) {

	std::string inp_file = "";
	std::string bb_inp_file = "";
	if( argc == 3 ) {
		inp_file = argv[1];
		bb_inp_file =argv[2];
	}
	else {
		std::cout << "Usage: ./exefile InpInputFile InpBBInputFile\n";
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	// Open the input file and read the two polygons from it.
	//std::string filename_("out_sm_ellipse.off");
	//const char* filename = (argc > 1) ? argv[1] : filename_.c_str();

	//const std::string inp_file (filename_ + ".inp");
	const char* filename = inp_file.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	int num_vertices, num_cells, dummy, vertices_per_cell;

	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.clear();

		//std::cout << line << std::endl;
		ist >> num_vertices >> num_cells >> dummy;
	}

	if (num_cells == 0){
		num_cells = 1;
		vertices_per_cell = num_vertices;
	}

	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));
	//std::vector< std::vector<int> > cells_(num_cells, std::vector<int>(3));

	//double minX, maxX, minY, maxY, minZ, maxZ;	

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
			ist >> dummy >> vertices_[ivert][0] >> vertices_[ivert][1] >> vertices_[ivert][2];
		}
	}
	else{
		//std::cout << "Starting to read vertice data" << std::endl;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> dummy >> vertices_[0][0] >> vertices_[0][1] >> vertices_[0][2];

		for (int ivert=1; ivert!=num_vertices; ivert++)
		{
			std::getline(infile, line);
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> dummy >> vertices_[ivert][0] >> vertices_[ivert][1] >> vertices_[ivert][2];
			/*std::cout << "Line " << line << std::endl;
			std::cout << "vertices_[ivert][0] = " << vertices_[ivert][0] << std::endl;
			*/
		}
	}

	//std::cout << "Starting to read face data" << std::endl;

	std::getline(infile, line);
	std::istringstream ist(line);

	int dummy_1, dummy_2;
	std::string elem_type;

	ist >> dummy_1 >> dummy_2 >> elem_type;

	if ( elem_type == std::string("tri")){
		vertices_per_cell = 3;
	}
	else if (elem_type == std::string("quad")) {
		vertices_per_cell = 4;
	}

	std::cout << "Number of vertices, num_vertices  = " << num_vertices << std::endl;
	std::cout << "Number of cells, num_cells  = " << num_cells << std::endl;
	std::cout << "Number of vertices per cell, vertices_per_cell = " << vertices_per_cell << std::endl;

	std::vector< std::vector<int> > cells_(num_cells, std::vector<int>(vertices_per_cell));

	for (int i =0; i < vertices_per_cell; i++)
		ist >> cells_[0][i];

	for (int icells=1; icells < num_cells; icells++)
	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist >> dummy_1 >> dummy_2 >> elem_type;
		for (int i =0; i < vertices_per_cell; i++)
			ist >> cells_[icells][i];
	}

	infile.close();

	std::ifstream    bb_infile(bb_inp_file.c_str());

	if (! bb_infile.is_open()) {
		std::cerr << "Failed to open the bb input file." << std::endl;
		return -1;
	}

	int bb_num_vertices, bb_num_cells, bb_vertices_per_cell;

	{
		std::getline(bb_infile, line);
		std::istringstream bb_ist(line);
		bb_ist.clear();

		//std::cout << line << std::endl;
		bb_ist >> bb_num_vertices >> bb_num_cells >> dummy;
	}

	if (bb_num_cells == 0){
		bb_num_cells = 1;
		bb_vertices_per_cell = bb_num_vertices;
	}

	std::vector< std::vector<double> > bb_vertices_(bb_num_vertices, std::vector<double>(3));

	std::getline(bb_infile, line);
	if (line.size() == 0)
	{
		for (int ivert=0; ivert!=bb_num_vertices; ivert++)
		{
			std::getline(bb_infile, line);
			std::istringstream bb_ist(line);
			bb_ist.str(line);
			bb_ist.clear();
			bb_ist >> dummy >> bb_vertices_[ivert][0] >> bb_vertices_[ivert][1] >> bb_vertices_[ivert][2];
		}
	}
	else{
		//std::cout << "Starting to read vertice data" << std::endl;
		std::istringstream bb_ist(line);
		bb_ist.str(line);
		bb_ist.clear();
		bb_ist >> dummy >> bb_vertices_[0][0] >> bb_vertices_[0][1] >> bb_vertices_[0][2];

		for (int ivert=1; ivert!=bb_num_vertices; ivert++)
		{
			std::getline(bb_infile, line);
			std::istringstream bb_ist(line);
			bb_ist.str(line);
			bb_ist.clear();
			bb_ist >> dummy >> bb_vertices_[ivert][0] >> bb_vertices_[ivert][1] >> bb_vertices_[ivert][2];
			/*std::cout << "Line " << line << std::endl;
			std::cout << "vertices_[ivert][0] = " << vertices_[ivert][0] << std::endl;
			*/
		}
	}

	//std::cout << "Starting to read face data" << std::endl;

	std::getline(bb_infile, line);
	std::istringstream bb_ist(line);

	std::string bb_elem_type;

	bb_ist >> dummy_1 >> dummy_2 >> bb_elem_type;

	if ( bb_elem_type == std::string("tri")){
		bb_vertices_per_cell = 3;
	}
	else if (bb_elem_type == std::string("quad")) {
		bb_vertices_per_cell = 4;
	}

	std::cout << "BB - Number of vertices, bb_num_vertices  = " << bb_num_vertices << std::endl;
	std::cout << "BB - Number of cells, num_cells  = " << bb_num_cells << std::endl;
	std::cout << "BB - Number of vertices per cell, bb_vertices_per_cell = " << bb_vertices_per_cell << std::endl;

	std::vector< std::vector<int> > bb_cells_(bb_num_cells, std::vector<int>(bb_vertices_per_cell));

	for (int i =0; i < bb_vertices_per_cell; i++)
		bb_ist >> bb_cells_[0][i];

	for (int icells=1; icells < bb_num_cells; icells++)
	{
		std::getline(bb_infile, line);
		std::istringstream bb_ist(line);
		bb_ist >> dummy_1 >> dummy_2 >> bb_elem_type;
		for (int i =0; i < bb_vertices_per_cell; i++)
			bb_ist >> bb_cells_[icells][i];
	}

	bb_infile.close();

	std::string  outfile_(inp_file.substr(0,inp_file.length()-4)+ ".poly");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	// header
	outfile << "# Part 1 - node list\n";
  	outfile << "# node count, 3 dim, no attribute, no boundary marker\n";
	outfile << num_vertices + bb_num_vertices << " 3 0 0\n";
	outfile << "# Node index, node coordinates\n"; 

	// Write vertices

	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	int i = 0;
	double eps = 1E-20;
	for (int it = 0; it != num_vertices; ++it)
	{
		outfile << std::setprecision(12) << std::uppercase << std::scientific << it+1 << " " << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1] << " " << vertices_[it][2] << std::endl;
	}

	for (int it = 0; it != bb_num_vertices; ++it)
	{
		outfile << std::setprecision(12) << std::uppercase << std::scientific << it+num_vertices+1 << " " << std::setw(20) << 
				bb_vertices_[it][0] << " " << bb_vertices_[it][1] << " " << bb_vertices_[it][2] << std::endl;
	}


	outfile << "\n";
 	outfile << "# Part 2 - facet list\n";
   	outfile << "# facet count, no boundary marker\n";
   	outfile << num_cells + bb_num_cells << "  0\n";
   	outfile << "# facets\n";

        for (int it = 0; it != num_cells; ++it)
        {
		outfile << "1\n";
                outfile << vertices_per_cell << " ";
                for (i = 0; i < vertices_per_cell; i++){
                        outfile << cells_[it][i] << " ";
                }
                outfile << std::endl;
        }
	
       for (int it = 0; it != bb_num_cells; ++it)
        {
		outfile << "1\n";
                outfile << bb_vertices_per_cell << " ";
                for (i = 0; i < bb_vertices_per_cell; i++){
                        outfile << bb_cells_[it][i] + num_vertices << " ";
                }
                outfile << std::endl;
        }

	outfile << "\n";
 	outfile << "# Part 3 - hole list\n";
	outfile << "0            # no hole\n";
	
	outfile << "\n";
 	outfile << "# Part 4 - region list\n";
	outfile << "0            # no region\n";
	outfile.close();

	std::cout << "Writing *.poly file: Complete." << std::endl;

	return 0;
}
