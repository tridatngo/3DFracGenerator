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
	if( argc == 2 ) {
		inp_file = argv[1];
	}
	else {
		std::cout << "Usage: ./exefile InpInputFile \n";
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

	std::string  outfile_(inp_file.substr(0,inp_file.length()-4)+ ".off");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	// header
	outfile << "OFF\n";
	outfile << num_vertices << " " << num_cells << " 0\n";
	// Write vertices
	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	int i = 0;
	double eps = 1E-20;
	for (int it = 0; it != num_vertices; ++it)
	{	
		/*
		for (int j = 0; j != 3; j++) {
			if (abs(vertices_[it][j]) < eps) {
				vertices_[it][j] = 0;
			}
		}
		*/
		outfile << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1] << " " << vertices_[it][2] << std::endl;
	}

	// Write cells
	//outfile << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell + 1) << std::endl;

	if (num_cells == 1){
		outfile << vertices_per_cell << " ";
		for (i = 0; i < vertices_per_cell; i++){
			outfile << i << " ";
		}
	}
	else{
		for (int it = 0; it != num_cells; ++it)
		{
			outfile << vertices_per_cell << " ";
			for (i = 0; i < vertices_per_cell; i++){
				outfile << cells_[it][i] - 1 << " ";
			}
			outfile << std::endl;
		}
	}
	/*
	// cell types (type 5 is a three-node triangle)

	outfile << "\n";
	outfile << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3) {
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "5" << std::endl;
		}
	}
	else {
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "7" << std::endl;
		}
	}
	*/
	outfile.close();

	std::cout << "Writing *.off file: Complete." << std::endl;

	return 0;
}
