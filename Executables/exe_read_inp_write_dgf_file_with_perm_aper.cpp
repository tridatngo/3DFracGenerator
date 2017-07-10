/*
 * exe_read_inp_write_dgf_file_with_perm_aper.cpp
 *
 *  Created on: Nov 03, 2016
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
#include <list>


int getnumline (const std::string filename){

	int numline(0);
	std::string line;
	std::ifstream    infile(filename.c_str());

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file \""<< filename << "\"."<< std::endl;
		exit(-1);
	}

	while (!infile.eof()){
		std::getline(infile, line);
		numline++;
	}
	return (numline - 1);
}

int main(int argc, char** argv) {
	
	bool runtime_status(true);
        std::string filename_ = "";
        if( argc == 2 ) {
                filename_ = argv[1];
        }
        else {
                std::cout << "Usage: ./exefile InpInputFile \n";
                std::cerr << "Failed to open the input file." << std::endl;
                return -1;
        }

	std::string apertureFile(std::string("aperture.dat")), permFile(std::string("permeability.dat"));

	// ===========================================================================
	// Read geometric data
	// ===========================================================================

	const std::string inp_file (filename_.substr(0,filename_.length()-4)+ ".inp");
	const char* filename = inp_file.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line, linenospace;
	int num_vertices, num_cells, dummy, vertices_per_cell;

	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.clear();

		ist >> num_vertices >> num_cells >> dummy;
	}

	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));
	std::list<int> material_ID;

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
		}
	}

	//std::cout << "Starting to read face data" << std::endl;

	std::getline(infile, line);
	std::istringstream ist(line);

	int dummy_1, dummy_2;
	std::string elem_type;

	ist >> dummy_1 >> dummy_2 >> elem_type;

	material_ID.push_back(dummy_2);

	if ( elem_type == std::string("tri")){
		vertices_per_cell = 3;
	}
	else if (elem_type == std::string("quad")) {
		vertices_per_cell = 4;
	}

	//std::cout << "Number of vertices per cell, vertices_per_cell = " << vertices_per_cell << std::endl;

	std::vector< std::vector<int> > cells_(num_cells, std::vector<int>(vertices_per_cell));

	for (int i =0; i < vertices_per_cell; i++)
		ist >> cells_[0][i];

	for (int icells=1; icells < num_cells; icells++)
	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist >> dummy_1 >> dummy_2 >> elem_type;

		material_ID.push_back(dummy_2);

		for (int i =0; i < vertices_per_cell; i++)
			ist >> cells_[icells][i];
	}

	infile.close();

	// ===========================================================================
	// Read parameter data
	// ===========================================================================

	std::list<double> permeability_X, permeability_Y, permeability_Z;
	std::ifstream    paramfile1(permFile.c_str());

	if (! paramfile1.is_open()) {
		std::cerr << "Error: Failed to open the input file." << std::endl;
		return -1;
	}

	int num_fracs;
	num_fracs = getnumline(permFile);
	for (int i=0; i!= num_fracs; i++){
		std::getline(paramfile1,line);
		std::istringstream ist(line);
		double dummy_X, dummy_Y, dummy_Z;
		ist >> dummy_X >> dummy_Y >> dummy_Z;
		permeability_X.push_back(dummy_X);
		permeability_Y.push_back(dummy_Y);
		permeability_Z.push_back(dummy_Z);
	}

	std::list<double> aperture;
	std::ifstream    paramfile2(apertureFile.c_str());

	if (! paramfile2.is_open()) {
		std::cerr << "Error: Failed to open the input file." << std::endl;
		return -1;
	}

	for (int i=0; i!= num_fracs; i++){
		std::getline(paramfile2,line);
		std::istringstream ist(line);
		double dummy;
		ist >> dummy;
		aperture.push_back(dummy);
	}

	// ===========================================================================
	// Write dgf file
	// ===========================================================================

	std::string  outfile_(filename_.substr(0,filename_.length()-4)+ ".dgf");
	std::ofstream  outfile; outfile.open(outfile_.c_str());

	// header
	outfile << "DGF\n";
	outfile << "% Elements = " << num_vertices << "  |  Vertices = " << num_cells << "\n";
	outfile << "\n";
	outfile << "VERTEX\n";

	// Write vertices
	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	int i=0;
	double eps = 1E-20;
	for (int it=0; it != num_vertices; ++it)
	{
		for (int j = 0; j!=3; j++){
			if (abs( vertices_[it][j]) < eps){
				vertices_[it][j] = 0;
			}
		}

		outfile  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1]  << " " << vertices_[it][2]  << std::endl;
		//outfile  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1]  << std::endl;
	}

	// Write cells
	outfile << "#\n";
	outfile << "\n";
	outfile << "SIMPLEX\n";
	outfile << "parameters 4\n";

	for (int it=0; it != num_cells; ++it)
	{
		std::list<int>::iterator mat_ID(material_ID.begin());
		std::advance(mat_ID, it);

		//outfile << vertices_per_cell  << " ";
		std::list<double>::iterator perm_i_X(permeability_X.begin());
		std::advance(perm_i_X, *mat_ID-1);

		std::list<double>::iterator perm_i_Y(permeability_Y.begin());
		std::advance(perm_i_Y, *mat_ID-1);

		std::list<double>::iterator perm_i_Z(permeability_Z.begin());
		std::advance(perm_i_Z, *mat_ID-1);

		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] - 1 << " ";
		outfile << *perm_i_X << " " << *perm_i_Y << " " << *perm_i_Z << " ";

		std::list<double>::iterator apert_i(aperture.begin());
		std::advance(apert_i, *mat_ID-1);
		outfile << *apert_i << std::endl;
	}

	outfile << "#\n";
	/*
	outfile << "\n";
	outfile << "BOUNDARYSEGMENTS\n";


	const std::string twin_file (filename_ + ".dgf_");
	std::ifstream    twinfile(twin_file.c_str());

	if (! twinfile.is_open()) {
		std::cerr << "Error: Failed to open the input file *.dgf_." << std::endl;
		return -1;
	}

	std::getline(twinfile, line);
	linenospace = line;
	//linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
	//linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

	while (linenospace.find("BOUNDARYSEGMENTS") == std::string::npos){
		std::getline(twinfile, line);
		linenospace = line;
		//linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		//linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
	}

	while (!twinfile.eof()){
		std::getline(twinfile, line);
		linenospace = line;
		//linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		//linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
		if (linenospace.find("#") != std::string::npos ){
			break;
		}
		std::istringstream ist(line);
		ist >> dummy >> dummy_1 >> dummy_2;
		outfile << "1   " << dummy_1 << " " << dummy_2 << "\n";

	}

	outfile << "\n";
	outfile << "#\n";
	*/
	outfile.close();

	return 0;
}
