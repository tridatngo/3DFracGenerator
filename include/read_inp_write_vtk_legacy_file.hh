/*
 * read_inp_write_vtk_legacy_file.hh
 *
 *  Created on: Oct 3, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_READ_INP_WRITE_VTK_LEGACY_FILE_HH_
#define INCLUDE_READ_INP_WRITE_VTK_LEGACY_FILE_HH_

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
#include "myglobal_functions.hh"

namespace CGAL {
//template<typename VtkMultiWriter>

bool read_inp_write_vtk_legacy_file(const std::string &filename_)
{

	const std::string inp_file (filename_ + ".inp");
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

		ist >> num_vertices >> num_cells >> dummy;
	}

	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));

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
		for (int i =0; i < vertices_per_cell; i++)
			ist >> cells_[icells][i];
	}

	infile.close();

	std::string  outfile_(filename_ + ".vtk");
	std::ofstream    outfile;outfile.open(outfile_.c_str());

	// header

	outfile << "# vtk DataFile Version 4.0\n";
	outfile << "vtk output\n";
	outfile << "ASCII\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";
	outfile << "POINTS " << num_vertices << " float\n";

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
	}

	// Write cells
	outfile << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfile << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] - 1 << " ";
		outfile << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfile << "\n";
	outfile << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "7" << std::endl;
		}
	}

	outfile.close();

	return true;
}

bool read_inp_write_vtk_legacy_file_with_matID(const std::string &filename_)
{

	const std::string inp_file (filename_ + ".inp");
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

	//std::cout  << "dummy_2 = " << dummy_2 << std::endl;

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

	std::string  outfile_(filename_ + ".vtk");
	std::ofstream  outfile; outfile.open(outfile_.c_str());

	// header

	outfile << "# vtk DataFile Version 4.0\n";
	outfile << "vtk output\n";
	outfile << "ASCII\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";
	outfile << "POINTS " << num_vertices << " float\n";

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
	}

	// Write cells
	outfile << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfile << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] - 1 << " ";
		outfile << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfile << "\n";
	outfile << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "7" << std::endl;
		}
	}

	// Write material ID
	outfile << "CELL_DATA     " << num_cells << std::endl;
	outfile << "SCALARS Material_ID int 1\n";
	outfile << "LOOKUP_TABLE	default\n";

	//outfile << "Material_ID " << num_vertices << "  int\n";

	for (std::list<int>::iterator it=material_ID.begin(); it != material_ID.end(); ++it)
	{
		outfile  << "  " << *it << std::endl;
	}

	outfile.close();

	return true;
}

bool read_inp_write_vtk_legacy_file_with_matID_dimScaling(const std::string &filename_, const int &dimScaling)
{

	const std::string inp_file (filename_ + ".inp");
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

	//std::cout  << "dummy_2 = " << dummy_2 << std::endl;

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

	std::string  outfile_(filename_ + ".vtk");
	std::ofstream  outfile; outfile.open(outfile_.c_str());

	// header

	outfile << "# vtk DataFile Version 4.0\n";
	outfile << "vtk output\n";
	outfile << "ASCII\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";
	outfile << "POINTS " << num_vertices << " float\n";

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
	}

	// Write cells
	outfile << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfile << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] - 1 << " ";
		outfile << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfile << "\n";
	outfile << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "7" << std::endl;
		}
	}

	// Write material ID
	outfile << "CELL_DATA     " << num_cells << std::endl;
	outfile << "SCALARS Material_ID int 1\n";
	outfile << "LOOKUP_TABLE	default\n";

	//outfile << "Material_ID " << num_vertices << "  int\n";

	for (std::list<int>::iterator it=material_ID.begin(); it != material_ID.end(); ++it)
	{
		outfile  << "  " << *it << std::endl;
	}

	outfile.close();

	// Outfile with dimScaling
	std::string  outfile_dim(filename_ + "_orig.vtk");
	std::ofstream  outfileDim; outfileDim.open(outfile_dim.c_str());

	// header
	outfileDim << "# vtk DataFile Version 4.0\n";
	outfileDim << "vtk output\n";
	outfileDim << "ASCII\n";
	outfileDim << "DATASET UNSTRUCTURED_GRID\n";
	outfileDim << "POINTS " << num_vertices << " float\n";

	// Write vertices
	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	for (int it=0; it != num_vertices; ++it)
	{
		for (int j = 0; j!=3; j++){
			if (abs( vertices_[it][j]) < eps){
				vertices_[it][j] = 0;
			}
		}

		outfileDim  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] * dimScaling << " " << vertices_[it][1] * dimScaling << " " << vertices_[it][2] * dimScaling << std::endl;
	}

	// Write cells
	outfileDim << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfileDim << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfileDim << cells_[it][i] - 1 << " ";
		outfileDim << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfileDim << "\n";
	outfileDim << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfileDim << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfileDim << "7" << std::endl;
		}
	}

	// Write material ID
	outfileDim << "CELL_DATA     " << num_cells << std::endl;
	outfileDim << "SCALARS Material_ID int 1\n";
	outfileDim << "LOOKUP_TABLE	default\n";

	//outfileDim << "Material_ID " << num_vertices << "  int\n";

	for (std::list<int>::iterator it=material_ID.begin(); it != material_ID.end(); ++it)
	{
		outfileDim  << "  " << *it << std::endl;
	}

	outfileDim.close();

	return true;
}


bool read_inp_write_vtk_legacy_file_with_matID_dimScaling_round(const std::string &filename_, const int &dimScaling, const int &intRound)
{

	const std::string inp_file (filename_ + ".inp");
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

	//std::cout  << "dummy_2 = " << dummy_2 << std::endl;

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

	std::string  outfile_(filename_ + ".vtk");
	std::ofstream  outfile; outfile.open(outfile_.c_str());

	// header
	outfile << "# vtk DataFile Version 4.0\n";
	outfile << "vtk output\n";
	outfile << "ASCII\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";
	outfile << "POINTS " << num_vertices << " float\n";

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

		outfile  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20)
				<< CGAL::pround(vertices_[it][0],intRound) << " "
				<< CGAL::pround(vertices_[it][1],intRound) << " "
				<< CGAL::pround(vertices_[it][2],intRound)  << std::endl;
	}

	// Write cells
	outfile << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfile << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] - 1 << " ";
		outfile << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfile << "\n";
	outfile << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << "7" << std::endl;
		}
	}

	// Write material ID
	outfile << "CELL_DATA     " << num_cells << std::endl;
	outfile << "SCALARS Material_ID int 1\n";
	outfile << "LOOKUP_TABLE	default\n";

	//outfile << "Material_ID " << num_vertices << "  int\n";

	for (std::list<int>::iterator it=material_ID.begin(); it != material_ID.end(); ++it)
	{
		outfile  << "  " << *it << std::endl;
	}

	outfile.close();

	// Outfile with dimScaling
	std::string  outfile_dim(filename_ + "_orig.vtk");
	std::ofstream  outfileDim; outfileDim.open(outfile_dim.c_str());

	// header
	outfileDim << "# vtk DataFile Version 4.0\n";
	outfileDim << "vtk output\n";
	outfileDim << "ASCII\n";
	outfileDim << "DATASET UNSTRUCTURED_GRID\n";
	outfileDim << "POINTS " << num_vertices << " float\n";

	// Write vertices
	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	for (int it=0; it != num_vertices; ++it)
	{
		for (int j = 0; j!=3; j++){
			if (abs( vertices_[it][j]) < eps){
				vertices_[it][j] = 0;
			}
		}

		outfileDim  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20)
			<< CGAL::pround(vertices_[it][0],intRound) * dimScaling << " "
			<< CGAL::pround(vertices_[it][1],intRound) * dimScaling << " "
			<< CGAL::pround(vertices_[it][2],intRound) * dimScaling << std::endl;
	}

	// Write cells
	outfileDim << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfileDim << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfileDim << cells_[it][i] - 1 << " ";
		outfileDim << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfileDim << "\n";
	outfileDim << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfileDim << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfileDim << "7" << std::endl;
		}
	}

	// Write material ID
	outfileDim << "CELL_DATA     " << num_cells << std::endl;
	outfileDim << "SCALARS Material_ID int 1\n";
	outfileDim << "LOOKUP_TABLE	default\n";

	//outfileDim << "Material_ID " << num_vertices << "  int\n";

	for (std::list<int>::iterator it=material_ID.begin(); it != material_ID.end(); ++it)
	{
		outfileDim  << "  " << *it << std::endl;
	}

	outfileDim.close();

	return true;
}
} // end namespace CGAL
#endif /* INCLUDE_READ_INP_WRITE_VTK_LEGACY_FILE_HH_ */
