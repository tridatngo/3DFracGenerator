//
//  write_mesh_to_vtk_xml_file.h
//
//  Created by David Bernstein on 5/1/13.
//

#ifndef read_off_write_vtk_xml_file_h
#define read_off_write_vtk_xml_file_h

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

namespace CGAL {
//template<typename VtkMultiWriter>

bool read_off_write_vtk_xml_file(const std::string &filename_)
{
	// Open the input file and read the two polygons from it.
	//std::string filename_("out_sm_ellipse.off");
	//const char* filename = (argc > 1) ? argv[1] : filename_.c_str();

	const std::string off_file (filename_ + ".off");
	const char* filename = off_file.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}
	//std::cout << "Length of char* " << sizeof(filename.str()) << std::endl;

	//std::cout << "File name is " << filename.substr(0,filename.length()-4)+".vtu" << std::endl;

	//std::ifstream infile(file_name);
	//std::cout << "Read permeability data from " << concd << std::endl;

	std::string line;
	int num_vertices, num_cells, dummy, vertices_per_cell;

	std::getline(infile, line);
	//std::cout << line << std::endl;
	//std::cout << "" << std::endl;

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

	infile.close();

	//std::string  outfile_(filename_.substr(0,filename_.length()-4)+".vtu");
	std::string  outfile_(filename_ + ".vtu");
	std::ofstream    outfile;
	//outfile.open(filename.substr(0,file_name.length()-4)+".vtu", std::ios_base::app);
	//outfile.open("tridatngo.vtu", std::ios_base::app);
	//outfile.open("tridatngo.vtu");
	outfile.open(outfile_.c_str());

	// header
	outfile << "<VTKFile type=\"UnstructuredGrid\" ";
	outfile << "version=\"0.1\" ";
	outfile << "byte_order=\"BigEndian\">" << std::endl;

	int indent_size = 2;
	std::string indent_unit(indent_size, ' ');
	std::string indent = indent_unit;
	outfile << indent + "<UnstructuredGrid>" << std::endl;

	// write mesh
	indent += indent_unit;
	outfile << indent + "<Piece NumberOfPoints=\"" << num_vertices << "\" ";
	outfile << "NumberOfCells=\"" << num_cells << "\">" << std::endl;

	// Write vertices
	indent += indent_unit;
	outfile << indent + "<Points>" << std::endl;

	indent += indent_unit;
	outfile << indent;
	outfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;

	int i=0;
	indent += indent_unit;

	for (int it=0; it != num_vertices; ++it)
	{
		outfile << indent;
		outfile << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20)
				<< vertices_[it][0] << " " << vertices_[it][1]  << " " << vertices_[it][2]  << std::endl;
	}

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</DataArray>" << std::endl;

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</Points>" << std::endl;

	// Write cells
	outfile << indent << "<Cells>" << std::endl;

	indent += indent_unit;
	outfile << indent;
	outfile << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">";
	outfile << std::endl;

	indent += indent_unit;
	for (int it=0; it != num_cells; ++it)
	{
		outfile << indent;
		for (i=0; i < vertices_per_cell; i++)
			outfile << cells_[it][i] << " ";
		outfile << std::endl;
	}

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</DataArray>" << std::endl;

	// offsets
	// every element is a three node triangle so all offsets are multiples of 3
	outfile << indent;
	outfile << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
	outfile << std::endl;
	//i = 3;

	indent += indent_unit;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << indent << i << std::endl;
			i += 3;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << indent << vertices_per_cell << std::endl;
			i += 3;
		}
	}

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</DataArray>" << std::endl;

	// cell types (type 5 is a three node triangle)

	outfile << indent;
	outfile << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">";
	outfile << std::endl;
	indent += indent_unit;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << indent << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile << indent << "7" << std::endl;
		}
	}

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</DataArray>" << std::endl;

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</Cells>" << std::endl;

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</Piece>" << std::endl;

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << indent + "</UnstructuredGrid>" << std::endl;

	indent.erase(indent.length()-indent_size, indent_size);
	outfile << "</VTKFile>" << std::endl;

	outfile.close();

	return true;
}
} // end namespace CGAL

#endif

