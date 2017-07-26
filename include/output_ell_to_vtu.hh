/*
 * output_ell_to_vtu.hh
 *
 *  Created on: Jan 16, 2017
 *      Author: ngotr
 */

#ifndef INCLUDE_OUTPUT_ELL_TO_VTU_HH_
#define INCLUDE_OUTPUT_ELL_TO_VTU_HH_

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>

#include "read_off_write_vtk_xml_file.h"

namespace CGAL{

template <typename K, typename NT>
void output_ell_to_vtu(std::string outputName, CGAL::Surface_mesh<CGAL::Point_3<Kernel> > &sm){

	typedef typename CGAL::Point_3<Kernel>                    	Point_3;
	typedef typename CGAL::Surface_mesh<Point_3>				SurfaceMesh_;


	typedef typename SurfaceMesh_::Vertex_index 		   vertex_descriptor;
	typedef typename SurfaceMesh_::Vertex_iterator vertex_iterator;

	// Write mesh to vtu file
	std::string line;
	int num_vertices(sm.number_of_vertices()), num_cells(1), dummy, vertices_per_cell(sm.number_of_vertices());
	std::vector< std::vector<NT> > vertices_(num_vertices, std::vector<NT>(3));

	int iter_vert(0);
	for (vertex_iterator it = sm.vertices_begin(); it != sm.vertices_end(); ++it){
		//std::cout << "Point of poly_1_child_2 : " << sm.point(*it) << std::endl;
		Point_3 vertex_i = sm.point(*it);
		vertices_[iter_vert][0] = vertex_i[0];
		vertices_[iter_vert][1] = vertex_i[1];
		vertices_[iter_vert][2] = vertex_i[2];
		iter_vert++;
	}

	std::string  outfile_(outputName + ".vtu");
	std::ofstream    outfile;
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
		outfile << std::setprecision(15) << std::uppercase << std::scientific << std::setw(20)
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
			outfile << i << " ";
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
}

}

#endif /* INCLUDE_OUTPUT_ELL_TO_VTU_HH_ */
