/*
 * write_full_mesh_for_DFNWorks.hh
 *
 *  Created on: Oct 26, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_WRITE_FULL_MESH_FOR_DFNWORKS_HH_
#define INCLUDE_WRITE_FULL_MESH_FOR_DFNWORKS_HH_

/*
 * Read full_mesh.inp file and write full_mesh_DFNWorks.inp to adapt for DFNWorks simulation.
 * This function update the parent name of the fracture in the *.inp file.
 */
#include <list>
#include "myglobal_functions.hh"
#include <iostream>
#include <sstream>
#include <stdio.h>      /* printf */
#include <math.h>       /* round, floor, ceil, trunc */
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>

namespace CGAL{

int return_index_inf(const std::list<int> &intList, const int a){
	// return the index of the maximal element of the integer list  which is inferior than a

	int index(0);

	if (*intList.begin() > a){
		return index;
	}
	else{
		unsigned int inum(0);
		for (; inum < intList.size(); inum++){
			std::list<int>::const_iterator int_iter = intList.begin();
			std::advance(int_iter,inum);
			//std::cout << "*int_iter = " << *int_iter << std::endl;

			if(*int_iter > a){
				index = inum;
				break;
			}
		}

		if (*boost::prior(intList.end()) < a){
			index = intList.size();
		}

		//return (index+1);
		return index;
	}
}

template <typename NT>
bool write_full_mesh_for_DFNWorks(const std::string &infilename_, const std::string &outfilename_,
									const std::list<int> &intList, const int &dimScaling, const std::string &lagrit_path_){

	//const std::string inp_file (infilename_ + ".inp");
	std::ifstream    infile(infilename_.c_str());

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file." << std::endl;
		return -1;
	}
	/*
	for (std::list<int>::const_iterator it = intList.begin(); it != intList.end(); it++){
		std::cout << "intList = " << *it << std::endl;
	}
	 */

	//const std::string out_file (outfilename_ + "_DFNWorks.inp");
	std::ofstream    outfile(outfilename_.c_str());

	// Vertices
	std::string line, linenospace;
	int precis(12), num_vertices, num_cells, num_data, dummy, vertices_per_cell;
	double matID;

	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.clear();

		ist >> num_vertices >> num_cells >> dummy;
	}

	outfile << line << std::endl;

	for (int i =0; i< num_vertices; i++){
		std::getline(infile, line);

		//outfile << line << std::endl;

		int cellNum;
		double vertex_0, vertex_1, vertex_2;
		std::istringstream ist(line);
		ist >> cellNum >> std::setprecision(precis) >> vertex_0 >> vertex_1 >> vertex_2;

		outfile << int2str_setw(cellNum, CGAL::count_digits(num_vertices)+1) << std::setprecision(precis) << std::scientific
				<< std::setw(precis+8) <<  vertex_0 * dimScaling
				<< std::setw(precis+8) <<  vertex_1 * dimScaling
				<< std::setw(precis+8) <<  vertex_2 * dimScaling << std::endl;
	}

	std::getline(infile, line);
	linenospace = line;
	linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
	linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
	unsigned first_ = line.find("1");

	// Cells
	int digits(first_+1), numfrac, vertice_, cell_;
	//std::cout << "digits = " << digits << std::endl;

	{
		std::istringstream ist(line);
		ist.clear();
		ist >> cell_ >> numfrac;
	}

	unsigned first = line.find("t");
	unsigned last = line.size();

	std::string strNew = line.substr (first+1,last-first);

	//std::cout << "CGAL::return_index_inf = " << CGAL::return_index_inf(intList, 1) << std::endl;

	outfile <<  CGAL::int2str_setw(cell_, digits) << std::setw(first-first_-4) << numfrac - CGAL::return_index_inf(intList,numfrac);
	outfile << "   t" << strNew << std::endl;

	for (int i = 0; i< num_cells -1; i++){
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		std::istringstream ist(line);
		ist.clear();
		ist >> cell_ >> numfrac;

		first = line.find("t");
		last = line.size();
		std::string strNew = line.substr (first+1,last-first);
		outfile <<  CGAL::int2str_setw(cell_, digits) << std::setw(first-first_-4) << numfrac - CGAL::return_index_inf(intList,numfrac);
		outfile << "   t" << strNew << std::endl;
	}

	std::getline(infile, line);
	{
		std::istringstream ist(line);
		ist.clear();
		ist >> num_data;
	}

	outfile << line << std::endl;
	for (int i =0; i< num_data; i++){
		std::getline(infile, line);
		outfile << line << std::endl;
	}
	std::getline(infile, line);

	linenospace = line;
	linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
	linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
	first_ = line.find("1");
	digits = first_ + 1;

	{
		std::istringstream ist(line);
		ist.clear();
		ist >> vertice_ >> matID;
	}

	first = line.find("E");
	last = line.size();

	strNew = line.substr (first+4,last-first-3);

	outfile <<  CGAL::int2str_setw(vertice_, digits) << "  "<<  std::setprecision(12)  << std::uppercase << std::scientific <<
			CGAL::to_double(round(matID) - CGAL::return_index_inf(intList,round(matID)));
	outfile << strNew << std::endl;

	for (int i = 0; i< num_vertices -1; i++){
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

		std::istringstream ist(line);
		ist.clear();
		ist >> vertice_ >> matID;

		first = line.find("E");
		last = line.size();
		strNew = line.substr (first+4,last-first-3);

		/*
		std::cout << "matID = " << matID << std::endl;
		std::cout << "CGAL::return_index_inf = " << CGAL::return_index_inf(intList, matID) << std::endl;
		*/

		outfile <<  CGAL::int2str_setw(vertice_, digits) << "  " <<  std::setprecision(12)  << std::uppercase << std::scientific <<
				CGAL::to_double(round(matID) - CGAL::return_index_inf(intList,round(matID)));
		outfile << strNew << std::endl;
	}

	for (int i = 0; i< num_cells +2; i++){
		std::getline(infile, line);
		outfile <<  line << std::endl;
	}

	infile.close();
	outfile.close();

	/*
	std::string outfile_1("write_new_tri_fracture.lgi");
	std::ofstream    outfile1;
	outfile1.open(outfile_1.c_str());


	outfile1 << "read / avs / full_mesh_DFNWorks.inp / mo\n";
	outfile1 << "dump / stor / tri_fracture / mo / ascii\n";
	outfile1 << "\n";
	outfile1 << "finish\n";
	outfile1.close();



	std::string cmd;
	std::string lagrit_path(lagrit_path_);
	cmd = lagrit_path + " < write_new_tri_fracture.lgi > write_new_tri_fracture.txt ";
	std::system(cmd.c_str());
	//std::system("rm -rf write_new_tri_fracture*");
	 */

	return true;
}
}

#endif /* INCLUDE_WRITE_FULL_MESH_FOR_DFNWORKS_HH_ */
