/*
 * exe_read_vtu_write_txt_for_ADFNE.cpp
 *
 *  Created on: Jan 16, 2017
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

std::string int2str (const int &a){
	std::string out;
	std::ostringstream out_ss;
	out_ss << a;

	out = out_ss.str();

	return out;
}

std::string int2str_setw (const int &a,const int &w){
	std::string out;
	std::ostringstream out_ss;
	out_ss << std::setfill('0') << std::setw(w) << a;
	out = out_ss.str();

	return out;
}

int main(int argc, char** argv) {

	bool runtime_status(true);
	int fileNumber = 1E6;
	if( argc == 2 ) {
		fileNumber = atoi(argv[1]);
	}
	else {
		std::cout << "Usage: ./exefile fileNumber \n";
		std::cerr << "Introduce the number of *.vtk files." << std::endl;
		return -1;
	}


	// ===========================================================================
	// Read geometric data
	// ===========================================================================
	double  xmin(-1.5), xmax(1.5),
			ymin(-1.5), ymax(1.5),
			zmin(-1.5), zmax(1.5);
	double dx(xmax - xmin), dy(ymax - ymin), dz(zmax - zmin);

	system("rm -rf ADFNE");
    system("mkdir ADFNE");

    for (int it_ell=0; it_ell!=fileNumber; it_ell++){
    	std::cout << "Ellipse number " << it_ell << std::endl;
    	const std::string inp_file ("ellipse_within_dfnbb_" + int2str(it_ell)+ ".vtu");
    	std::cout << "Ellipse : " << inp_file << std::endl;
    	const char* filename = inp_file.c_str();
    	std::ifstream    infile(filename);

    	if (! infile.is_open()) {
    		std::cerr << "Error: Failed to open the input file." << std::endl;
    		return -1;
    	}

    	std::string line, linenospace;
    	int num_vertices;

    	for (int it=0; it!=3;it++){
    		std::getline(infile, line);
    	}
    	//<Piece NumberOfPoints="142" NumberOfCells="1">

    	{
    		unsigned first = line.find("Points=");
    		unsigned last = line.find("NumberOfCells");
    		//std::cout << "first = "<< first << std::endl;
    		//std::cout << "last = "<< last << std::endl;

    		std::string strNew = line.substr (first+8,last-first-3);
    		std::istringstream sstream(strNew);
    		sstream >> num_vertices;
    		//std::cout << "num_vertices = "<< num_vertices << std::endl;
    	}


    	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));

    	for (int it=0; it!=2;it++){
    		std::getline(infile, line);
    	}

    	for (int ivert=0; ivert!=num_vertices; ivert++)
    	{
    		std::getline(infile, line);
    		std::istringstream ist(line);
    		ist.str(line);
    		ist.clear();
    		ist >> vertices_[ivert][0] >> vertices_[ivert][1] >> vertices_[ivert][2];
    	}

    	infile.close();

    	std::string  outfile_("ADFNE/"+inp_file.substr(0,inp_file.length()-4)+".txt");
    	std::cout << "outfile_ = " << outfile_ << std::endl;
    	std::ofstream    outfile;
    	outfile.open(outfile_.c_str());
    	for (int it=0; it!=num_vertices; it++)
    	{
    		outfile << std::setprecision(15) << std::uppercase << std::scientific << std::setw(20)
    		<< (vertices_[it][0]-xmin)/dx << " " << (vertices_[it][1]-ymin)/dy  << " " << (vertices_[it][2]-zmin)/dz  << std::endl;
    	}
    	outfile.close();
    }
	return 0;
}


