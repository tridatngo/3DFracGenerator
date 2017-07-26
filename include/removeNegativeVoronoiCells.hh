/*
 * removeNegativeVoronoiCells.hh
 *
 *  Created on: Jan 5, 2017
 *      Author: ngotr
 */

#ifndef INCLUDE_REMOVENEGATIVEVORONOICELLS_HH_
#define INCLUDE_REMOVENEGATIVEVORONOICELLS_HH_

#include <list>
#include <stdio.h>      /* printf */
#include <math.h>       /* round, floor, ceil, trunc */
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>
#include "myglobal_functions.hh"
#include "getNumLines.hh"

namespace CGAL{

template<typename Ell_>
bool removeNegativeVoronoiCells(const std::list<Ell_> ell_list){

	int num_frac = ell_list.size();
	double MinVoronoiArea;

	std::string line, linenospace;
	std::list<int> fracNameList;
	int num_count_frac = 0;

	std::cout << "Writing : merge_poly_new.lgi \n";
	std::string outfile_("merge_poly_new.lgi");
	std::ofstream    outfile;
	//outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);
	outfile.open(outfile_.c_str());

	int digits(CGAL::count_digits(num_frac));

	for (int iter = 1; iter!=num_frac+1; iter++){
		std::string files("lagrit_logs/log_lagrit_" + CGAL::int2str_setw(iter,digits));
		//std::cout << "file : " <<  files << std::endl;

		std::ifstream    infile(files.c_str());
		if (! infile.is_open()) {
			std::cerr << "Failed to open the file "<< files << std::endl;
			exit(-1);
		}

		while (!infile.eof())
		{
			std::getline(infile, line);
			//std::cout << "line : " <<  line << std::endl;

			if (line.find("Minimum Voronoi area        =") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> MinVoronoiArea;
				//std::cout << "MinVoronoiArea = " << std::scientific << std::setprecision(6) << MinVoronoiArea << std::endl;
				if (MinVoronoiArea < 0) {
					fracNameList.push_back(MinVoronoiArea);
				}
				else{

					if (CGAL::getNumLines(std::string("ellipse_cutoff_" + CGAL::int2str_setw(iter,digits) + ".inp")) > 3){
						outfile << "read / mesh_" << CGAL::int2str_setw(iter,digits) << ".inp / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
						outfile << "define / MO_NAME_FRAC / mo_" << CGAL::int2str_setw(iter,digits) << " \n";

						outfile << "cmo / addatt / MO_NAME_FRAC / volume / evol_onen \n";
						outfile << "math / sum / MO_NAME_FRAC / evol_sum / 1 0 0 / MO_NAME_FRAC / evol_one \n";

						outfile << "addmesh / merge / cmo_tmp / cmo_tmp / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
						outfile << "cmo / delete / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
						outfile << "\n";
					}
				}
			}
		}

		infile.close();
	}

	outfile << " # Writing out merged fractures\n";

	outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
	outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
	outfile << "cmo select cmo_tmp\n";
	outfile << "dump lagrit part_1.lg cmo_tmp\n";
	outfile << "finish \n";

	//std::cout << "fracNameList.size() = " << fracNameList.size() << std::endl;

	/*
		fracNameList.sort();
		fracNameList.unique();

		std::ofstream    outfile(filename2.c_str());
		for (std::list<int>::iterator it = fracNameList.begin(); it != fracNameList.end(); it++){
			std::cout << "fracNameList = " << *it << std::endl;
			outfile << *it << std::endl;
		}
	 */

	/*
	// Path to lagrit executable
	//std::string lagrit_path( std::string("/work/irlin104_1/ngotr/Outils/lagrit/lagrit_ulin3.2") );
	std::string lagrit_path(lagrit_path_);

	std::string cmd;
	cmd = lagrit_path + std::string(" < merge_poly_new.lgi > lagrit_logs/log_merge_poly_new");
	std::system(cmd.c_str());

	//std::system("cp merge_rmpts.lgi merge_rmpts_new.lgi");
	//std::system("sed -i 's/full_mesh/full_mesh_new/' merge_rmpts_new.lgi");
	//cmd = lagrit_path + std::string(" < merge_rmpts_new.lgi > lagrit_logs/log_merge_all_new");

	cmd = lagrit_path + std::string(" < merge_rmpts.lgi > lagrit_logs/log_merge_all_new");
	std::system(cmd.c_str());
	*/

	return 1;
}

template<typename Ell_>
bool removeNegativeVoronoiCells_multi(const std::list<Ell_> ell_list, const std::list<int> endis){

	int num_frac = ell_list.size();
	double MinVoronoiArea;

	std::string line, linenospace;
	std::list<int> fracNameList;
	int num_count_frac = 0;

	for (int j = 1; j!=endis.size(); j++){
		typedef	typename std::list<int>::const_iterator Cst_Int_iterator;

		Cst_Int_iterator begin_frac = endis.begin();
		Cst_Int_iterator end_frac = endis.begin();

		std::advance(begin_frac,j-1);
		std::advance(end_frac,j);

		std::cout << "Writing : merge_poly_part_"<< j <<"_new.lgi \n";
		std::string outfile_("merge_poly_part_"+ int2str(j) +"_new.lgi");

		std::ofstream    outfile;
		//outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);
		outfile.open(outfile_.c_str());

		int digits(CGAL::count_digits(num_frac));

		for (int iter = 1; iter!=num_frac+1; iter++){
			std::string files("lagrit_logs/log_lagrit_" + CGAL::int2str_setw(iter,digits));
			//std::cout << "file : " <<  files << std::endl;

			std::ifstream    infile(files.c_str());
			if (! infile.is_open()) {
				std::cerr << "Failed to open the file "<< files << std::endl;
				exit(-1);
			}

			while (!infile.eof())
			{
				std::getline(infile, line);
				//std::cout << "line : " <<  line << std::endl;

				if (line.find("Minimum Voronoi area        =") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> MinVoronoiArea;
					//std::cout << "MinVoronoiArea = " << std::scientific << std::setprecision(6) << MinVoronoiArea << std::endl;
					if (MinVoronoiArea < 0) {
						fracNameList.push_back(MinVoronoiArea);
					}
					else{

						if (CGAL::getNumLines(std::string("ellipse_cutoff_" + CGAL::int2str_setw(iter,digits) + ".inp")) > 3){
							outfile << "read / mesh_" << CGAL::int2str_setw(iter,digits) << ".inp / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
							outfile << "define / MO_NAME_FRAC / mo_" << CGAL::int2str_setw(iter,digits) << " \n";

							outfile << "cmo / addatt / MO_NAME_FRAC / volume / evol_onen \n";
							outfile << "math / sum / MO_NAME_FRAC / evol_sum / 1 0 0 / MO_NAME_FRAC / evol_one \n";

							outfile << "addmesh / merge / cmo_tmp / cmo_tmp / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
							outfile << "cmo / delete / mo_" << CGAL::int2str_setw(iter,digits) << " \n";
							outfile << "\n";
						}
					}
				}
			}

			infile.close();
		}

		outfile << " # Writing out merged fractures\n";

		outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
		outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
		outfile << "cmo select cmo_tmp\n";
		outfile << "dump lagrit part_1.lg cmo_tmp\n";
		outfile << "finish \n";
	}

	return 1;
}
}

#endif /* INCLUDE_REMOVENEGATIVEVORONOICELLS_HH_ */
