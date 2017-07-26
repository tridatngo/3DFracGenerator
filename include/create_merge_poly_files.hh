/*
 * create_merge_poly_files.hh
 *
 *  Created on: Sep 13, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CREATE_MERGE_POLY_FILES_HH_
#define INCLUDE_CREATE_MERGE_POLY_FILES_HH_

#include<list>
#include<string>
#include<iostream>
#include<sstream>
#include "myglobal_functions.hh"
#include "getNumLines.hh"

struct stat info_file;

namespace CGAL{

void create_merge_poly_files(const std::string &outfile_, const std::string &infile_, const int &numPoly, const int &digits){

	/*
	 * Create merge_poly file :
	 * 1. Creates a lagrit script that reads in each mesh, appends it to the main mesh, and then deletes that mesh object.
	 * 2. Then duplicate points are removed from the main mesh using EPS_FILTER.
	 * 3. The points are compressed, and then written in the files full_mesh.gmv, full_mesh.inp, and an FEHM dump is performed.
	 */

	std::ofstream    outfile;
	outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);

	outfile << "read / " << infile_ << " / mo_" << CGAL::int2str_setw(numPoly, digits) << " \n";
	outfile << "define / MO_NAME_FRAC / mo_" << CGAL::int2str_setw(numPoly, digits) << " \n";

	outfile << "cmo / addatt / MO_NAME_FRAC / volume / evol_onen \n";
	outfile << "math / sum / MO_NAME_FRAC / evol_sum / 1 0 0 / MO_NAME_FRAC / evol_one \n";

	outfile << "addmesh / merge / cmo_tmp / cmo_tmp / mo_" << CGAL::int2str_setw(numPoly, digits) << " \n";
	outfile << "cmo / delete / mo_" << CGAL::int2str_setw(numPoly, digits) << "\n";
	outfile << "\n";
	outfile.close();

	/*
	// This specific lines are used for parallel run

	outfile << " # Writing out merged fractures\n";

	outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
	outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
	outfile << "cmo select cmo_tmp\n";

	//outfile << "dump lagrit part%d.lg cmo_tmp\n";
	//outfile << "finish \n";

	//outfile << "read / lagrit / part%d.lg / junk / binary\n";

	outfile << "addmesh / merge / mo_all / mo_all / cmo_tmp\n";
	outfile << "cmo / delete / cmo_tmp\n";
	 */
}

template<typename Ell_>
void create_merge_poly_files_list(std::list<Ell_> ell_list){

	typedef	typename std::list<Ell_>::iterator Ell_Iterator;
	Ell_			cEll_;			// Current ellipse
	int 			digits( CGAL::count_digits(ell_list.size()) );

	std::cout << "Writing : merge_poly.lgi \n";
	std::string outfile_("merge_poly.lgi");

	for (int it = 0; it !=ell_list.size(); it++){
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;
		std::string tmp("mesh_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp");
		if (CGAL::getNumLines(std::string("ellipse_cutoff_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp")) > 3){
			CGAL::create_merge_poly_files(outfile_, tmp, cEll_.E_name + 1, digits);
		}
	}

	std::ofstream    outfile;
	outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);

	outfile << " # Writing out merged fractures\n";

	outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
	outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
	outfile << "cmo select cmo_tmp\n";

	outfile << "dump lagrit part_1.lg cmo_tmp\n";

	outfile << "finish \n";
	outfile.close();

	std::cout << "Writing : merge_poly.lgi : Complete. \n";
}

void create_merge_rmpts_file(double &EPS_INT, double &EPS_FILTER){

	std::string outfile_("merge_rmpts.lgi");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	/*
	// This specific lines are used for parallel run

	// Write lg file
	outfile << "read / " << outfile_ << " / mo_" << CGAL::int2str_setw(cEll_.E_name + 1,digits) << "\n";
	outfile << "define / MO_NAME_FRAC / mo_%d\n";
	outfile << "cmo / addatt / MO_NAME_FRAC / volume / evol_onen\n";
	outfile << "math / sum / MO_NAME_FRAC / evol_sum / 1 0 0 / MO_NAME_FRAC / evol_one \n";

	outfile << "addmesh / merge / cmo_tmp / cmo_tmp / mo_%d\n";
	outfile << "cmo / delete / mo_%d \n";

	// Read lg file
	for (int j = 1; j!=len(endis)+1; j++){
		outfile << "read / lagrit / part_" << CGAL::int2str(j)<< ".lg / junk / binary \n";
	}
	 */

	outfile << "read / lagrit / part_1.lg / junk / binary \n";
	outfile << "addmesh / merge / mo_all / mo_all / cmo_tmp\n";
	outfile << "cmo / delete / cmo_tmp\n";
	outfile << "\n";

	// Append meshes complete
	outfile << "# Appending the meshes complete\n";
	outfile << "# LaGriT Code to remove duplicates and output the mesh\n";
	outfile << "cmo / select / mo_all \n";
	outfile << "#recon 1\n";
	//outfile << "define / EPS / 1.e-5\n";
	outfile << "define / EPS / " << std::setprecision(3) <<std::scientific << EPS_INT << "\n";
	outfile << "define / EPS_FILTER / " << std::setprecision(3) <<std::scientific << EPS_FILTER << "\n";
	//outfile << "define / EPS_FILTER / 1.e-15\n";
	outfile << "pset / pinter / attribute / dfield / 1,0,0 / lt / EPS\n";
	outfile << "filter / pset get pinter / EPS_FILTER \n";
	outfile << "#filter / 1 0 0 / EPS_FILTER \n";
	outfile << "rmpoint / compress \n";
	outfile << "# SORT can affect a_b attribute\n";
	outfile << "sort / mo_all / index / ascending / ikey / imt xic yic zic\n";
	outfile << "reorder / mo_all / ikey \n";
	outfile << "cmo / DELATT / mo_all / ikey \n";
	outfile << "resetpts / itp \n";
	outfile << "boundary_components \n";

	// Export the final mesh
	outfile << "dump / full_mesh.gmv / mo_all \n";
	outfile << "dump / full_mesh.inp / mo_all \n";
	outfile << "dump / stor / tri_fracture / mo_all / ascii \n";

	outfile << "quality \n";
	outfile << "finish \n";

	//outfile << "dump / stor / tri_fracture / mo_all / ascii\n";
	outfile.close();
}

void run_merge_mesh_scripts(std::string &lagrit_path_){
	/*
		Triangulate each polygon. This python script runs lagrit to mesh each individual fractures
		Line Format std::system( "lagrit < merge_poly.lgi ")
	 */

	std::string cmd;
	std::string lagrit_path(lagrit_path_);

	cmd = lagrit_path + std::string(" < merge_poly.lgi > lagrit_logs/log_merge_poly");
	std::system(cmd.c_str());

	cmd = lagrit_path + std::string(" < merge_rmpts.lgi > lagrit_logs/log_merge_all");
	std::system(cmd.c_str());
}

//============================================================//
//============    Using multi-processing   ===================//
//============================================================//


template<typename Ell_>
void create_merge_poly_files_list_multi(std::list<Ell_> ell_list, const std::list<int> &endis){

	typedef	typename std::list<Ell_>::iterator Ell_Iterator;
	Ell_			cEll_;			// Current ellipse
	int 			digits( CGAL::count_digits(ell_list.size()) );

	for (int j = 1; j!=endis.size(); j++){
		std::cout << "Writing : merge_poly_part_"<< j <<".lgi \n";
		std::string outfile_("merge_poly_part_"+ int2str(j) +".lgi");

		typedef	typename std::list<int>::const_iterator Cst_Int_iterator;

		Cst_Int_iterator begin_frac = endis.begin();
		Cst_Int_iterator end_frac = endis.begin();

		std::advance(begin_frac,j-1);
		std::advance(end_frac,j);

		for (int it = *begin_frac; it != *end_frac; it++){
			Ell_Iterator	cEll_iter = ell_list.begin();
			std::advance(cEll_iter, it);
			cEll_ = *cEll_iter;
			std::string tmp("mesh_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp");
			if (CGAL::getNumLines(std::string("ellipse_cutoff_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp")) > 3){
				CGAL::create_merge_poly_files(outfile_, tmp, cEll_.E_name + 1, digits);
			}
		}

		std::ofstream    outfile;
		outfile.open(outfile_.c_str(), std::ofstream::out | std::ofstream::app);

		outfile << " # Writing out merged fractures\n";

		outfile << "mo / addatt/ cmo_tmp / volume / evol_all\n";
		outfile << "math / sum / cmo_tmp / evol_sum / 1 0 0 / cmo_tmp / evol_all\n";
		outfile << "cmo select cmo_tmp\n";

		outfile << "dump lagrit part_"<< j << ".lg cmo_tmp\n";

		outfile << "finish \n";
		outfile.close();

	}

	std::cout << "Writing : merge_poly.lgi : Complete. \n";
}

void create_merge_rmpts_file_multi(double &EPS_INT, double &EPS_FILTER, const std::list<int> &endis){

	std::string outfile_("merge_rmpts_multi.lgi");
	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	for (int j = 1; j!=endis.size(); j++){
		outfile << "read / lagrit / part_" << CGAL::int2str(j)<< ".lg / junk / binary \n";
		outfile << "addmesh / merge / mo_all / mo_all / cmo_tmp\n";
		outfile << "cmo / delete / cmo_tmp\n";
		outfile << "\n";
	}

	// Append meshes complete
	outfile << "# Appending the meshes complete\n";
	outfile << "# LaGriT Code to remove duplicates and output the mesh\n";
	outfile << "cmo / select / mo_all \n";
	outfile << "#recon 1\n";
	//outfile << "define / EPS / 1.e-5\n";
	outfile << "define / EPS / " << std::setprecision(3) <<std::scientific << EPS_INT << "\n";
	outfile << "define / EPS_FILTER / " << std::setprecision(3) <<std::scientific << EPS_FILTER << "\n";
	//outfile << "define / EPS_FILTER / 1.e-15\n";
	outfile << "pset / pinter / attribute / dfield / 1,0,0 / lt / EPS\n";
	outfile << "filter / pset get pinter / EPS_FILTER \n";
	outfile << "#filter / 1 0 0 / EPS_FILTER \n";
	outfile << "rmpoint / compress \n";
	outfile << "# SORT can affect a_b attribute\n";
	outfile << "sort / mo_all / index / ascending / ikey / imt xic yic zic\n";
	outfile << "reorder / mo_all / ikey \n";
	outfile << "cmo / DELATT / mo_all / ikey \n";
	outfile << "resetpts / itp \n";
	outfile << "boundary_components \n";

	// Export the final mesh
	outfile << "dump / full_mesh.gmv / mo_all \n";
	outfile << "dump / full_mesh.inp / mo_all \n";
	outfile << "dump / stor / tri_fracture / mo_all / ascii \n";

	outfile << "quality \n";
	outfile << "finish \n";

	//outfile << "dump / stor / tri_fracture / mo_all / ascii\n";
	outfile.close();
}


void run_merge_mesh_scripts_multi(std::string &lagrit_path_, const std::list<int> &endis){
	/*
     Merges all the meshes together, deletes duplicate points,
	 */

	std::string cmd;
	std::string lagrit_path(lagrit_path_);

	printf("\nMerging triangulated polygon meshes\n");

	int n_jobs = endis.size();
	pid_t pid;

	/*
	// should be converted to using multiprocessing
	for (int j =1; j != n_jobs + 1; j++){
		//pid = os.fork()
		if (pid == 0){ // clone a child job
			cmd = lagrit_path_ + " < merge_poly_part_" << CGAL::int2str(j) <<".lgi > log_merge_poly_part"<< CGAL::int2str(j);
			std::system(cmd.c_str());
			exit(0);
		}
		else{
			std::cout << "Merging part " << j << " of " << n_jobs << " job." << std::endl;
		}
	}

	// wait for all child processes to complete
	int j = 0;
	while (j < n_jobs){
		//(pid, status) = os.waitpid(0,os.WNOHANG);
		if (pid > 0){
			std::cout << "Process " + CGAL::int2str(j+1) + " finished." << std::endl;
			j += 1;
		}
	}
	*/
	for (int j = 1; j!=endis.size(); j++){
		cmd = lagrit_path + std::string(" < merge_poly_part_" + int2str(j) + ".lgi > lagrit_logs/log_merge_poly_part_" + int2str(j));
		std::system(cmd.c_str());
	}

	printf("Starting Final Merge\n");
	cmd = lagrit_path + std::string(" < merge_rmpts_multi.lgi > lagrit_logs/log_merge_all");
	std::system(cmd.c_str());

	// Check log_merge_all.txt for LaGriT complete successfully
	if (stat("full_mesh.inp", &info_file) == 0){
		printf("Final Merge Complete.\n");
		printf("Merging triangulated polygon meshes: Complete.\n");
	}
	else{
		printf("Final Merge Failed");
		exit(1);
	}
}

} // end of namespace CGAL

#endif /* INCLUDE_CREATE_MERGE_POLY_FILES_HH_ */
