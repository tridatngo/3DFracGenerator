/*
 * create_parameter_mlgi_files.hh
 *
 *  Created on: Sep 9, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CREATE_PARAMETER_MLGI_FILES_HH_
#define INCLUDE_CREATE_PARAMETER_MLGI_FILES_HH_

#include<list>
#include <sys/types.h>
#include <sys/stat.h>
#include "myglobal_functions.hh"

//using namespace std;

struct stat info_1;

namespace CGAL{

template <typename Ell_, typename NT>
bool create_parameter_mlgi_files(std::list<Ell_> & ell_list, const std::string &filename_) {

	typedef typename std::list<Ell_>::iterator Ell_Iterator;

	/*  #Section 2 : Outputs parameter_i.mlgi files used in running LaGriT Script */
	std::cout<< "\nCreating parameter_i.mlgi files" << std::endl;

	int 	E_name, num_ell, digits(CGAL::count_digits(ell_list.size()) );
	NT theta, dL0, x1, y1, z1, x2, y2, z2;

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	std::istringstream ist(line);
	//std::getline(infile, line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	for (int it =0; it!=num_ell; it++){
		E_name = it;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> E_name >> theta >> dL0 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;

		//std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::ofstream   outfile;
		outfile.open(outfile_.c_str());

		if (! outfile.is_open()) {
			std::cerr << "Failed to open the input file " << outfile_ << "." << std::endl;
			return -1;
		}

		// header
		outfile << "define / ID / " << E_name << "\n";
		outfile << "define / OUTFILE_GMV / mesh_" << CGAL::int2str_setw(E_name,digits) <<  ".gmv\n";
		outfile << "define / OUTFILE_AVS / mesh_" <<  CGAL::int2str_setw(E_name,digits) << ".inp\n";
		//outfile << "define / POLY_FILE / poly_" <<  CGAL::int2str_setw(E_name,digits) <<  ".inp\n";
		//outfile << "define / QUAD_FILE / tmp_quad_" <<  E_name <<  ".inp\n";
		//outfile << "define / EXCAVATE_FILE / tmp_excavate_" <<  E_name <<  ".inp\n";
		//outfile << "define / PRE_FINAL_FILE / tmp_pre_final_"<< E_name <<  ".inp\n";
		//outfile << "define / PRE_FINAL_MASSAGE / tmp_pre_final_massage_" <<  E_name << ".gmv\n";
		
		// NTD 03/11/2016 [IMP --> change fracture connectivity when using std::setprecision(20), to be improve]
		outfile << "define / THETA  / " << -theta << "\n";
		//outfile << "define / THETA  / " << std::setprecision(20) << -theta << "\n";
		
		/*
		std::cout << "theta = " << theta << std::endl;
		if (theta < -90.1 && theta>-100){
			outfile << "define / THETA  / " << std::setprecision(20) << theta << "\n";
		}
		else{
			outfile << "define / THETA  / " << std::setprecision(20) << -theta << "\n";
		}
		*/
		outfile << "define / dL0  / " << std::setprecision(20) << dL0 << "\n";
		outfile << "define / X1 / " << std::setprecision(20) << x1 << "\n";
		outfile << "define / Y1 / " << std::setprecision(20) << y1 << "\n";
		outfile << "define / Z1 / " << std::setprecision(20) << z1 << "\n";
		outfile << "define / X2 / " << std::setprecision(20) << x2 << "\n";
		outfile << "define / Y2 / " << std::setprecision(20) << y2 << "\n";
		outfile << "define / Z2 / " << std::setprecision(20) << z2 << "\n";
		outfile << "finish \n";

		outfile.close();

		std::getline(infile, line);
	}

	infile.close();

	std::cout << "Creating parameter_i.mlgi files: Complete." << std::endl;

	return true;
}


template <typename Ell_, typename NT>
bool create_parameter_mlgi_files_vercors(std::list<Ell_> & ell_list, const std::string &filename_) {
	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	/*  #Section 2 : Outputs parameter_i.mlgi files used in running LaGriT Script */
	std::cout<< "\nCreating parameter_i.mlgi files" << std::endl;

	int 	E_name, num_ell, digits(CGAL::count_digits(ell_list.size()) );
	NT theta, dL0, x1, y1, z1, x2, y2, z2;

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	std::istringstream ist(line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	for (int it =0; it!=num_ell; it++){
		E_name = it;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> E_name >> theta >> dL0 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;

		//std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::ofstream   outfile;
		outfile.open(outfile_.c_str());

		if (! outfile.is_open()) {
			std::cerr << "Failed to open the input file " << outfile_ << "." << std::endl;
			return -1;
		}

		// header
		outfile << "define / ID / " << E_name << "\n";
		outfile << "define / OUTFILE_GMV / mesh_" << CGAL::int2str_setw(E_name,digits) <<  ".gmv\n";
		outfile << "define / OUTFILE_AVS / mesh_" <<  CGAL::int2str_setw(E_name,digits) << ".inp\n";
		// NTD 03/11/2016 [IMP --> change fracture connectivity when using std::setprecision(20), to be improve]
		//outfile << "define / THETA  / " << std::setprecision(20) << theta << "\n";
		outfile << "define / THETA  / " << theta << "\n";
		outfile << "define / dL0  / " << std::setprecision(20) << dL0 << "\n";
		outfile << "define / X1 / " << std::setprecision(20) << x1 << "\n";
		outfile << "define / Y1 / " << std::setprecision(20) << y1 << "\n";
		outfile << "define / Z1 / " << std::setprecision(20) << z1 << "\n";
		outfile << "define / X2 / " << std::setprecision(20) << x2 << "\n";
		outfile << "define / Y2 / " << std::setprecision(20) << y2 << "\n";
		outfile << "define / Z2 / " << std::setprecision(20) << z2 << "\n";
		outfile << "finish \n";

		outfile.close();

		std::getline(infile, line);
	}

	infile.close();

	std::cout << "Creating parameter_i.mlgi files: Complete." << std::endl;

	return true;
}


template <typename Ell_, typename NT>
bool create_parameter_mlgi_files_setprecision(std::list<Ell_> & ell_list, const std::string &filename_) {

	typedef typename std::list<Ell_>::iterator Ell_Iterator;

	/*  #Section 2 : Outputs parameter_i.mlgi files used in running LaGriT Script */
	std::cout<< "\nCreating parameter_i.mlgi files" << std::endl;

	int 	E_name, num_ell, digits(CGAL::count_digits(ell_list.size()) );
	NT theta, dL0, x1, y1, z1, x2, y2, z2;

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	std::istringstream ist(line);
	//std::getline(infile, line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	for (int it =0; it!=num_ell; it++){
		E_name = it;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> E_name >> theta >> dL0 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;

		//std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::ofstream   outfile;
		outfile.open(outfile_.c_str());

		if (! outfile.is_open()) {
			std::cerr << "Failed to open the input file " << outfile_ << "." << std::endl;
			return -1;
		}

		// header
		outfile << "define / ID / " << E_name << "\n";
		outfile << "define / OUTFILE_GMV / mesh_" << CGAL::int2str_setw(E_name,digits) <<  ".gmv\n";
		outfile << "define / OUTFILE_AVS / mesh_" <<  CGAL::int2str_setw(E_name,digits) << ".inp\n";
		//outfile << "define / POLY_FILE / poly_" <<  CGAL::int2str_setw(E_name,digits) <<  ".inp\n";
		//outfile << "define / QUAD_FILE / tmp_quad_" <<  E_name <<  ".inp\n";
		//outfile << "define / EXCAVATE_FILE / tmp_excavate_" <<  E_name <<  ".inp\n";
		//outfile << "define / PRE_FINAL_FILE / tmp_pre_final_"<< E_name <<  ".inp\n";
		//outfile << "define / PRE_FINAL_MASSAGE / tmp_pre_final_massage_" <<  E_name << ".gmv\n";

		// NTD 03/11/2016 [IMP --> change fracture connectivity when using std::setprecision(20), to be improve]
		//outfile << "define / THETA  / " << -theta << "\n";
		outfile << "define / THETA  / " << std::setprecision(20) << -theta << "\n";
		outfile << "define / dL0  / " << std::setprecision(20) << dL0 << "\n";
		outfile << "define / X1 / " << std::setprecision(20) << x1 << "\n";
		outfile << "define / Y1 / " << std::setprecision(20) << y1 << "\n";
		outfile << "define / Z1 / " << std::setprecision(20) << z1 << "\n";
		outfile << "define / X2 / " << std::setprecision(20) << x2 << "\n";
		outfile << "define / Y2 / " << std::setprecision(20) << y2 << "\n";
		outfile << "define / Z2 / " << std::setprecision(20) << z2 << "\n";
		outfile << "finish \n";

		outfile.close();

		std::getline(infile, line);
	}

	infile.close();

	std::cout << "Creating parameter_i.mlgi files: Complete." << std::endl;

	return true;
}


template <typename Ell_, typename NT>
bool create_parameter_mlgi_files_vercors_setprecision(std::list<Ell_> & ell_list, const std::string &filename_) {
	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	/*  #Section 2 : Outputs parameter_i.mlgi files used in running LaGriT Script */
	std::cout<< "\nCreating parameter_i.mlgi files" << std::endl;
	//std::system("rm -rf output/parameters");
	//std::system("mkdir output/parameters");

	int 	E_name, num_ell, digits(CGAL::count_digits(ell_list.size()) );
	NT theta, dL0, x1, y1, z1, x2, y2, z2;

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		return -1;
	}

	std::string line;
	std::istringstream ist(line);

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	for (int it =0; it!=num_ell; it++){

		E_name = it;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> E_name >> theta >> dL0 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;

		if (VERY_VERBOSE == 1){
			std::cout << E_name << " " << theta << " " << dL0 << " " << x1 << " " << y1 << " " << z1
					<< " " << x2 << " " << y2 << " " << z2 << std::endl;
		}

		if (stat("output/parameters", &info_1) != 0){
			std::system("mkdir output/parameters");
		}

		//std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::string  	outfile_("output/parameters/parameters_" + CGAL::int2str_setw(E_name,digits) + ".mlgi");
		std::ofstream   outfile;
		outfile.open(outfile_.c_str());

		if (! outfile.is_open()) {
			std::cerr << "Failed to open the input file " << outfile_ << "." << std::endl;
			return -1;
		}

		// header
		outfile << "define / ID / " << E_name << "\n";
		outfile << "define / OUTFILE_GMV / mesh_" << CGAL::int2str_setw(E_name,digits) <<  ".gmv\n";
		outfile << "define / OUTFILE_AVS / mesh_" <<  CGAL::int2str_setw(E_name,digits) << ".inp\n";
		// NTD 03/11/2016 [IMP --> change fracture connectivity when using std::setprecision(20), to be improve]
		outfile << "define / THETA  / " << std::setprecision(20) << theta << "\n";
		//outfile << "define / THETA  / " << theta << "\n";
		outfile << "define / dL0  / " << std::setprecision(20) << dL0 << "\n";
		outfile << "define / X1 / " << std::setprecision(20) << x1 << "\n";
		outfile << "define / Y1 / " << std::setprecision(20) << y1 << "\n";
		outfile << "define / Z1 / " << std::setprecision(20) << z1 << "\n";
		outfile << "define / X2 / " << std::setprecision(20) << x2 << "\n";
		outfile << "define / Y2 / " << std::setprecision(20) << y2 << "\n";
		outfile << "define / Z2 / " << std::setprecision(20) << z2 << "\n";
		outfile << "finish \n";

		outfile.close();

		std::getline(infile, line);
	}

	infile.close();

	std::cout << "Creating parameter_i.mlgi files: Complete." << std::endl;

	return true;
}
} // end of namespace

#endif /* INCLUDE_CREATE_PARAMETER_MLGI_FILES_HH_ */
