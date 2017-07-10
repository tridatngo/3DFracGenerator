/*
 * exe_read_inp_write_dgf_file_hete.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: ngotr
 */

//#include <boost/fusion/iterator/next.hpp>
//#include <boost/fusion/include/next.hpp>
//#include <boost/fusion/iterator/prior.hpp>
//#include <boost/fusion/include/prior.hpp>

#include <string>     // std::string, std::to_string
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iomanip> // setprecision
#include <unistd.h>
#include <math.h>       /* fabs */

#include <list>
#include <algorithm>
#include <functional>
#include <numeric>

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

struct Apert_ID
{
	double Apert_;
	int ID_;
};

bool compareID(const Apert_ID &lhs, const Apert_ID &rhs)
{
	return lhs.ID_ < rhs.ID_;
}

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


template < typename T>
T MyNorm (const std::vector<T> & vect)
{
	return sqrt( vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
}

template <typename T>
std::vector<T> MyNormalize (std::vector<T> const &vect)
{

	/*
	 if (MyNorm<T>(vect) == 0) {
		std::cout << "Vector " << vect << " is a zero vector. Return the unit vector (0,0,1)." <<std::endl;
	}
	else{
	 */
	std::vector<T> vect_norm(3);
	if (MyNorm<T>(vect) != 0) {
		vect_norm[0] = vect[0] / MyNorm<T>(vect);
		vect_norm[1] = vect[1] / MyNorm<T>(vect);
		vect_norm[2] = vect[2] / MyNorm<T>(vect);
	}
	return vect_norm;
}


template <typename T>
std::vector<T> MyCrossProduct(std::vector<T> const &a, std::vector<T> const &b)
{
	std::vector<T> r(a.size());
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
	T len;
	len = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

	r[0] = r[0]/len;
	r[1] = r[1]/len;
	r[2] = r[2]/len;

	return r;
}

template < typename T>
T MyAngle (std::vector<T> const &v1_, std::vector<T> const &v2_)
{
	/*
	NT angle_ = acos (CGAL::to_double(CGAL::MyNormalize<K, NT>(v1_) *  CGAL::MyNormalize<K, NT>(v2_)));
	return angle_;
	 */

	// NTD 10/11/2016
	if ( MyNorm<T>( MyCrossProduct<T>(v1_, v2_) ) == 0 ){
		return 0;
	}
	else{
		//return acos (MyNormalize<T>(v1_) * MyNormalize<T>(v2_));

		std::vector<T> a(v1_.size()), b(v2_.size());

		a = MyNormalize<T>(v1_);
		b = MyNormalize<T>(v2_);
		
		//std::cout << "======"<< std::endl;
		//std::cout << "MyNormalize<T>(v1_) = " << MyNormalize<T>(v1_)[0] << " " << MyNormalize<T>(v1_)[1] << " " << MyNormalize<T>(v1_)[2] << std::endl;
		//std::cout << "MyNormalize<T>(v2_) = " << MyNormalize<T>(v2_)[0] << " " << MyNormalize<T>(v2_)[1] << " " << MyNormalize<T>(v2_)[2] << std::endl;
		//std::cout << "std::inner_product(a.begin(), a.end(), b.begin(), 0) = " << std::inner_product(a.begin(), a.end(), b.begin(), 0.0) << std::endl;
		//std::cout << "MyNormalize<T>(v2_) = " << MyNormalize<T>(v2_) << std::endl;

		return acos(std::inner_product(a.begin(), a.end(), b.begin(),0.0));
	}
}

// 3D rotation

typedef struct {
	double x;
	double y;
	double z;
}Point;
Point points;

float rotationMatrix[4][4];
float inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
float outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};

void multiplyMatrix()
{
    for(int i = 0; i < 4; i++ ){
        for(int j = 0; j < 1; j++){
            outputMatrix[i][j] = 0;
            for(int k = 0; k < 4; k++){
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
}

void setUpRotationMatrix(const double &theta, const std::vector<double> &rotvect)
{
	double u(rotvect[0]);
	double v(rotvect[1]);
	double w(rotvect[2]);

    double L = (u*u + v * v + w * w);
    //theta = theta * M_PI / 180.0; //converting to radian value
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;
	/*
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(theta)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(theta)) - w * sqrt(L) * sin(theta)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(theta)) + v * sqrt(L) * sin(theta)) / L;
    rotationMatrix[0][3] = 0.0;

    rotationMatrix[1][0] = (u * v * (1 - cos(theta)) + w * sqrt(L) * sin(theta)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(theta)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(theta)) - u * sqrt(L) * sin(theta)) / L;
    rotationMatrix[1][3] = 0.0;

    rotationMatrix[2][0] = (u * w * (1 - cos(theta)) - v * sqrt(L) * sin(theta)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(theta)) + u * sqrt(L) * sin(theta)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(theta)) / L;
    rotationMatrix[2][3] = 0.0;
	*/
	double c = cos(theta);
	double s = sin(theta);
	
	rotationMatrix[0][0] = u * u *(1-c) + c;
	rotationMatrix[0][1] = v * u *(1-c) - w*s;
	rotationMatrix[0][2] = w * u *(1-c) + v*s;
	rotationMatrix[0][3] = 0.0;
	
	rotationMatrix[1][0] = u * v *(1-c) + w*s;
	rotationMatrix[1][1] = v * v *(1-c) + c;
	rotationMatrix[1][2] = w * v *(1-c) - u*s;
	rotationMatrix[1][3] = 0.0;

	rotationMatrix[2][0] = u * w *(1-c) - v*s;
	rotationMatrix[2][1] = v * w *(1-c) + u*s;
	rotationMatrix[2][2] = w * w *(1-c) + c;
	rotationMatrix[2][3] = 0.0;
	
    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;
}

void showPoint(){
    std::cout<<"("<<outputMatrix[0][0]<<","<<outputMatrix[1][0]<<","<<outputMatrix[2][0]<<")"<<std::endl;
}

/*
template <typename K, typename NT>
CGAL::Aff_transformation_3<K> MyRotationMatrix (CGAL::Vector_3<K> &normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis
	Vector_3 rotvect = CGAL::MyCross<K>(normal_vect, e3);
	Vector_3 u = CGAL::Normalize<K,NT>(rotvect);

	// NTD - 25/10/2016
	//NT theta = CGAL::MyAngle<Kernel, NT_MP>(normal_vect,e3);
	NT theta = CGAL::MyAngle<K, NT>(normal_vect,e3);
	NT c = cos(CGAL::to_double(theta));
	NT s = sin(CGAL::to_double(theta));

	Aff_transformation_3 rotMat(
			u[0]*u[0]*(1-c) +      c,  u[1]*u[0]*(1-c) - u[2]*s,  u[2]*u[0]*(1-c) + u[1]*s,  transvect[0],
			u[0]*u[1]*(1-c) + u[2]*s,  u[1]*u[1]*(1-c) +      c,  u[2]*u[1]*(1-c) - u[0]*s,  transvect[1],
			u[0]*u[2]*(1-c) - u[1]*s,  u[1]*u[2]*(1-c) + u[0]*s,  u[2]*u[2]*(1-c) +      c,  transvect[2],
			1);
	return rotMat;
}
*/

// Compute permeability from the fracture aperture
double compute_K_from_apert(const double &b){
	/*
	 * K = b^2/12
	 */

	return pow(b,2)/12;
}

double compute_K_from_apert_mul(const double &b){
	/*
	 * K = alpha * b^2/12
	 */
	double alpha =1.2E-4; // alpha : multiplicative factor

	return alpha * pow(b,2)/12;
}

// Compute transmissivity from the fracture aperture
double compute_sigma_from_apert(const double &b){
	double rho = 1E3; // [kg/m^3] fluid density
	double mu = 1E-3; // [Pa.s] fluid viscosity
	double g = 9.81; // [m/s^2] gravity acceleration

	return pow(b,3)*rho*g/12/mu;
}

double compute_sigma_from_apert_mul(const double &b){
	double rho = 1E3; // [kg/m^3] fluid density
	double mu = 1E-3; // [Pa.s] fluid viscosity
	double g = 9.81; // [m/s^2] gravity acceleration
	double alpha = 1.2E-4; // alpha : multiplicative factor

	return alpha * pow(b,3)*rho*g/12/mu;
}


void find_and_replace(std::string& source, std::string const& find, std::string const& replace)
{
	for(std::string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
	{
		source.replace(i, find.length(), replace);
		i += replace.length();
	}
}
//=============================================================================================
// Main program
//=============================================================================================

int main(int argc, char** argv) {

	//char* a_cwd = _getcwd(NULL, 0);
	bool bool_linux(true);
	
	bool runtime_status(true);
	std::string filename_ = "";
	std::string inputfile_ = "";
	std::string cwd = "";
	std::string s_cwd= "";

	if( argc == 3 ) {
		filename_ = argv[1];
		inputfile_ = argv[2];
	}
	else if ( argc == 4 ){
		filename_ = argv[1];
		inputfile_ = argv[2];
		cwd = argv[3];
	}
	else {
		std::cout << "Usage: \"./exefile InpInputFile DataInputFile\" or \"./exefile InpInputFile DataInputFile CurrentWorkingDirectory\"\n";
		std::cerr << "Failed to open the input files." << std::endl;
		return -1;
	}

	if( argc == 3 ) {
		char* a_cwd = getcwd(NULL, 0);
		std::string str_cwd(a_cwd);
		s_cwd = str_cwd;
		std::replace( s_cwd.begin(), s_cwd.end(), '\\', '\/');
		std::cout << "Current working direction:" << s_cwd << std::endl;
	}
	else{
		s_cwd = cwd;
		std::cout << "Current working direction:" << s_cwd << std::endl;
	}

	std::string apertureFile(std::string(s_cwd + "/aperture.dat")), permFile(std::string(s_cwd + "/permeability.dat"));
	find_and_replace( apertureFile, "//", "/");
	find_and_replace( permFile, "//", "/");

	/*
	 * hete = 0: homogeneous
	 * hete = 1: fracture-to-fracture heterogeneity
	 * hete = 2: cell-to-cell heterogeneity
	 */
	int hete(2);
	std::string input = "";

	while (true) {
		std::cout << "Please enter a valid number of hete: "<< std::endl;
		std::cout << "* hete = 0: homogeneous medium." << std::endl;
		std::cout << "* hete = 1: fracture-to-fracture heterogeneity." << std::endl;
		std::cout << "* hete = 2: cell-to-cell heterogeneity." << std::endl;
		std::cout << "*** hete = ";
		getline(std::cin, input);

		// This code converts from string to number safely.
		std::stringstream myStream(input);
		if (myStream >> hete){
			if (hete > 2){
				std::cout << "Invalid option: hete > 2. ";
			}
			else if(hete <0){
				std::cout << "Invalid option: hete < 0. ";
			}
			else
				break;
		}
		std::cout << "Invalid number, please try again" << std::endl;
	}

	// ===========================================================================
	// Read input data about fracture characteristics
	// ===========================================================================
	int 	E_name(0), num_fracs(0), num_ell(0);								// Name of the ellipse
	double 	E_center_0, E_center_1, E_center_2 /*Center of the ellipse*/,
	E_xradius, E_yradius, /* x_radius, y_radius */
	E_normalvect_0, E_normalvect_1, E_normalvect_2 /* normal vector */;

	std::list<int>						E_name_list;		// Name of the ellipse
	std::list< std::vector<double> >	E_center_list;		// Center of the ellipse
	std::list<double>					E_radius_list;		// x_radius, y_radius
	std::list< std::vector<double> >	E_normalvect_list;	// normal vector
	
	std::ifstream    data_infile(inputfile_.c_str());
	if (! data_infile.is_open()) {
		std::cerr << "Error: Failed to open the data input file." << std::endl;
		return -1;
	}

	std::string line, linenospace;
	std::cout << "Reading data input file \"" << inputfile_ << "\" ..."<< std::endl;
	//std::cout << "Line : " << line << std::endl;
	std::getline(data_infile, line);
	
	while(!data_infile.eof()){

		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> 	E_center_0 >> E_center_1 >> E_center_2 >>
		E_xradius >> E_yradius >>
		E_normalvect_0 >> E_normalvect_1 >> E_normalvect_2;

		std::vector<double> eCenter(3), enormVec(3);
		eCenter[0] = E_center_0; eCenter[1] = E_center_1; eCenter[2] = E_center_2;
		enormVec[0] = E_normalvect_0; enormVec[1] = E_normalvect_1; enormVec[2] = E_normalvect_2;

		//std::cout << "Line : " << line << std::endl;
		//std::cout << "E_center_0 : " << E_center_0 << std::endl;
		//std::cout << "E_normalvect_0 : " << E_normalvect_0 << std::endl;

		E_name_list.push_back(num_fracs);
		E_center_list.push_back(eCenter);
		E_radius_list.push_back( std::min(E_xradius, E_yradius) );
		E_normalvect_list.push_back(enormVec);
		
		std::getline(data_infile, line);
		num_fracs++;
	}
	std::cout << "num_fracs = " << num_fracs << std::endl;

	/*
	for (int ifrac=0; ifrac < num_fracs; ifrac++){
		std::list<double>::iterator E_radius = E_radius_list.begin();
		std::advance(E_radius, ifrac);
		std::cout << "E_radius["<< ifrac <<"].size() = " << *E_radius << std::endl;
	}
	*/

	std::vector< std::list<std::vector<double> > > 	Points_list_withID(num_fracs);
	std::vector< std::list<double> >				Aper_list_withID(num_fracs); // list of aperture for each fracture
	std::vector< std::list<int> > 					Global_cellID(num_fracs); // list of global ID for each fracture
	
	std::cout << "Points_list_withID.size() = " << Points_list_withID.size() << std::endl;

	// ===========================================================================
	// Read parameter data (aperture and permeability)
	// ===========================================================================

	std::list<double> permeability_X, permeability_Y, permeability_Z;
	std::ifstream    paramfile1(permFile.c_str());

	if (! paramfile1.is_open()) {
		std::cerr << "Error: Failed to open the permeability input file." << std::endl;
		return -1;
	}

	//int num_fracs;
	//num_fracs = getnumline(permFile);
	for (int i=0; i!= num_fracs; i++){
		std::getline(paramfile1,line);
		std::istringstream ist(line);
		double dummy_X, dummy_Y, dummy_Z;
		ist >> dummy_X >> dummy_Y >> dummy_Z;
		permeability_X.push_back(dummy_X);
		permeability_Y.push_back(dummy_Y);
		permeability_Z.push_back(dummy_Z);
	}

	std::vector<double> aperture(num_fracs);
	std::ifstream    paramfile2(apertureFile.c_str());

	if (! paramfile2.is_open()) {
		std::cerr << "Error: Failed to open the aperture input file." << std::endl;
		return -1;
	}

	for (int i=0; i!= num_fracs; i++){
		std::getline(paramfile2,line);
		std::istringstream ist(line);
		double dummy;
		ist >> dummy;
		aperture[i] = dummy;
	}
	// ===========================================================================
	// Read geometric data
	// ===========================================================================

	std::string inp_file (s_cwd + "/" + filename_.substr(0,filename_.length()-4)+ ".inp");
	find_and_replace( inp_file, "//", "/");
	const char* filename = inp_file.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file." << std::endl;
		return -1;
	}

	int num_vertices, num_cells, dummy, vertices_per_cell;

	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist.clear();

		ist >> num_vertices >> num_cells >> dummy;
	}

	std::vector< std::vector<double> > vertices_(num_vertices, std::vector<double>(3));
	std::vector<int> material_ID(num_cells), check_material_ID(num_cells);

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

	int int_dummy_1, int_dummy_2;
	std::string elem_type;

	ist >> int_dummy_1 >> int_dummy_2 >> elem_type;

	material_ID[0] = int_dummy_2;

	if ( elem_type == std::string("tri")){
		vertices_per_cell = 3;
	}
	else if (elem_type == std::string("quad")) {
		vertices_per_cell = 4;
	}
	//==================================================================================================
	//  Define new variables
	//==================================================================================================

	//std::cout << "Number of vertices per cell, vertices_per_cell = " << vertices_per_cell << std::endl;

	std::vector< std::vector<int> > 		cells_(num_cells, std::vector<int>(vertices_per_cell));
	std::vector< std::vector<double> > 		cell_centers(num_cells, std::vector<double>(3));
	std::vector< int > 						Arr_Global_cellID(num_cells); // Array of Global ID for each cell
	std::vector< double >					Arr_Aper_Hete(num_cells); // Array of aperture for each cell
	std::vector< double >					Arr_Perm_Hete(num_cells); // Array of aperture for each cell
	
	//std::vector< std::vector<double> > cell_centers(num_cells, std::vector<double>(3));

	double cc_dummy_x(0), cc_dummy_y(0), cc_dummy_z(0); // X, Y, Z coordianate of the current cell center

	for (int i =0; i < vertices_per_cell; i++){

		ist >> cells_[0][i];

		int icell = cells_[0][i]-1;
		cc_dummy_x = cc_dummy_x + vertices_[icell][0];
		cc_dummy_y = cc_dummy_y + vertices_[icell][1];
		cc_dummy_z = cc_dummy_z + vertices_[icell][2];

		//std::cout << "vertices_[icell][0] = " << vertices_[icell][0] << "\n";
		//std::cout << "vertices_[icell][1] = " << vertices_[icell][1] << "\n";
		//std::cout << "vertices_[icell][2] = " << vertices_[icell][2] << "\n";
	}

	cell_centers[0][0] = cc_dummy_x/vertices_per_cell;
	cell_centers[0][1] = cc_dummy_y/vertices_per_cell;
	cell_centers[0][2] = cc_dummy_z/vertices_per_cell;
	cc_dummy_x = 0; cc_dummy_y = 0; cc_dummy_y =0;

	//std::cout << "** cell_centers(0):" << cell_centers[0][0] << " " << cell_centers[0][1] << " " << cell_centers[0][2] << "\n";

	for (int icells=1; icells < num_cells; icells++)
	{
		std::getline(infile, line);
		std::istringstream ist(line);
		ist >> int_dummy_1 >> int_dummy_2 >> elem_type;

		material_ID[icells] = int_dummy_2;

		for (int i =0; i < vertices_per_cell; i++)
		{
			ist >> cells_[icells][i];

			int icell = cells_[icells][i]-1;
			cc_dummy_x = cc_dummy_x + vertices_[icell][0];
			cc_dummy_y = cc_dummy_y + vertices_[icell][1];
			cc_dummy_z = cc_dummy_z + vertices_[icell][2];
		}
		cell_centers[icells][0] = cc_dummy_x/vertices_per_cell;
		cell_centers[icells][1] = cc_dummy_y/vertices_per_cell;
		cell_centers[icells][2] = cc_dummy_z/vertices_per_cell;
		cc_dummy_x = 0; cc_dummy_y = 0; cc_dummy_z =0;

		//std::cout << "** cell_centers["<< icells <<"]:" << cell_centers[icells][0] << " " << cell_centers[icells][1] << " " << cell_centers[icells][2] << "\n";
	}
	
	check_material_ID = material_ID;
	sort( check_material_ID.begin(), check_material_ID.end() );
	check_material_ID.erase( unique( check_material_ID.begin(), check_material_ID.end() ), check_material_ID.end() );
	std::cout << "check_material_ID.size() = " << check_material_ID.size() << std::endl;

	/*
	for (int i = 0; i!=check_material_ID.size(); i++){
		std::cout << "check_material_ID[" << patch::to_string(i) <<"] = " << check_material_ID[i] << std::endl;
	}
	*/

	if (check_material_ID.size() > num_fracs) {
		std::cerr << "Error: Incompatible between the numbers of fractures from the final *.inp file and the data input file." << std::endl;
		std::cerr << "Please correct the data input file \"" << inputfile_ << "\"." << std::endl;
		return -1;
	}
	
	/*
	for (int icells=0; icells < num_cells; icells++)
	{
		std::cout << "** material_ID["<< icells <<"]:" << material_ID[icells] << "\n";
		std::cout << "** cell_centers["<< icells <<"]:" << cell_centers[icells][0] << " " << cell_centers[icells][1] << " " << cell_centers[icells][2] << "\n";
	}
	*/
	
	infile.close();

	// Filter element cells according to its material_ID
	
	for (int icells=0; icells < num_cells; icells++)
	{		
		//std::cout << "material_ID["<< icells << "] = " << material_ID[icells]<< std::endl;
		Global_cellID[ material_ID[icells] -1 ].push_back(icells) ;
		Points_list_withID[ material_ID[icells] -1 ].push_back(cell_centers[icells]);
	}

	for (int ifrac=0; ifrac < num_fracs; ifrac++)
	{
		std::string  s_data_outfile("Vertices_of_Frac_" + patch::to_string(ifrac) + ".dat");
		std::ofstream  data_outfile(s_data_outfile.c_str());
		
		if (Points_list_withID[ifrac].size() > 0){
			std::cout << "***\n";
			//std::cout << "Global_cellID["<< ifrac <<"].size() = " << Global_cellID[ifrac].size() << std::endl;
			std::cout << "Points_list_withID["<< ifrac <<"].size() = " << Points_list_withID[ifrac].size() << std::endl;

			std::list<double>::iterator E_radius = E_radius_list.begin();
			std::advance(E_radius, ifrac);
			std::cout << "E_radius["<< ifrac <<"].size() = " << *E_radius << std::endl;

			//--------------------------------------------------------------------------------
			// Rotation fracture i to the xy plane
			//--------------------------------------------------------------------------------
			//std::vector<double>	frac_normal_vect_dum(3);
			//std::vector<double>	frac_center_dum(3);

			std::list< std::vector<double> >::iterator frac_center_dum(E_center_list.begin());
			std::advance(frac_center_dum,ifrac);

			std::list< std::vector<double> >::iterator frac_normal_vect_dum(E_normalvect_list.begin());
			std::advance(frac_normal_vect_dum,ifrac);

			std::vector< std::vector<double> > arr_toZero_pts_frac_dummy (Points_list_withID[ifrac].size(), std::vector<double>(3));
			std::vector< std::vector<double> > arr_rot_pts_frac_dummy (Points_list_withID[ifrac].size(), std::vector<double>(3));


			for (int i =0; i != Points_list_withID[ifrac].size(); i++){
				std::list< std::vector<double> >::iterator iter_pts_list(Points_list_withID[ifrac].begin());
				std::advance(iter_pts_list, i);

				arr_toZero_pts_frac_dummy[i][0] = (*iter_pts_list)[0] -(*frac_center_dum)[0];
				arr_toZero_pts_frac_dummy[i][1] = (*iter_pts_list)[1] -(*frac_center_dum)[1];
				arr_toZero_pts_frac_dummy[i][2] = (*iter_pts_list)[2] -(*frac_center_dum)[2];
			
				//std::cout << "(*iter_pts_list) = " << (*iter_pts_list)[0] << " " << (*iter_pts_list)[1] << " " << (*iter_pts_list)[2] << std::endl;
				//std::cout << "(*frac_center_dum) = " << (*frac_center_dum)[0] << " " << (*frac_center_dum)[1] << " " << (*frac_center_dum)[2] << std::endl;
				//std::cout << "arr_toZero_pts_frac_dummy = " << arr_toZero_pts_frac_dummy[i][0] << " " << arr_toZero_pts_frac_dummy[i][1] << " " << arr_toZero_pts_frac_dummy[i][2] << std::endl;
					
			}
			std::vector<double> e3(3);
			e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;

			//std::cout << "(*frac_normal_vect_dum) = " << (*frac_normal_vect_dum)[0] << " " << (*frac_normal_vect_dum)[1] << " " << (*frac_normal_vect_dum)[2] << std::endl;
			//std::cout << "arr_toZero_pts_frac_dummy = " << (arr_toZero_pts_frac_dummy)[ifrac][0] << " " << (arr_toZero_pts_frac_dummy)[ifrac][1] << " " << (arr_toZero_pts_frac_dummy)[ifrac][2] << std::endl;

			if( pow((*frac_normal_vect_dum)[0] - e3[0],2) + pow((*frac_normal_vect_dum)[1] - e3[1],2) + pow((*frac_normal_vect_dum)[2] - e3[2],2) < 1E-16){
				std::cout << "Fracture has already been in XY plane" << std::endl;
				arr_rot_pts_frac_dummy = arr_toZero_pts_frac_dummy;
			}
			else{

				std::vector<double> rotvect(3), rotvect_norm(3);
				rotvect = MyCrossProduct<double>((*frac_normal_vect_dum), e3);
				rotvect_norm = MyNormalize<double>(rotvect);

				//std::cout << "(rotvect_norm) = " << (rotvect_norm)[0] << " " << (rotvect_norm)[1] << " " << (rotvect_norm)[2] << std::endl;
			
				double theta = MyAngle<double>((*frac_normal_vect_dum),e3);
				//std::cout << "** theta = " << theta *180/M_PI <<std::endl;
				//double theta = MyAngle<double>((*frac_normal_vect_dum),e3);
				//double c = cos(theta);
				//double s = sin(theta);
				setUpRotationMatrix(theta, rotvect_norm);

				for (int i =0; i != Points_list_withID[ifrac].size(); i++){
					arr_rot_pts_frac_dummy[i][0] = 0.0;
					arr_rot_pts_frac_dummy[i][1] = 0.0;
					arr_rot_pts_frac_dummy[i][2] = 0.0;
					
					inputMatrix[0][0] = arr_toZero_pts_frac_dummy[i][0];
					inputMatrix[1][0] = arr_toZero_pts_frac_dummy[i][1];
					inputMatrix[2][0] = arr_toZero_pts_frac_dummy[i][2];
					inputMatrix[3][0] = 1.0;
					
					//std::cout << "arr_toZero_pts_frac_dummy = " << arr_toZero_pts_frac_dummy[i][0] << " " << arr_toZero_pts_frac_dummy[i][1] << " " << arr_toZero_pts_frac_dummy[i][2] << std::endl;
					//std::cout << "inputMatrix = " << (inputMatrix)[0][0] << " " << (inputMatrix)[1][0] << " " << (inputMatrix)[2][0] << std::endl;
					
					multiplyMatrix();

					arr_rot_pts_frac_dummy[i][0] = outputMatrix[0][0];
					arr_rot_pts_frac_dummy[i][1] = outputMatrix[0][1];
					arr_rot_pts_frac_dummy[i][2] = outputMatrix[0][2];
					
					outputMatrix[0][0] = 0.0;
					outputMatrix[0][1] = 0.0;
					outputMatrix[0][2] = 0.0;
					outputMatrix[0][3] = 0.0;
				}
			}

			// Write data of the rotated fracture to text file for generating random field using R packages
			/*
			for (int ipts = 0; ipts!= Points_list_withID[ifrac].size(); ipts++){
				std::list<std::vector<double> >::iterator iter_pts(Points_list_withID[ifrac].begin());
				std::advance(iter_pts, ipts);
				data_outfile << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) <<
								(*iter_pts)[0] << " " << (*iter_pts)[1] << " " << (*iter_pts)[0] << std::endl;
			}
			 */

			for (int ipts = 0; ipts!= arr_rot_pts_frac_dummy.size(); ipts++){
				data_outfile << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) <<
						arr_rot_pts_frac_dummy[ipts][0] << " " << arr_rot_pts_frac_dummy[ipts][1] << " " << arr_rot_pts_frac_dummy[ipts][2] << std::endl;
			}
		}
	}
	
	/*
	for (int icells=0; icells < num_cells; icells++)
	{		
		std::cout << "Arr_Global_cellID[ " << icells << " ] = " <<  Arr_Global_cellID[ icells ] << std::endl;
	}
	*/
	//=======================================================================================================
	// Generate Gaussian's random fields of apertures using R-packages (gstat,fields)
	//=======================================================================================================
	//std::system("rm -rf R_scripts");
	//std::system("mkdir R_scripts");

	std::cout << "Generating aperture data... "<< std::endl;

	// Create scripts to run R packages
	for (int ifrac=0; ifrac < num_fracs; ifrac++)
	{
		if (Points_list_withID[ifrac].size() > 0){
			std::list<double>::iterator frac_radius_dum(E_radius_list.begin());
			std::advance(frac_radius_dum,ifrac);

			std::string  outfile_Rscr_("generatingSCRF_Frac_" + patch::to_string(ifrac) + ".R");
			std::ofstream  outfile_Rscr(outfile_Rscr_.c_str());

			outfile_Rscr << "# Generating spatially correlated random fields with R\n";
			outfile_Rscr << "# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/\n";

			outfile_Rscr << "# unconditional simulations on a 100 x 100 grid using gstat\n";
			if(bool_linux){
				outfile_Rscr << "library(gstat, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";
				outfile_Rscr << "library(sp, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";
				outfile_Rscr << "library(maptools, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";
				outfile_Rscr << "library(spam, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";
				outfile_Rscr << "library(maps, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";
				outfile_Rscr << "library(fields, lib.loc='/work/irlin104_1/ngotr/Outils/R/InstalledPackages')\n";			
			}
			else{				
				outfile_Rscr << "library(gstat)\n";
				outfile_Rscr << "library(sp)\n";
				outfile_Rscr << "library(maptools)\n";
				outfile_Rscr << "library(spam)\n";
				outfile_Rscr << "library(maps)\n";
				outfile_Rscr << "library(fields)\n";
			}
			outfile_Rscr << "\n";

			outfile_Rscr << "# set current working directory\n";
			outfile_Rscr << "setwd(\""<< s_cwd << "\")\n";

			outfile_Rscr << "\n";

			outfile_Rscr << "# Read data\n";
			outfile_Rscr << "data_xyz = read.table('Vertices_of_Frac_" << patch::to_string(ifrac) << ".dat', header = FALSE)\n";
			outfile_Rscr << "names(data_xyz) <- c(\"x\",\"y\",\"z\")\n";

			outfile_Rscr << "xy = data_xyz[,1:2]\n";
			outfile_Rscr << "names(xy) <- c(\"x\",\"y\")\n";

			outfile_Rscr << "# define the gstat object (spatial model)\n";
			outfile_Rscr << "g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,"
						<< "model=vgm(psill=0.025,model=\"Exp\",range=" << (*frac_radius_dum) * 0.7 << "), nmax=30)\n";

			outfile_Rscr << "# make simulation based on the stat object\n";
			outfile_Rscr << "yy <- predict(g.dummy, newdata=xy, nsim=1)\n";

			outfile_Rscr << "scale = " << aperture[ifrac] << "\n";

			std::string infile_apert_gen_("R_aperture_Frac_" + patch::to_string(ifrac) + ".txt");

			outfile_Rscr << "write.table(yy$sim1 * scale, file = \"" << infile_apert_gen_ <<"\", quote = FALSE, col.names = FALSE, row.names = FALSE)\n";

			outfile_Rscr.close();

			std::string cmd("R CMD BATCH " + outfile_Rscr_);
			std::system(cmd.c_str());

			std::vector<double> list_aper_dummy;

			std::ifstream infile_apert_gen(infile_apert_gen_.c_str());
			for (int ipts = 0; ipts!= Points_list_withID[ifrac].size(); ipts++){
				std::getline(infile_apert_gen, line);
				std::istringstream ist(line);
				double db_dummy_1;
				ist >> db_dummy_1;
				Aper_list_withID[ifrac].push_back(db_dummy_1);
			}
			//std::cout << "Aper_list_withID[" << ifrac <<"][" << ipts <<"= " << Aper_list_withID[ifrac][ipts]<< std::endl;
		}
	}
	std::cout << "Generating aperture data... Complete. \n";

	int cell_globalID = 0;
	for (int ifrac=0; ifrac < num_fracs; ifrac++)
	{
		if (Points_list_withID[ifrac].size() > 0){
			for (int ipts = 0; ipts!= Points_list_withID[ifrac].size(); ipts++){
				std::list<int>::iterator iter_GlobalID(Global_cellID[ifrac].begin());
				std::advance(iter_GlobalID,ipts);

				std::list<double>::iterator iter_Aper_list_withID(Aper_list_withID[ifrac].begin());
				std::advance(iter_Aper_list_withID,ipts);

				Arr_Global_cellID[cell_globalID] = *iter_GlobalID;
				Arr_Aper_Hete[cell_globalID] = *iter_Aper_list_withID;

				cell_globalID++;
			}
		}
	}

	// Sort array of aperture values based on global index
	std::vector<Apert_ID> arr_Apert_ID(num_cells);
	for (int i =0; i!= num_cells; i++){
		(arr_Apert_ID[i]).Apert_ =  Arr_Aper_Hete[i];
		(arr_Apert_ID[i]).ID_ =  Arr_Global_cellID[i];
	}

	std::sort(arr_Apert_ID.begin(), arr_Apert_ID.end(), compareID);

	/*
	for (int icell=0; icell < num_cells; icell++){
		std::cout << "Cell number " << Arr_Global_cellID[icell] << " - Aperture = " << Arr_Aper_Hete[icell] << " (m).\n";
	}
	*/

	// Clean the directory
	if(1){
		std::string cmd("");
		cmd = "rm -rf R_aperture_Frac_*.txt";
		std::system(cmd.c_str());

		cmd = "rm -rf generatingSCRF_Frac_*";
		std::system(cmd.c_str());

		cmd = "rm -rf Vertices_of_Frac_*.dat";
		std::system(cmd.c_str());
	}

	// ===========================================================================
	// Write dgf file
	// ===========================================================================

	std::string  outfile_(s_cwd + "/" + filename_.substr(0,filename_.length()-4)+ ".dgf");
	find_and_replace( outfile_, "//", "/");
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
			if (fabs(vertices_[it][j]) < eps){
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
		// Frac-to-frac heterogeneity
		if (hete == 1){
			int mat_ID(material_ID[it]-1);

			//outfile << vertices_per_cell  << " ";
			std::list<double>::iterator perm_i_X(permeability_X.begin());
			std::advance(perm_i_X, mat_ID);

			std::list<double>::iterator perm_i_Y(permeability_Y.begin());
			std::advance(perm_i_Y, mat_ID);

			std::list<double>::iterator perm_i_Z(permeability_Z.begin());
			std::advance(perm_i_Z, mat_ID);

			for (i=0; i < vertices_per_cell; i++)
				outfile << cells_[it][i] - 1 << " ";
			outfile << *perm_i_X << " " << *perm_i_Y << " " << *perm_i_Z << " ";

			outfile << aperture[mat_ID] << std::endl;
		}
		// Cell-to-cell heterogeneity
		else if (hete == 2){
			for (i=0; i < vertices_per_cell; i++)
				outfile << cells_[it][i] - 1 << " ";
			outfile << compute_K_from_apert_mul((arr_Apert_ID[it]).Apert_) << " ";
			outfile << compute_K_from_apert_mul((arr_Apert_ID[it]).Apert_) << " ";
			outfile << compute_K_from_apert_mul((arr_Apert_ID[it]).Apert_) << " ";
			outfile << (arr_Apert_ID[it]).Apert_ << std::endl;
		}
	}

	outfile << "#\n";
	//std::cout << "material_ID[0] = " << material_ID[0] << std::endl;

	outfile.close();

	// ===========================================================================
	// Write vtk file
	// ===========================================================================

	std::string  outfile_vtk_(filename_.substr(0,filename_.length()-4)+ "_hete.vtk");
	std::ofstream  outfile_vtk; outfile_vtk.open(outfile_vtk_.c_str());

	// header
	outfile_vtk << "# vtk DataFile Version 4.0\n";
	outfile_vtk << "vtk output\n";
	outfile_vtk << "ASCII\n";
	outfile_vtk << "DATASET UNSTRUCTURED_GRID\n";
	outfile_vtk << "POINTS " << num_vertices << " float\n";

	// Write vertices
	// Sometime the rotation causes Points to be eps close of zero, set those Points to 0.
	i=0;
	for (int it=0; it != num_vertices; ++it)
	{
		for (int j = 0; j!=3; j++){
			if (fabs( vertices_[it][j]) < eps){
				vertices_[it][j] = 0;
			}
		}

		outfile_vtk  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1]  << " " << vertices_[it][2]  << std::endl;
		//outfile_vtk  << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20) << vertices_[it][0] << " " << vertices_[it][1]  << std::endl;
	}

	// Write cells
	outfile_vtk << "CELLS " << num_cells << " " << num_cells * (vertices_per_cell+1) << std::endl;

	for (int it=0; it != num_cells; ++it)
	{
		outfile_vtk << vertices_per_cell  << " ";
		for (i=0; i < vertices_per_cell; i++)
			outfile_vtk << cells_[it][i] - 1 << " ";
		outfile_vtk << std::endl;
	}

	// cell types (type 5 is a three-node triangle)

	outfile_vtk << "\n";
	outfile_vtk << "CELL_TYPES " << num_cells << std::endl;
	if (vertices_per_cell == 3){
		for (int j = 0; j < num_cells; ++j)
		{
			outfile_vtk << "5" << std::endl;
		}
	}
	else{
		for (int j = 0; j < num_cells; ++j)
		{
			outfile_vtk << "7" << std::endl;
		}
	}

	// Write material ID
	outfile_vtk << "CELL_DATA     " << num_cells << std::endl;
	outfile_vtk << "SCALARS Material_ID int 1\n";
	outfile_vtk << "LOOKUP_TABLE	default\n";

	//outfile_vtk << "Material_ID " << num_vertices << "  int\n";
	
	for (int it =0; it != material_ID.size(); ++it)
	{
		outfile_vtk  << "  " << material_ID[it] << std::endl;
	}

	// Write aperture
	//outfile_vtk << "CELL_DATA     " << num_cells << std::endl;
	outfile_vtk << "SCALARS Aperture float\n";
	outfile_vtk << "LOOKUP_TABLE	default\n";

	//outfile_vtk << "Material_ID " << num_vertices << "  int\n";

	for (int it =0; it != arr_Apert_ID.size(); ++it)
	{
		//outfile_vtk  << "  " << (arr_Apert_ID[it]).Apert_ << std::endl;
		outfile_vtk  << (arr_Apert_ID[it]).Apert_ << std::endl;
	}

	outfile_vtk << "SCALARS Permeability float\n";
	outfile_vtk << "LOOKUP_TABLE	default\n";

	for (int it =0; it != arr_Apert_ID.size(); ++it)
	{
		if (hete == 1){
			int mat_ID(material_ID[it]-1);

			std::list<double>::iterator perm_i_X(permeability_X.begin());
			std::advance(perm_i_X, mat_ID);
			outfile_vtk << *perm_i_X << std::endl;
		}
		// Cell-to-cell heterogeneity
		else if (hete == 2){
			outfile_vtk << compute_K_from_apert_mul((arr_Apert_ID[it]).Apert_) << std::endl;
		}
	}
	outfile_vtk.close();

	std::cout << "Aperture = " << 1E-3 << ", the corresponding permeability K = " <<  compute_K_from_apert(1E-3) << std::endl;
	std::cout << "Aperture = " << 1E-3 << ", the corresponding transmissivity sigma = " << compute_sigma_from_apert(1E-3) << std::endl;
	std::cout << "Aperture = " << 1E-3 << ", the corresponding permeability K_mul = " <<  compute_K_from_apert_mul(1E-3) << std::endl;
	std::cout << "Aperture = " << 1E-3 << ", the corresponding transmissivity sigma_mul = " << compute_sigma_from_apert_mul(1E-3) << std::endl;

	std::cout << std::endl <<"Program terminated with success." << std::endl << std::endl;

	return 0;
}

