/*
 * read_fracs_params_file.hh
 *
 *  Created on: Sep 5, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_READ_FRACS_PARAMS_FILE_HH_
#define INCLUDE_READ_FRACS_PARAMS_FILE_HH_

#include <list>
#include <algorithm>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>
#include "Vector_3_Operations.hh"

#ifndef VERY_VERBOSE
#define VERY_VERBOSE 0
#endif

using namespace std;

namespace CGAL{

template < typename K, typename NT, typename Ell_list_>
Ell_list_ read_fracs_params_file_ell(const std::string &filename_)
{
	typedef typename CGAL::Point_3<K>	Point_3;
	typedef typename CGAL::Vector_3<K>	Vector_3;
	Ell_list_ Ell_list;

	int 	E_name(0), num_frac(0), num_ell(0);								// Name of the ellipse
	float 	E_center_0, E_center_1, E_center_2 /*Center of the ellipse*/,
	E_xradius, E_yradius, /* x_radius, y_radius */
	E_normalvect_0, E_normalvect_1, E_normalvect_2 /* normal vector */;

	std::list<int> 			E_name_list;						// Name of the ellipse
	std::list<int> 			E_family_list;						// Family of the ellipse
	std::list<Point_3> 		E_center_list;						// Center of the ellipse
	std::list<NT> 			E_xradius_list, E_yradius_list;		// x_radius, y_radius
	std::list<Vector_3> 	E_normalvect_list;					// normal vector

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename_ << "\"."<< std::endl;
		exit(-1);
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
		ist >> num_frac;
	}

	//std::cout << "Number of fractures = " << num_frac << std::endl;

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	//for (int it =0; it!=num_frac; it++){
	int num_count_frac(0);

	//while (!infile.eof()){
	while(num_count_frac < num_frac && !infile.eof()){
		// Read fracture name and family
		int name, family, numPts(0);

		std::istringstream ist_nf(line);
		ist_nf.str(line);
		ist_nf.clear();
		ist_nf >> name >> family;


		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
		}

		if (family == 1){
			num_ell++;
			E_name = name;
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> 	E_center_0 		>> E_center_1 		>> E_center_2 >>
			E_xradius	 	>> E_yradius    	>>
			E_normalvect_0 	>> E_normalvect_1 	>> E_normalvect_2;

			E_name_list.push_back(E_name);
			E_center_list.push_back(Point_3(E_center_0, E_center_1, E_center_2));
			E_xradius_list.push_back(E_xradius);
			E_yradius_list.push_back(E_yradius);
			E_normalvect_list.push_back(Vector_3(E_normalvect_0, E_normalvect_1, E_normalvect_2));

			if (VERY_VERBOSE == 1){
				std::cout << E_center_0	<< " " << E_center_1 << " " << E_center_2 <<std::endl;
			}
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			if (infile.eof()){
				break;
			}
		}

		num_count_frac++;
	}

	if (VERY_VERBOSE == 1){
		for (std::list<int>::iterator iter = E_name_list.begin(); iter!= E_name_list.end(); iter++){
			std::cout << "E_name_list : " << *iter << std::endl;
		}
	}
	Ell_list.num_Ell		= num_ell;
	Ell_list.E_name 		= E_name_list;
	Ell_list.E_family		= E_family_list;
	Ell_list.E_center		= E_center_list;
	Ell_list.E_xradius		= E_xradius_list;
	Ell_list.E_yradius		= E_yradius_list;
	Ell_list.E_normalvect	= E_normalvect_list;

	return Ell_list;
	//return true;
}


template < typename K, typename NT, typename Ell_list_>
Ell_list_ read_fracs_params_file_ell_desactfracs(const std::string &filename_, const std::string &desactfracsfilename_)
{
	typedef typename CGAL::Point_3<K>	Point_3;
	typedef typename CGAL::Vector_3<K>	Vector_3;
	Ell_list_ Ell_list;

	int 	E_name(0), num_frac(0), num_ell(0);								// Name of the ellipse
	float 	E_center_0, E_center_1, E_center_2 /*Center of the ellipse*/,
	E_xradius, E_yradius, /* x_radius, y_radius */
	E_normalvect_0, E_normalvect_1, E_normalvect_2 /* normal vector */;

	std::list<int> 			desactfracs;		// Name of the ellipse
	std::list<int> 			E_name_list;						// Name of the ellipse
	std::list<int> 			E_family_list;						// Family of the ellipse
	std::list<Point_3> 		E_center_list;						// Center of the ellipse
	std::list<NT> 			E_xradius_list, E_yradius_list;		// x_radius, y_radius
	std::list<Vector_3> 	E_normalvect_list;					// normal vector


	std::ifstream    data_infile(desactfracsfilename_.c_str());

	std::string line, linenospace;
	std::istringstream ist(line);

	if (! data_infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< desactfracsfilename_ << "\"."<< std::endl;
		//exit(-1);
	}
	else{
		while (!data_infile.eof()){
			std::getline(data_infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			int dummy_int;
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> dummy_int;
			if (linenospace.size() > 0){
				desactfracs.push_back(dummy_int);
			}
		}
	}

	//std::cout << "desactfracs.size() = "<< desactfracs.size() << std::endl;

	if (desactfracs.size() > 0){
		desactfracs.sort();
		desactfracs.unique();
	}

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename_ << "\"."<< std::endl;
		exit(-1);
	}

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_frac;
	}

	//std::cout << "Number of fractures = " << num_frac << std::endl;

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	//for (int it =0; it!=num_frac; it++){
	int num_count_frac(0);

	//while (!infile.eof()){
	while(num_count_frac < num_frac && !infile.eof()){
		// Read fracture name and family
		int name, family, numPts(0);

		std::istringstream ist_nf(line);
		ist_nf.str(line);
		ist_nf.clear();
		ist_nf >> name >> family;


		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
		}

		if (family == 1){
			num_ell++;

			bool found(false);
			if(desactfracs.size() > 0){
				found = (std::find(desactfracs.begin(), desactfracs.end(), name) != desactfracs.end());
			}

			if (found){
				num_ell = num_ell -1;
			}
			else{
				E_name = name;
				std::istringstream ist(line);
				ist.str(line);
				ist.clear();
				ist >> 	E_center_0 		>> E_center_1 		>> E_center_2 >>
				E_xradius	 	>> E_yradius    	>>
				E_normalvect_0 	>> E_normalvect_1 	>> E_normalvect_2;

				E_name_list.push_back(E_name);
				E_center_list.push_back(Point_3(E_center_0, E_center_1, E_center_2));
				E_xradius_list.push_back(E_xradius);
				E_yradius_list.push_back(E_yradius);
				E_normalvect_list.push_back(Vector_3(E_normalvect_0, E_normalvect_1, E_normalvect_2));

				if (VERY_VERBOSE == 1){
					std::cout << E_center_0	<< " " << E_center_1 << " " << E_center_2 <<std::endl;
				}
			}
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			if (infile.eof()){
				break;
			}
		}

		num_count_frac++;
	}

	if (VERY_VERBOSE == 1){
		for (std::list<int>::iterator iter = E_name_list.begin(); iter!= E_name_list.end(); iter++){
			std::cout << "E_name_list : " << *iter << std::endl;
		}
	}
	Ell_list.num_Ell		= num_ell;
	Ell_list.E_name 		= E_name_list;
	Ell_list.E_family		= E_family_list;
	Ell_list.E_center		= E_center_list;
	Ell_list.E_xradius		= E_xradius_list;
	Ell_list.E_yradius		= E_yradius_list;
	Ell_list.E_normalvect	= E_normalvect_list;

	return Ell_list;
	//return true;
}

template < typename K, typename NT, typename Poly_list_>
Poly_list_ read_fracs_params_file_poly(const std::string &filename_)
{
	typedef typename CGAL::Point_3<K>	Point_3;
	typedef typename CGAL::Plane_3<K>	Plane_3;
	typedef typename CGAL::Vector_3<K>	Vector_3;
	typedef typename CGAL::Line_3<K>	Line_3;
	typedef typename K::Intersect_3 	Intersect_3;

	typedef typename std::list<Point_3>::iterator	Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator	Vect_Iterator;

	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type	result_type_lp;

	Poly_list_ Poly_list;
	int num_Poly(0);		// Number of polygons
	int P_name(0), num_frac(0), num_poly(0);				// Name of the ellipse

	std::list<int> 					P_name_list;		// Name of the ellipse
	std::list<int> 					P_family_list;		// Family of the ellipse
	std::list< std::list<Point_3> > P_point_list;		// List of points
	std::list<Vector_3> 			P_normalvect_list;	// normal vector

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename_ << "\"."<< std::endl;
		exit(-1);
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
		ist >> num_frac;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	//for (int it =0; it!=num_frac; it++){
	int num_count_frac(0);

	//while (!infile.eof()){
	while(num_count_frac < num_frac && !infile.eof()){
		// Read fracture name and family
		int name, family, numPts(0);
		std::list<Point_3> points, new_points;

		std::istringstream ist_nf(line);
		ist_nf.str(line);
		ist_nf.clear();
		ist_nf >> name >> family;

		if (family != 1){
			ist_nf >> numPts;
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
		}

		if (family != 1){
			num_poly++;
			P_name = name;

			for (int it = 0; it < numPts; it++){

				std::istringstream ist(line);
				ist.str(line);
				ist.clear();

				float dummy_1, dummy_2, dummy_3;
				ist >> dummy_1 >> dummy_2 >> dummy_3;
				points.push_back(Point_3(dummy_1, dummy_2, dummy_3));

				std::getline(infile, line);
			}

			/*
			 *  Sometime, the fourth point p4 does not lie on the plane (P) of three first ones (p1,p2,p3)
			 *  We must project p4 to (P)
			 */
			if (points.size() == 4){
				Plane_3 pl_poly;
				Line_3 pp_line;

				for (int i =0; i!=3;i++){
					Pts_Iterator it = points.begin();
					std::advance(it,i);
					new_points.push_back(*it);
				}
				// Plane (P) of the fracture
				pl_poly = Plane_3(*new_points.begin(), *(boost::next(points.begin())), *(boost::next(boost::next(points.begin()))));

				// Perpendicular line to and passing p4
				pp_line = pl_poly.perpendicular_line(*(boost::next(boost::next(boost::next(points.begin())))));

				result_type_lp result_lp = CGAL::intersection(pp_line, pl_poly);

				if (result_lp) {
					if (const Point_3* p = boost::get<Point_3 >(&*result_lp)) {
						//cout << "Fourth point: " << *p << endl;
						new_points.push_back(*p);
					}
				}
				points.clear();
				points = new_points;
				new_points.clear();
			}
			/*
			for (Pts_Iterator it = points.begin(); it != points.end();it++){
				cout <<  "points = " << *it << endl;
			}
			*/

			Point_3 P_c_cumul, P_c;
			for (Pts_Iterator it = points.begin(); it != points.end(); it++){
				P_c_cumul = Point_3(P_c_cumul.x() + (*it).x(), P_c_cumul.y() + (*it).y(), P_c_cumul.z() + (*it).z() );
			}

			NT ratio = 1.0 / points.size();
			P_c = Point_3(P_c_cumul.x()*ratio, P_c_cumul.y()*ratio,P_c_cumul.z()*ratio);

			Pts_Iterator P_A = points.begin();
			Pts_Iterator P_B = points.begin();
			std::advance(P_B, 1);

			if (VERY_VERBOSE == 1){
				std::cout << "points.size() = " << points.size() << std::endl;
				std::cout << "ratio = " << ratio <<std::endl;
				std::cout << "ratio = " << 1 / CGAL::to_double(points.size()) <<std::endl;

				std::cout << "P_c = " << P_c <<std::endl;
				std::cout << "P_A = " << *P_A <<std::endl;
				std::cout << "P_B = " << *P_B <<std::endl;
			}

			Vector_3 vect_A(Vector_3(P_c, *P_A) ), vect_B(Vector_3(P_c, *P_B) ), normal_vector;

			normal_vector = CGAL::Normalize<K,NT>(CGAL::MyCross<K>(vect_A,vect_B));

			P_name_list.push_back(name);
			P_point_list.push_back(points);
			P_family_list.push_back(family);
			P_normalvect_list.push_back(normal_vector);
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			if (infile.eof()){
				break;
			}
		}
		points.clear();

		num_count_frac++;
	}


	if (VERY_VERBOSE == 1){
		for (std::list<int>::iterator iter = P_name_list.begin(); iter!= P_name_list.end(); iter++){
			std::cout << "P_name_list : " << *iter << std::endl;
		}

		for (Vect_Iterator iter = P_normalvect_list.begin(); iter!= P_normalvect_list.end(); iter++){
			std::cout << "P_normalvect_list : " << *iter << std::endl;
		}
	}

	Poly_list.num_Poly		= num_poly;
	Poly_list.P_name 		= P_name_list;
	Poly_list.P_family		= P_family_list;
	Poly_list.P_points		= P_point_list;
	Poly_list.P_normalvect	= P_normalvect_list;

	return Poly_list;
}

template < typename K, typename NT, typename Poly_list_>
Poly_list_ read_fracs_params_file_poly_insidebbox(const std::string &filename_, const CGAL::Point_3<K> p_min, const CGAL::Point_3<K> p_max)
{
	typedef typename CGAL::Point_3<K>	Point_3;
	typedef typename CGAL::Plane_3<K>	Plane_3;
	typedef typename CGAL::Vector_3<K>	Vector_3;
	typedef typename CGAL::Line_3<K>	Line_3;
	typedef typename K::Intersect_3 	Intersect_3;

	typedef typename std::list<Point_3>::iterator	Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator	Vect_Iterator;

	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type	result_type_lp;

	Poly_list_ Poly_list;
	int num_Poly(0);		// Number of polygons
	int P_name(0), num_frac(0), num_poly(0);				// Name of the ellipse

	std::list<int> 					P_name_list;		// Name of the ellipse
	std::list<int> 					P_family_list;		// Family of the ellipse
	std::list< std::list<Point_3> > P_point_list;		// List of points
	std::list<Vector_3> 			P_normalvect_list;	// normal vector

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename_ << "\"."<< std::endl;
		exit(-1);
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
		ist >> num_frac;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	//for (int it =0; it!=num_frac; it++){
	int num_count_frac(0);

	//while (!infile.eof()){
	while(num_count_frac < num_frac && !infile.eof()){
		// Read fracture name and family
		int name, family, numPts(0);
		std::list<Point_3> points, new_points;

		std::istringstream ist_nf(line);
		ist_nf.str(line);
		ist_nf.clear();
		ist_nf >> name >> family;

		if (family != 1){
			ist_nf >> numPts;
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
		}

		if (family != 1){
			num_poly++;
			P_name = name;

			for (int it = 0; it < numPts; it++){

				std::istringstream ist(line);
				ist.str(line);
				ist.clear();

				float dummy_1, dummy_2, dummy_3;
				ist >> dummy_1 >> dummy_2 >> dummy_3;
				points.push_back(Point_3(dummy_1, dummy_2, dummy_3));

				std::getline(infile, line);
			}

			/*
			 *  Sometime, the fourth point p4 does not lie on the plane (P) of three first ones (p1,p2,p3)
			 *  We must project p4 to (P)
			 */
			if (points.size() == 4){
				Plane_3 pl_poly;
				Line_3 pp_line;

				for (int i =0; i!=3;i++){
					Pts_Iterator it = points.begin();
					std::advance(it,i);
					new_points.push_back(*it);
				}
				// Plane (P) of the fracture
				pl_poly = Plane_3(*new_points.begin(), *(boost::next(points.begin())), *(boost::next(boost::next(points.begin()))));

				// Perpendicular line to and passing p4
				pp_line = pl_poly.perpendicular_line(*(boost::next(boost::next(boost::next(points.begin())))));

				result_type_lp result_lp = CGAL::intersection(pp_line, pl_poly);

				if (result_lp) {
					if (const Point_3* p = boost::get<Point_3 >(&*result_lp)) {
						//cout << "Fourth point: " << *p << endl;
						new_points.push_back(*p);
					}
				}
				points.clear();
				points = new_points;
				new_points.clear();
			}
			/*
			for (Pts_Iterator it = points.begin(); it != points.end();it++){
				cout <<  "points = " << *it << endl;
			}
			*/

			Point_3 P_c_cumul, P_c;
			for (Pts_Iterator it = points.begin(); it != points.end(); it++){
				P_c_cumul = Point_3(P_c_cumul.x() + (*it).x(), P_c_cumul.y() + (*it).y(), P_c_cumul.z() + (*it).z() );
			}

			NT ratio = 1.0 / points.size();
			P_c = Point_3(P_c_cumul.x()*ratio, P_c_cumul.y()*ratio,P_c_cumul.z()*ratio);

			Pts_Iterator P_A = points.begin();
			Pts_Iterator P_B = points.begin();
			std::advance(P_B, 1);

			if (VERY_VERBOSE == 1){
				std::cout << "points.size() = " << points.size() << std::endl;
				std::cout << "ratio = " << ratio <<std::endl;
				std::cout << "ratio = " << 1 / CGAL::to_double(points.size()) <<std::endl;

				std::cout << "P_c = " << P_c <<std::endl;
				std::cout << "P_A = " << *P_A <<std::endl;
				std::cout << "P_B = " << *P_B <<std::endl;
			}

			Vector_3 vect_A(Vector_3(P_c, *P_A) ), vect_B(Vector_3(P_c, *P_B) ), normal_vector;
			normal_vector = CGAL::Normalize<K,NT>(CGAL::MyCross<K>(vect_A,vect_B));

			NT xmin((*points.begin())[0]), xmax((*points.begin())[0]), ymin((*points.begin())[1]), ymax((*points.begin())[1]),
			   zmin((*points.begin())[2]), zmax((*points.begin())[2]);

			for (Pts_Iterator it = points.begin(); it != points.end(); it++){
				if (xmin > (*it)[0]){
					xmin = (*it)[0];
				}
				if (xmax < (*it)[0]){
					xmax = (*it)[0];
				}
				if (ymin > (*it)[1]){
					ymin = (*it)[1];
				}
				if (ymax < (*it)[1]){
					ymax = (*it)[1];
				}
				if (zmin > (*it)[2]){
					zmin = (*it)[2];
				}
				if (zmax < (*it)[2]){
					zmax = (*it)[2];
				}
			}

			if ( xmin > p_max[0] || (xmin == p_max[0] && xmin < xmax) || xmax < p_min[0] || (xmax == p_min[0] && xmax > xmin)||
					ymin > p_max[1] || (ymin == p_max[1] && ymin < ymax) || ymax < p_min[1] || (ymax == p_min[1] && ymax > ymin)||
					zmin > p_max[2] || (zmin == p_max[2] && zmin < zmax) || zmax < p_min[2] || (zmax == p_min[2] && zmax > zmin)){
				if (VERY_VERBOSE == 1){
					std::cout<< "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;
					std::cout<< "p_min[0] = " << p_min[0]<< ", p_max[0] = " << p_max[0]<< std::endl;
					std::cout<< "ymin = " << ymin<< ", ymax = " << ymax<< std::endl;
					std::cout<< "p_min[1] = " << p_min[1]<< ", p_max[1] = " << p_max[1]<< std::endl;
					std::cout<< "zmin = " << zmin<< ", zmax = " << zmax<< std::endl;
					std::cout<< "p_min[2] = " << p_min[2]<< ", p_max[2] = " << p_max[2]<< std::endl;
				}

				num_poly = num_poly -1;
			}
			else{
				P_name_list.push_back(name);
				P_point_list.push_back(points);
				P_family_list.push_back(family);
				P_normalvect_list.push_back(normal_vector);
			}
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			if (infile.eof()){
				break;
			}
		}
		points.clear();

		num_count_frac++;
	}


	if (VERY_VERBOSE == 1){
		for (std::list<int>::iterator iter = P_name_list.begin(); iter!= P_name_list.end(); iter++){
			std::cout << "P_name_list : " << *iter << std::endl;
		}

		for (Vect_Iterator iter = P_normalvect_list.begin(); iter!= P_normalvect_list.end(); iter++){
			std::cout << "P_normalvect_list : " << *iter << std::endl;
		}
	}

	Poly_list.num_Poly		= num_poly;
	Poly_list.P_name 		= P_name_list;
	Poly_list.P_family		= P_family_list;
	Poly_list.P_points		= P_point_list;
	Poly_list.P_normalvect	= P_normalvect_list;

	return Poly_list;
}


template < typename K, typename NT, typename Poly_list_>
Poly_list_ read_fracs_params_file_poly_insidebbox_desactfracs(const std::string &filename_, const std::string &desactfracsfilename_, const CGAL::Point_3<K> p_min, const CGAL::Point_3<K> p_max)
{
	typedef typename CGAL::Point_3<K>	Point_3;
	typedef typename CGAL::Plane_3<K>	Plane_3;
	typedef typename CGAL::Vector_3<K>	Vector_3;
	typedef typename CGAL::Line_3<K>	Line_3;
	typedef typename K::Intersect_3 	Intersect_3;

	typedef typename std::list<Point_3>::iterator	Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator	Vect_Iterator;

	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type	result_type_lp;

	Poly_list_ Poly_list;
	int num_Poly(0);		// Number of polygons
	int P_name(0), num_frac(0), num_poly(0);				// Name of the ellipse

	std::list<int> 					desactfracs;		// Name of the ellipse
	std::list<int> 					P_name_list;		// Name of the ellipse
	std::list<int> 					P_family_list;		// Family of the ellipse
	std::list< std::list<Point_3> > P_point_list;		// List of points
	std::list<Vector_3> 			P_normalvect_list;	// normal vector

	std::ifstream    data_infile(desactfracsfilename_.c_str());

	std::string line, linenospace;
	std::istringstream ist(line);

	if (! data_infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< desactfracsfilename_ << "\"."<< std::endl;
		//exit(-1);
	}
	else{
		while (!data_infile.eof()){
			std::getline(data_infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			int dummy_int;
			std::istringstream ist(line);
			ist.str(line);
			ist.clear();
			ist >> dummy_int;
			if (linenospace.size() > 0){
				desactfracs.push_back(dummy_int);
			}
		}
	}

	std::cout << "desactfracs.size() = "<< desactfracs.size() << std::endl;
	//std::cout << "desactfracs.begin() = "<< *desactfracs.begin() << std::endl;

	if (desactfracs.size() > 0){
		desactfracs.sort();
		desactfracs.unique();
	}

	const char* filename = filename_.c_str();
	std::ifstream    infile(filename);

	if (! infile.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename_ << "\"."<< std::endl;
		exit(-1);
	}

	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_frac;
	}

	std::getline(infile, line);
	while (line.size() == 0 || line.substr(0,1) == "#"){
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
	}

	//for (int it =0; it!=num_frac; it++){
	int num_count_frac(0);

	//while (!infile.eof()){
	while(num_count_frac < num_frac && !infile.eof()){
		// Read fracture name and family
		int name, family, numPts(0);
		std::list<Point_3> points, new_points;

		std::istringstream ist_nf(line);
		ist_nf.str(line);
		ist_nf.clear();
		ist_nf >> name >> family;

		if (family != 1){
			ist_nf >> numPts;
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
		}

		if (family != 1){
			num_poly++;
			P_name = name;

			for (int it = 0; it < numPts; it++){

				std::istringstream ist(line);
				ist.str(line);
				ist.clear();

				float dummy_1, dummy_2, dummy_3;
				ist >> dummy_1 >> dummy_2 >> dummy_3;
				points.push_back(Point_3(dummy_1, dummy_2, dummy_3));

				std::getline(infile, line);
			}

			/*
			 *  Sometime, the fourth point p4 does not lie on the plane (P) of three first ones (p1,p2,p3)
			 *  We must project p4 to (P)
			 */
			if (points.size() == 4){
				Plane_3 pl_poly;
				Line_3 pp_line;

				for (int i =0; i!=3;i++){
					Pts_Iterator it = points.begin();
					std::advance(it,i);
					new_points.push_back(*it);
				}
				// Plane (P) of the fracture
				pl_poly = Plane_3(*new_points.begin(), *(boost::next(points.begin())), *(boost::next(boost::next(points.begin()))));

				// Perpendicular line to and passing p4
				pp_line = pl_poly.perpendicular_line(*(boost::next(boost::next(boost::next(points.begin())))));

				result_type_lp result_lp = CGAL::intersection(pp_line, pl_poly);

				if (result_lp) {
					if (const Point_3* p = boost::get<Point_3 >(&*result_lp)) {
						//cout << "Fourth point: " << *p << endl;
						new_points.push_back(*p);
					}
				}
				points.clear();
				points = new_points;
				new_points.clear();
			}
			/*
			for (Pts_Iterator it = points.begin(); it != points.end();it++){
				cout <<  "points = " << *it << endl;
			}
			*/

			Point_3 P_c_cumul, P_c;
			for (Pts_Iterator it = points.begin(); it != points.end(); it++){
				P_c_cumul = Point_3(P_c_cumul.x() + (*it).x(), P_c_cumul.y() + (*it).y(), P_c_cumul.z() + (*it).z() );
			}

			NT ratio = 1.0 / points.size();
			P_c = Point_3(P_c_cumul.x()*ratio, P_c_cumul.y()*ratio,P_c_cumul.z()*ratio);

			Pts_Iterator P_A = points.begin();
			Pts_Iterator P_B = points.begin();
			std::advance(P_B, 1);

			if (VERY_VERBOSE == 1){
				std::cout << "points.size() = " << points.size() << std::endl;
				std::cout << "ratio = " << ratio <<std::endl;
				std::cout << "ratio = " << 1 / CGAL::to_double(points.size()) <<std::endl;

				std::cout << "P_c = " << P_c <<std::endl;
				std::cout << "P_A = " << *P_A <<std::endl;
				std::cout << "P_B = " << *P_B <<std::endl;
			}

			Vector_3 vect_A(Vector_3(P_c, *P_A) ), vect_B(Vector_3(P_c, *P_B) ), normal_vector;
			normal_vector = CGAL::Normalize<K,NT>(CGAL::MyCross<K>(vect_A,vect_B));

			NT xmin((*points.begin())[0]), xmax((*points.begin())[0]), ymin((*points.begin())[1]), ymax((*points.begin())[1]),
			   zmin((*points.begin())[2]), zmax((*points.begin())[2]);

			for (Pts_Iterator it = points.begin(); it != points.end(); it++){
				if (xmin > (*it)[0]){
					xmin = (*it)[0];
				}
				if (xmax < (*it)[0]){
					xmax = (*it)[0];
				}
				if (ymin > (*it)[1]){
					ymin = (*it)[1];
				}
				if (ymax < (*it)[1]){
					ymax = (*it)[1];
				}
				if (zmin > (*it)[2]){
					zmin = (*it)[2];
				}
				if (zmax < (*it)[2]){
					zmax = (*it)[2];
				}
			}

			bool found(false);
			if(desactfracs.size() > 0){
				found = (std::find(desactfracs.begin(), desactfracs.end(), name) != desactfracs.end());
			}

			//std::cout << "name = " << name << ", found = " << found << std::endl;

			if ( xmin >= p_max[0] || xmax <= p_min[0] || ymin >= p_max[1] || ymax <= p_min[1] || zmin >= p_max[2] || zmax <= p_min[2]){
				num_poly = num_poly -1;
			}
			else{
				if (found){
					num_poly = num_poly -1;
				}
				else{
					P_name_list.push_back(name);
					P_point_list.push_back(points);
					P_family_list.push_back(family);
					P_normalvect_list.push_back(normal_vector);
				}
			}
		}

		std::getline(infile, line);
		while (line.size() == 0 || line.substr(0,1) == "#"){
			std::getline(infile, line);
			if (VERY_VERBOSE == 1){
				std::cout << "Line : " << line << std::endl;
			}
			if (infile.eof()){
				break;
			}
		}
		points.clear();

		num_count_frac++;
	}


	if (VERY_VERBOSE == 1){
		for (std::list<int>::iterator iter = P_name_list.begin(); iter!= P_name_list.end(); iter++){
			std::cout << "P_name_list : " << *iter << std::endl;
		}

		for (Vect_Iterator iter = P_normalvect_list.begin(); iter!= P_normalvect_list.end(); iter++){
			std::cout << "P_normalvect_list : " << *iter << std::endl;
		}
	}

	Poly_list.num_Poly		= num_poly;
	Poly_list.P_name 		= P_name_list;
	Poly_list.P_family		= P_family_list;
	Poly_list.P_points		= P_point_list;
	Poly_list.P_normalvect	= P_normalvect_list;

	return Poly_list;
}

} // end of namespace

#endif /* INCLUDE_READ_FRACS_PARAMS_FILE_HH_ */
