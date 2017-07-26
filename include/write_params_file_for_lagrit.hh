/*
 * write_params_file_for_lagrit.hh
 *
 *  Created on: Sep 8, 2016
 *      Author: ngotr
 */

#ifndef WRITE_PARAMS_FILE_FOR_LAGRIT_HH_
#define WRITE_PARAMS_FILE_FOR_LAGRIT_HH_

#include <include/points_are_close.hh>
#include "myglobal_functions.hh"
#include "MyRotation.hh"
#include "Vector_3_Operations.hh"
#include <include/reserve_orientation_of_clockwise_polygon.hh>

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>

namespace CGAL{
template <typename Ell_, typename K, typename NT>
void write_params_file_for_lagrit(
		std::list<Ell_> & ell_list,
		CGAL::Point_3<K> p_min,
		CGAL::Point_3<K> p_max,
		const std::string &filename_){

	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	typedef typename CGAL::Vector_3<K> Vector_3;
	typedef typename CGAL::Point_3<K> Point_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis

	const char* filename = filename_.c_str();
	std::ofstream    outfile(filename);
	Vector_3 vect_;
	NT x1,y1,z1,x2,y2,z2;

	//std::ofstream    outfile;
	//outfile.open(outfile_.c_str());

	// header
	outfile << "# Number of polygons" << std::endl;
	outfile << ell_list.size() << std::endl;
	outfile << "\n";

	outfile << "# npoly ## Theta ## dL0 ## X1 ## Y1 ## Z1 ## X2 ## Y2 ## Z2" << std::endl;

	int inum(1);
	for (Ell_Iterator it = ell_list.begin(); it != ell_list.end(); it++){
		vect_ = (*it).E_normalvect;
		vect_ = CGAL::Normalize<K, NT>(vect_);
		/*  Find Angle of Rotation*/

		NT theta = CGAL::MyAngle<K, NT>(vect_,e3);
		NT c = cos(CGAL::to_double(theta));
		NT s = sin(CGAL::to_double(theta));

		/* Convert angle to degrees*/
		theta = CGAL::Rad2Deg(theta);

		/* Find angle to Rotate into xy plane */
		Vector_3 rotvect = CGAL::MyCross<K>(vect_, e3);
		Vector_3 u = CGAL::Normalize<K,NT>(rotvect);
		/*
		NT dk(CGAL::max(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-12){
			dk = CGAL::min(p_min[0], p_max[0]);
		}
		 */

		NT dk(CGAL::min(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-12){
			dk = CGAL::max(p_min[0], p_max[0]);
		}

		x1 = - dk * u[0]; x2 = dk * u[0];
		y1 = - dk * u[1]; y2 = dk * u[1];
		z1 = - dk * u[2]; z2 = dk * u[2];

		outfile << inum << " " << std::setprecision(18) << -theta << " " << (*it).E_target_edge_length << " " << x1 << " " << y1 << " " << z1
				<< " " << x2 << " " << y2 << " " << z2 << "\n";

		//outfile << inum << " " << vect_[0] << " " <<  vect_[1] << " " << vect_[2] << std::endl;

		inum++;


	}
	outfile.close();
}

template <typename K, typename NT>
CGAL::Surface_mesh<CGAL::Point_3<K> > remove_near_points_sm(CGAL::Surface_mesh<CGAL::Point_3<K> >  & poly_, NT &tol){

	typedef typename CGAL::Point_3<K> 				Point_3;
	typedef typename std::list<Point_3> 			list_point;
	typedef typename std::list<Point_3>::iterator	Pts_Iterator;
	typedef typename CGAL::Surface_mesh<Point_3> 	Surface_mesh;
	typedef typename Surface_mesh::Vertex_index 	vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 	vertex_iterator;
	typedef typename Surface_mesh::Vertex_range		vertex_range;

	Surface_mesh poly_new;
	list_point list_pts;

	vertex_descriptor vh(0);
	Point_3 point_to_add(Point_3(poly_.point(vh).x(),poly_.point(vh).y(), poly_.point(vh).z()) );
	list_pts.push_back(point_to_add);


	unsigned int iter = 1, iter_end = poly_.number_of_vertices();
	for( ; iter < iter_end; ++iter) {
		vertex_descriptor vh_prev(iter-1);
		vertex_descriptor vh(iter);
		if (iter < iter_end-1){
			if (!points_are_close_custom(poly_.point(vh), poly_.point(vh_prev), tol) ){
				list_pts.push_back( Point_3(poly_.point(vh)[0],poly_.point(vh)[1], poly_.point(vh)[2]) );
			}
		}
		else{
			if (!points_are_close_custom(poly_.point(vh), poly_.point(vh_prev), tol) &&
					!points_are_close_custom(poly_.point(vh), *list_pts.begin(), tol) ){
				list_pts.push_back( Point_3(poly_.point(vh)[0],poly_.point(vh)[1], poly_.point(vh)[2]) );
			}
		}
	}

	for (Pts_Iterator it = list_pts.begin(); it != list_pts.end(); it++){
		poly_new.add_vertex(*it);
	}

	vertex_range poly_new_range = poly_new.vertices();
	poly_new.add_face(poly_new_range);

	return poly_new;
}

template <typename Ell_, typename K, typename NT>
void export_inp_files_from_ell_list_xy_plane(std::list<Ell_> & ell_list,
		std::string fname_,
		CGAL::Point_3<K> p_min,
		CGAL::Point_3<K> p_max,
		const int &digits){

	typedef typename CGAL::Point_3<K> 				Point_3;
	typedef typename CGAL::Vector_3<K> 				Vector_3;
	typedef typename CGAL::Surface_mesh<Point_3> 	Surface_mesh;
	typedef	typename std::list<Ell_>::iterator 		Ell_Iterator;
	typedef typename Surface_mesh::Vertex_index 	vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 	vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 	vertex_range;
	typedef typename std::list<Point_3> 			list_points;
	typedef typename std::list<Point_3>::iterator 	Pts_Iterator;

	for (int it = 0; it !=ell_list.size(); it++){

		Ell_			cEll_;			// Current ellipse

		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		// Remove near points
		//NT tolerance(1E-16);
		NT tolerance(1E-7); // Important parameter (Be careful _ NTD)
		cEll_.E_mesh = CGAL::remove_near_points_sm(cEll_.E_mesh, tolerance);
		cEll_.E_mesh = CGAL::remove_near_points_sm(cEll_.E_mesh, tolerance);

		NT x1,y1,z1,x2,y2,z2;
		Vector_3 vect_, transvect, e3(CGAL::UnitVector_3<K>(3));
		vect_ = cEll_.E_normalvect;
		vect_ = CGAL::Normalize<K, NT>(vect_);

		/*  Find Angle of Rotation*/
		NT theta = CGAL::MyAngle<K, NT>(vect_,e3);
		NT c = cos(CGAL::to_double(theta));
		NT s = sin(CGAL::to_double(theta));

		/* Find angle to Rotate into xy plane */
		Vector_3 rotvect = CGAL::MyCross<K>(e3,vect_);
		Vector_3 u = CGAL::Normalize<K,NT>(rotvect);
		NT dk(CGAL::max(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-16){
			dk = CGAL::min(p_min[0], p_max[0]);
		}

		x1 = - dk * u[0]; x2 = dk * u[0];
		y1 = - dk * u[1]; y2 = dk * u[1];
		z1 = - dk * u[2]; z2 = dk * u[2];

		//transvect = Vector_3(x2-x1, y2-y1, z2-z1);
		//std::cout << " z2-z1 = " << z2-z1 << std::endl;

		std::string fname;
		fname = fname_ + CGAL::int2str_setw(it + 1,digits) + ".inp";
		std::ofstream cE_name_ss(fname.c_str()) ;

		// Convert all polygon into counterclockwise-oriented polygon
		//cEll_.E_mesh = CGAL::reverse_polygon_point_order< K, NT >(cEll_.E_mesh);
		//bool testTheta( (0 <= theta && theta <= 90) ||  (-180 <= theta && theta <= -90) );

		/* Convert angle to degrees*/
		//NT theta_ = CGAL::Rad2Deg(theta);
		//bool testTheta( 120 < theta_ && theta_ < 180) ;
		//cEll_.E_mesh = CGAL::reverse_polygon_point_order_old_2< K, NT >(cEll_.E_mesh, testTheta);

		Surface_mesh E_mesh_xy;
		list_points points;

		unsigned int i = 0, end = cEll_.E_mesh.number_of_vertices();
		for( ; i < end; ++i) {
			vertex_descriptor vh(i);
			points.push_back(CGAL::PointRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect));
		}

		for( Pts_Iterator iter=points.begin(); iter !=points.end(); ++iter) {
			E_mesh_xy.add_vertex(*iter);
		}

		vertex_range m_range = E_mesh_xy.vertices();
		E_mesh_xy.add_face(m_range);

		if ( !CGAL::isCounterClockwise<K,NT>(E_mesh_xy) ){
			E_mesh_xy = CGAL::reverse_polygon_point_order< K, NT >(E_mesh_xy);
		}

		// Export files
		cE_name_ss << cEll_.E_mesh.num_vertices() << " 0  0 0 0 \n";

		// There is a precision conflict between C++ and Fortran library (The points read by Fortran program do not lie in a plane)
		unsigned int iter = 0, iter_end = cEll_.E_mesh.number_of_vertices();
		for( ; iter < iter_end; ++iter) {
			vertex_descriptor vh(iter);
			vertex_descriptor v0(0);
			/*
				cE_name_ss << iter + 1 << "  " << std::setprecision(6) << std::uppercase << std::scientific <<
									CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect) << std::endl;
			 */
			/*
				cE_name_ss << iter + 1 << "  "  << std::setprecision(12) << std::uppercase << std::scientific <<
						CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[0]  << "  " <<
						CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[1]  << "  " <<
						CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[2]  << std::endl;
			 */

			/*
			cE_name_ss << iter + 1 << "  "  << std::setprecision(12) << std::uppercase << std::scientific <<
								CGAL::PointRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[0]  << "  " <<
								CGAL::PointRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[1]  << "  " <<
								CGAL::PointRotation<K,NT>(cEll_.E_mesh.point(vh), vect_,transvect)[2]  << std::endl;
			*/

			// NTD Modification 17/01/2017
			//cE_name_ss << iter + 1 << "  "  << std::setprecision(12) << std::uppercase << std::scientific <<
			cE_name_ss << iter + 1 << "  "  << std::setprecision(12) << std::uppercase << std::scientific <<
					E_mesh_xy.point(vh)[0]  << "  " <<
					E_mesh_xy.point(vh)[1]  << "  " <<
					E_mesh_xy.point(vh)[2]  << std::endl;

			// Force the Z-coordinate of the export point <~~~~~~~~~~~~~~~~~~~~~ [IMP]
			//CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(v0), vect_,transvect)[2]  << std::endl;

		}

		/*
		//int pindex(1);
		for (vertex_iterator iter = cEll_.E_mesh.vertices_begin(); iter != cEll_.E_mesh.vertices_end(); iter++){
			cE_name_ss << pindex << "  " << std::uppercase << std::scientific <<
					CGAL::PointInverseRotation<K,NT>(cEll_.E_mesh.point(*iter), vect_,transvect) <<std::endl;
			pindex++;
		}
		 */
		cE_name_ss.close();
	}
}

template <typename Ell_, typename K, typename NT>
void write_normal_vertors(
		std::list<Ell_> & ell_list,
		const std::string &filename_){

	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	typedef typename CGAL::Vector_3<K> Vector_3;
	typedef typename CGAL::Point_3<K> Point_3;

	const char* filename = filename_.c_str();
	std::ofstream    outfile(filename);
	Vector_3 vect_;

	// header
	outfile << "# Number of polygons" << std::endl;
	outfile << ell_list.size() << std::endl;
	outfile << "\n";

	outfile << "# Normal vectors ## npoly ## V[0] ## V[1] ## V[2]##" << std::endl;

	int inum(1);
	for (Ell_Iterator it = ell_list.begin(); it != ell_list.end(); it++){
		vect_ = (*it).E_normalvect;
		vect_ = CGAL::Normalize<K, NT>(vect_);

		outfile << inum << " " << std::setprecision(12) << vect_[0] << " " <<  vect_[1] << " " << vect_[2] << std::endl;

		inum++;
	}
	outfile.close();
}

}

#endif /* WRITE_PARAMS_FILE_FOR_LAGRIT_HH_ */
