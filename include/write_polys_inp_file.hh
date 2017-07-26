/*
 * write_polys_inp_file.hh
 *
 *  Created on: Oct 25, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_WRITE_POLYS_INP_FILE_HH_
#define INCLUDE_WRITE_POLYS_INP_FILE_HH_

#include <include/points_are_close.hh>
#include "myglobal_functions.hh"
#include "MyRotation.hh"
#include "write_params_file_for_lagrit.hh"
#include "Vector_3_Operations.hh"
#include <include/reserve_orientation_of_clockwise_polygon.hh>

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>

namespace CGAL{

template < typename K, typename Ell_, typename NT>
void write_polys_inp_file(
		const std::list<Ell_> &ell_list,
		const CGAL::Point_3<K> p_min,
		const CGAL::Point_3<K> p_max,
		const std::string &filename_) {

	typedef typename CGAL::Point_3<K> 				Point_3;
	typedef typename CGAL::Vector_3<K> 				Vector_3;
	typedef typename CGAL::Surface_mesh<Point_3> 	Surface_mesh;
	typedef	typename std::list<Ell_>::iterator 		Ell_Iterator;
	typedef typename Surface_mesh::Vertex_index 	vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 	vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 	vertex_range;
	typedef typename std::list<Point_3> 			list_points;
	typedef typename std::list<Point_3>::iterator 	Pts_Iterator;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis

	const char* filename = filename_.c_str();
	std::ofstream    fname_ofs(filename);

	std::list<Ell_> ell_list_ = ell_list;
	std::list<Ell_> ell_list_new;
	ell_list_new.clear();

	//=======================================================================================================================
	/*
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

		// Find Angle of Rotation
		NT theta = CGAL::MyAngle<K, NT>(vect_,e3);
		NT c = cos(CGAL::to_double(theta));
		NT s = sin(CGAL::to_double(theta));

		// Find angle to Rotate into xy plane
		Vector_3 rotvect = CGAL::MyCross<K>(e3,vect_);
		Vector_3 u = CGAL::Normalize<K,NT>(rotvect);
		NT dk(CGAL::max(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-16){
			dk = CGAL::min(p_min[0], p_max[0]);
		}

		x1 = - dk * u[0]; x2 = dk * u[0];
		y1 = - dk * u[1]; y2 = dk * u[1];
		z1 = - dk * u[2]; z2 = dk * u[2];

		Surface_mesh E_mesh_xy;
		Ell_ Ell_xy;
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

		Ell_xy.E_mesh = E_mesh_xy;

		ell_list_new.push_back(Ell_xy);
	}

	ell_list_.clear();
	ell_list_ = ell_list_new;
	ell_list_new.clear();
	*/

	//=======================================================================================================================
	int num_pts(0), curr_pts(0), curr_line(0);

	for (int i = 0; i !=ell_list_.size(); i++){
		Ell_			cEll_1;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list_.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;
		//std::cout << "Polygon name " << cEll_1.E_name << " -- number_of_vertices() = " << cEll_1.E_mesh.number_of_vertices()<< std::endl;
		num_pts = num_pts + cEll_1.E_mesh.number_of_vertices();
	}

	//std::cout << "num_pts " << num_pts << std::endl;

	fname_ofs << num_pts << " " << num_pts - ell_list_.size() << " 0 0 0\n";

	for (int i = 0; i !=ell_list_.size(); i++){
		Ell_			cEll_1;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list_.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;
		unsigned int it = 0, end = cEll_1.E_mesh.number_of_vertices();
		for( ; it < end; ++it) {
			vertex_descriptor vh(it);
			fname_ofs << curr_pts + it + 1 << " ";

			// NTD 12/01/2017
			//fname_ofs << cEll_1.E_mesh.point(vh)[0] << " " << cEll_1.E_mesh.point(vh)[1] << " "<< cEll_1.E_mesh.point(vh)[2];
			fname_ofs << std::setprecision(12) << std::uppercase << std::scientific << std::setw(20)
						<< cEll_1.E_mesh.point(vh)[0] << " " << cEll_1.E_mesh.point(vh)[1] << " "<< cEll_1.E_mesh.point(vh)[2];

			fname_ofs << "\n";
		}

		curr_pts = curr_pts + cEll_1.E_mesh.number_of_vertices();
	}

	for (int i = 0; i !=ell_list_.size(); i++){
		Ell_			cEll_1;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list_.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;

		unsigned int it = 0, end = cEll_1.E_mesh.number_of_vertices();

		for( ; it < end-1; ++it) {
			fname_ofs << curr_line + it + 1 << " ";
			fname_ofs << i + 1 << " line ";
			fname_ofs << curr_line + i + it + 1 << " " << curr_line + i + it + 2;
			fname_ofs << "\n";
		}

		curr_line = curr_line + cEll_1.E_mesh.number_of_vertices() - 1;
	}

	fname_ofs.close();
}
}

#endif /* INCLUDE_WRITE_POLYS_INP_FILE_HH_ */
