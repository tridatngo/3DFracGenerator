/*
 * check_conformity_newpoly.hh
 *
 *  Created on: 17 août 2016
 *      Author: ngotr
 */

/*
 * This function check_conformity_newpoly() to check conformity of the vertex map
 * between the father polygon and the points at the intersection with other ones
 * This ensures that the child polygons are simple (no self-intersecting)
 */

#ifndef CHECK_CONFORMITY_NEWPOLY_HH_
#define CHECK_CONFORMITY_NEWPOLY_HH_

//Transpose a list
#include "transpose_a_list.hh"
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <include/points_are_close.hh>

#define VERBOSE_NTD 0

namespace CGAL{

template <typename Kernel, typename SurfaceMesh_, typename Point_3>
void check_conformity_newpoly(SurfaceMesh_ &poly, std::list<Point_3> &list){

	typedef typename SurfaceMesh_::Vertex_iterator 		   	vertex_iterator;
	typedef typename SurfaceMesh_::Vertex_index 		   	vertex_descriptor;
	typedef typename std::list<Point_3>::iterator			List_iterator;

	std::list<Point_3> list_new;
	list_new =list;

	//unsigned int vd_poly_begin = poly.vertices_begin();
	unsigned int poly_end = poly.num_vertices() -1;

	vertex_descriptor ve_poly(poly_end);
	List_iterator list_end = list.begin();
	std::advance(list_end, list.size()-1);

	if (VERBOSE_NTD){
		std::cout << " " << std::endl;
		std::cout << "################################################" << std::endl;
		std::cout << "list.size() = " << list.size() << std::endl;
	}

	Point_3 checking_points[4] = {
			poly.point( *poly.vertices_begin() ),
			poly.point( ve_poly ),
			*list.begin(),
			*list_end };

	if (VERBOSE_NTD){
		std::cout << " " << std::endl;
		std::cout << "################################################" << std::endl;
		std::cout << "Checking points 1 : " << poly.point( *poly.vertices_begin()) << std::endl;
		std::cout << "Checking points 2 : " << poly.point( ve_poly ) << std::endl;
		std::cout << "Checking points 3 : " << *list.begin() << std::endl;
		std::cout << "Checking points 4 : " << *list_end << std::endl;
		std::cout << "################################################" << std::endl;
		std::cout << " " << std::endl;
	}
	/*
	if (points_are_close<Point_3, NT_MP>(poly.point( ve_poly ), *list_end)
			||  points_are_close<Point_3, NT_MP>(poly.point(*poly.vertices_begin()) , *list.begin())){
		std::cout << " Must transpose the list !" << std::endl;
		list_new = transpose_a_list< std::list<Point_3> >(list);
	}
	*/
	if (points_are_close<Point_3, NT_MP>(poly.point( ve_poly ),*list.begin())
			||  points_are_close<Point_3, NT_MP>(poly.point(*poly.vertices_begin()), *list_end)){
		//std::cout << "Must not transpose the list !" << std::endl;
		list_new = list;
	}
	else{

		bool bool_x =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_yz_3<Kernel>());

		bool bool_y =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_xz_3<Kernel>());

		bool bool_z =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_xy_3<Kernel>());

		if (VERBOSE_NTD){
			std::cout << " " << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << "Bool_x : " << bool_x << std::endl;
			std::cout << "Bool_y : " << bool_y << std::endl;
			std::cout << "Bool_z : " << bool_z << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;
		}

		if (!bool_x && !bool_y && !bool_z ){
			if (VERBOSE_NTD){
				std::cerr << "Error : Polygon is not simple." << std::endl;
			}
			list_new = transpose_a_list< std::list<Point_3> >(list);
		}
	}
	list = list_new;
}


template <typename Kernel, typename SurfaceMesh_, typename Point_3>
bool bool_conformity_newpoly(SurfaceMesh_ &poly, std::list<Point_3> &list){

	typedef typename SurfaceMesh_::Vertex_iterator 		   	vertex_iterator;
	typedef typename SurfaceMesh_::Vertex_index 		   	vertex_descriptor;
	typedef typename std::list<Point_3>::iterator			List_iterator;

	std::list<Point_3> list_new;
	list_new =list;
	bool conf_(true);

	//unsigned int vd_poly_begin = poly.vertices_begin();
	unsigned int poly_end = poly.num_vertices() -1;

	vertex_descriptor ve_poly(poly_end);
	List_iterator list_end = list.begin();
	std::advance(list_end, list.size()-1);

	Point_3 checking_points[4] = {
			poly.point( *poly.vertices_begin() ),
			poly.point( ve_poly ),
			*list.begin(),
			*list_end };

	if (points_are_close<Point_3, NT_MP>(poly.point( ve_poly ),*list.begin())
			||  points_are_close<Point_3, NT_MP>(poly.point(*poly.vertices_begin()), *list_end)){
		list_new = list;
	}
	else{

		bool bool_x =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_yz_3<Kernel>());

		bool bool_y =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_xz_3<Kernel>());

		bool bool_z =  CGAL::is_simple_2(checking_points,
				checking_points+4,
				CGAL::Projection_traits_xy_3<Kernel>());

		if (!bool_x && !bool_y && !bool_z ){
			conf_ = false;
		}
	}
	return conf_;
}

} // end off namspace

#endif /* CHECK_CONFORMITY_NEWPOLY_HH_ */
