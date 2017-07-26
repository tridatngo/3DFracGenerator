/*
 * check_conformity_poly_list.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 */


/*
 * This function check_conformity_poly_list() to check conformity of the vertex map
 * between the father polygon and the points at the intersection with other ones
 * This ensures that the child polygons are convex
 */

#ifndef CHECK_CONFORMITY_POLY_LIST_HH_
#define CHECK_CONFORMITY_POLY_LIST_HH_

// Minimum & maximum of a list of numbers
#include "Min_Max_List.hh"

//Transpose a list
#include "transpose_a_list.hh"

// 2D algorithm on 3D data
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Polygon_2_algorithms.h>


template <typename SurfaceMesh_, typename Point_3>
void check_conformity_poly_list(SurfaceMesh_ &poly, std::list<Point_3> &list, std::list<int> &seg_of_endpoints){

	typedef typename SurfaceMesh_::Vertex_index 		   vertex_descriptor;

	unsigned int begin = Min_List<int>(seg_of_endpoints);
	unsigned int end = Max_List<int>(seg_of_endpoints);
	std::list<Point_3> list_new;

	vertex_descriptor vb_seg(end);
	vertex_descriptor ve_seg(end+1);
	vertex_descriptor ve_seg_0(0);

	if (end==poly.num_vertices()-1)
		ve_seg = ve_seg_0;

	if (CGAL::collinear(*list.begin(), poly.point(vb_seg), poly.point(ve_seg)) == 0){
		list_new = transpose_a_list< std::list<Point_3> >(list);;
	}
	list = list_new;
}



#endif /* CHECK_CONFORMITY_POLY_LIST_HH_ */
