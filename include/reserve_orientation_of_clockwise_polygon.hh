/*
 * check_orientation_of_polygon.hh
 *
 *  Created on: 10 sept. 2016
 *      Author: ngo
 */

#ifndef INCLUDE_RESERVE_ORIENTATION_OF_CLOCKWISE_POLYGON_HH_
#define INCLUDE_RESERVE_ORIENTATION_OF_CLOCKWISE_POLYGON_HH_

#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <iostream>
#include <list>

//Transpose a list
#include "transpose_a_list.hh"


namespace CGAL{

template < typename K, typename NT_>
bool isCounterClockwise (CGAL::Surface_mesh< CGAL::Point_3<K> > &poly_) {

	typedef typename CGAL::Point_3<K> Point_3;
	typedef typename CGAL::Surface_mesh< Point_3> Surface_mesh;
	typedef typename Surface_mesh::Vertex_index 		   vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator vertex_iterator;

	typedef typename CGAL::Point_2<K> Point_2;
	typedef typename CGAL::Polygon_2<K> Polygon_2;

	Polygon_2 poly_2_xy, poly_2_xz;

	unsigned int i = 0, end = poly_.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vh(i);
		poly_2_xy.push_back( Point_2( (poly_.point(vh)).x(), (poly_.point(vh)).y()) );
		poly_2_xz.push_back( Point_2( (poly_.point(vh)).x(), (poly_.point(vh)).z()) );
	}

	//std::cout << "Polygon orientation : " << poly_2.orientation() << std::endl;
	// COUNTERCLOCKWISE == NEGATIVE; CLOCKWISE == POSITIVE

	//bool IsCounterClockwise(poly_2.orientation() == CGAL::COUNTERCLOCKWISE);
	bool IsCounterClockwise(false);
	if (!poly_2_xy.is_collinear_oriented ()){
		IsCounterClockwise = bool(poly_2_xy.orientation() == CGAL::COUNTERCLOCKWISE);
	}
	else{
		IsCounterClockwise = bool(poly_2_xz.orientation() == CGAL::COUNTERCLOCKWISE);
	}

	return IsCounterClockwise;
}


template < typename K, typename NT_>
CGAL::Surface_mesh< CGAL::Point_3<K> > reverse_polygon_point_order_old_1 (CGAL::Surface_mesh< CGAL::Point_3<K> > &poly_) {

	typedef typename CGAL::Point_3<K> 					Point_3;
	typedef typename CGAL::Surface_mesh< Point_3> 		Surface_mesh;
	typedef typename std::list<Point_3> 				list_points;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename Surface_mesh::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 		vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 		vertex_range;

	Surface_mesh poly_new;
	list_points points;

	unsigned int i = 0, end = poly_.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vh(i);
		points.push_back(poly_.point(vh));
	}

	if ( isCounterClockwise<K,NT_>(poly_) ){
		poly_new = poly_;
	}
	else
	{

		points = transpose_a_list< std::list<Point_3> >(points);

		for( Pts_Iterator iter=points.begin(); iter !=points.end(); ++iter) {
			poly_new.add_vertex(*iter);
		}

		vertex_range m_range = poly_new.vertices();
		poly_new.add_face(m_range);

	}
	return poly_new;
}

template < typename K, typename NT_>
CGAL::Surface_mesh< CGAL::Point_3<K> > reverse_polygon_point_order_old_2 (CGAL::Surface_mesh< CGAL::Point_3<K> > &poly_, bool &bool_) {

	typedef typename CGAL::Point_3<K> 					Point_3;
	typedef typename CGAL::Surface_mesh< Point_3> 		Surface_mesh;
	typedef typename std::list<Point_3> 				list_points;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename Surface_mesh::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 		vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 		vertex_range;

	Surface_mesh poly_new;
	list_points points;

	unsigned int i = 0, end = poly_.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vh(i);
		points.push_back(poly_.point(vh));
	}
	if (!bool_){ // 120 < theta < 180 where theta is the angle between the normal vector of the poly and e3

		if ( isCounterClockwise<K,NT_>(poly_) ){
			poly_new = poly_;
		}
		else
		{

			points = transpose_a_list< std::list<Point_3> >(points);

			for( Pts_Iterator iter=points.begin(); iter !=points.end(); ++iter) {
				poly_new.add_vertex(*iter);
			}

			vertex_range m_range = poly_new.vertices();
			poly_new.add_face(m_range);
		}
	}
	else{
		if ( isCounterClockwise<K,NT_>(poly_) )
		{

			points = transpose_a_list< std::list<Point_3> >(points);

			for( Pts_Iterator iter=points.begin(); iter !=points.end(); ++iter) {
				poly_new.add_vertex(*iter);
			}

			vertex_range m_range = poly_new.vertices();
			poly_new.add_face(m_range);
		}
		else{
			poly_new = poly_;
		}
	}

	return poly_new;
}


template < typename K, typename NT_>
CGAL::Surface_mesh< CGAL::Point_3<K> > reverse_polygon_point_order (CGAL::Surface_mesh< CGAL::Point_3<K> > &poly_) {

	typedef typename CGAL::Point_3<K> 					Point_3;
	typedef typename CGAL::Surface_mesh< Point_3> 		Surface_mesh;
	typedef typename std::list<Point_3> 				list_points;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename Surface_mesh::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 		vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 		vertex_range;

	Surface_mesh poly_new;
	list_points points;

	unsigned int i = 0, end = poly_.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vh(i);
		points.push_back(poly_.point(vh));
	}

	points = transpose_a_list< std::list<Point_3> >(points);

	for( Pts_Iterator iter=points.begin(); iter !=points.end(); ++iter) {
		poly_new.add_vertex(*iter);
	}

	vertex_range m_range = poly_new.vertices();
	poly_new.add_face(m_range);

	return poly_new;
}

} // end of namespace

#endif /* INCLUDE_RESERVE_ORIENTATION_OF_CLOCKWISE_POLYGON_HH_ */
