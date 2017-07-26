/*
 * my_line_3_intersection.hh
 *
 *  Created on: Sep 29, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_MY_LINE_3_INTERSECTION_HH_
#define INCLUDE_MY_LINE_3_INTERSECTION_HH_

#include "myglobal_functions.hh"
#include "Vector_3_Operations.hh"
namespace CGAL {

template <typename K, typename NT, typename Output_ll_Int>
Output_ll_Int my_line_3_intersection( CGAL::Line_3<K> &l1, CGAL::Line_3<K> &l2, NT eps_){

	typedef typename K::Point_3                    		Point_3;
	typedef typename K::Vector_3                   		Vector_3;
	typedef typename K::Line_3                     		Line_3;
	typedef typename K::Segment_3                  		Segment_3;
	typedef typename K::Plane_3                  		Plane_3;
	typedef typename K::Intersect_3 					Intersect_3;

	typedef typename std::list<Line_3> 					LineList;
	typedef typename std::list<Line_3>::iterator 		LineList_Iter;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Segment_3) >::type result_type;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Line_3) >::type result_type_ll;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type result_type_lp;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Plane_3, Plane_3) >::type result_type_pp;

	bool is_coincident(false), no_intersect(true);

	Output_ll_Int output;

	Point_3 PA, PB, Pout;


	Vector_3 v1( l1.to_vector ()), v2( l2.to_vector ());

	// get plane p1 that contains l1 and is parallel to l2
	Plane_3 p1( l1, l1.point(0) + l2.to_vector() );

	// get plane p2 that contains l1 and is perpendicular to p1
	Plane_3 p2( l1, l1.point(0) + p1.orthogonal_vector() );

	if ( l1 == l2){
		is_coincident = true;
	}
	else {

		// get intersection i1 between p2 and l2
		Point_3 i1;

		CGAL::Object result = CGAL::intersection( p2, l2 );
		if( !assign( i1, result ) ) {
			Line_3 il;
			if( assign( il, result ) ){
				std::cout << "Intersection between plane and line is a line --> l1 and l2 are parallel!" << std::endl;
			}
		}

		// get intersection i2 on l1
		Point_3 i2 = l1.projection( i1 );


		std::cout << "PA = " << i1 <<std::endl;
		std::cout << "PB = " << i2 <<std::endl;

		if ( CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(PA, PB))) < eps_){
			// We take the middle point: final point is (i1+i2)/2
			Pout = Point_3( (i1.x() + i2.x()) / 2.0, (i1.y() + i2.y()) / 2.0, (i1.z() + i2.z()) / 2.0 );
			no_intersect = false;
		}
	}

	/*
	Plane_3 plane_1_( Plane_3(l2.point(0), v1) );
	Plane_3 plane_2_( Plane_3(l1.point(0), v2) );


	if ( l1 == l2){
		is_coincident = true;
	}
	else {

		result_type_lp result_lpa = CGAL::intersection(l1, plane_1_);
		if (const Point_3* p = boost::get<Point_3 >(&*result_lpa)) {
			PA = *p;
		}

		result_type_lp result_lpb = CGAL::intersection(l2, plane_2_);
		if (const Point_3* p = boost::get<Point_3 >(&*result_lpb)) {
			PB = *p;
		}

		if ( CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(PA, PB))) < eps_){
			// We take the middle point
			Pout = Point_3( 0.5*(PA.x()+PB.x()), 0.5*(PA.y()+PB.y()), 0.5*(PA.z()+PB.z()) );
			no_intersect = false;
		}
	}
	 */

	output.is_coincident_ = is_coincident;
	output.no_intersect_ = no_intersect;
	output.pts_ = Pout;

	return output;
}
}
#endif /* INCLUDE_MY_LINE_3_INTERSECTION_HH_ */
