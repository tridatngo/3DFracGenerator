/*
 * check_poly_poly_intersection.hh
 *
 *  Created on: Sep 29, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CHECK_POLY_POLY_INTERSECTION_HH_
#define INCLUDE_CHECK_POLY_POLY_INTERSECTION_HH_

#include "list_seg_index_unique.hh"

namespace CGAL{
template <typename K, typename Output_PlPl>
Output_PlPl check_poly_poly_intersection(CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon1,
										CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon2,
										const CGAL::Line_3<K> &Line_Int,
										bool &VERY_VERBOSE, bool & VERY_VERY_VERBOSE) {

	typedef typename K::Point_3                    		Point_3;
	typedef typename K::Vector_3                   		Vector_3;
	typedef typename K::Line_3                     		Line_3;
	typedef typename K::Segment_3                  		Segment_3;
	typedef typename K::Plane_3                  		Plane_3;
	typedef typename K::Intersect_3 					Intersect_3;
	typedef typename CGAL::Surface_mesh<Point_3>        Surface_mesh_;
	typedef typename Surface_mesh_::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh_::Vertex_iterator		vertex_iterator;

	typedef typename std::list<Line_3> 					LineList;
	typedef typename std::list<Line_3>::iterator 		LineList_Iter;
	typedef typename std::list<Point_3> 				PointList;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;

	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Line_3) >::type	result_type_ll;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type	result_type_lp;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Plane_3, Plane_3) >::type	result_type_pp;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Segment_3) >::type result_type_ls;

	Output_PlPl output_;

	int num_intersec_1(0), num_intersec_2(0);
	PointList intersection_points_1, intersection_points_2;
	std::list<int> list_seg_interst1, list_seg_interst2;

	bool no_polygons_interst(false), point_interst(false), point_interst2(false);
	std::list<bool> no_ls_interst_vec1, no_ls_interst_vec2;


	if(VERY_VERBOSE == 1){
		std::cout << "###########################" << std::endl;
	}

	{
		unsigned int i = 0, end = polygon1.number_of_vertices();
		for( ; i < end; ++i) {
			vertex_descriptor vb_seg1(i);
			vertex_descriptor ve_seg1(i+1);
			vertex_descriptor ve_seg1_0(0);

			if (i==end-1)
				ve_seg1 = ve_seg1_0;

			if(VERY_VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "(3) Re-check intersection between two planes." << std::endl;
				std::cout << std::setprecision(30) << Line_Int << std::endl;
			}

			Segment_3 mesh_1_seg = Segment_3(polygon1.point(vb_seg1), polygon1.point(ve_seg1));
			result_type_ls result_sl = CGAL::intersection(mesh_1_seg, Line_Int);

			if(VERY_VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "(4) Re-check intersection between two planes." << std::endl;
				std::cout << std::setprecision(30) << Line_Int << std::endl;
			}

			if (result_sl) {
				if (const Segment_3* s = boost::get<Segment_3>(&*result_sl)) {
					if(VERY_VERBOSE == 1){
						std::cout << "The segment lies on the plane-plane intersection." << std::endl;
						//std::cout << *s << std::endl;
					}
				}
				else
					if (const Point_3* p = boost::get<Point_3 >(&*result_sl)) {
						//std::cout << "The segment and the plane-plane intersection intersect on a point." << std::endl;
						//std::cout << *p << std::endl;
						point_interst =true;
						if(VERY_VERBOSE == 1){
							std::cout << "The segment and the plane-plane intersection intersect on a point: "<<
									*p << " on the segment " << i << std::endl;
						}
						list_seg_interst1.push_back(i);
						intersection_points_1.push_back(*p);
					}
				no_ls_interst_vec1.push_back(false);
			}

			else {
				no_ls_interst_vec1.push_back(true); // No intersection between line-segment
				//std::cout << "The segment and the plane-plane intersection do not intersect." << std::endl;
			}

		}

		if(VERY_VERBOSE == 1){
			std::cout << "###########################" << std::endl;
		}
		if (point_interst) {

			// Remove duplicate points

			// NTD 13/09/2016
			list_seg_index_unique<Point_3, NT> (intersection_points_1, list_seg_interst1);

			if(VERY_VERBOSE == 1){
				std::cout << "Size of intersection_points_1 : " << intersection_points_1.size() <<std::endl;
			}
			int it_(0);

			for (Pts_Iterator it=intersection_points_1.begin(); it!=intersection_points_1.end(); ++it){

				std::list<int>::iterator it_seg = list_seg_interst1.begin();
				std::advance(it_seg, it_);

				if(VERY_VERBOSE == 1){
					//std::cout << ' ' << *it << std::endl;
					std::cout << "The segment and the plane-plane intersection intersect on a point: "<<
							*it << " on the segment " << *it_seg << std::endl;
				}
				num_intersec_1 ++;
				it_++;
			}

			if(VERY_VERBOSE == 1){
				std::cout << "Number of intersection between the polygon 1 and the intersection line = " << num_intersec_1 << std::endl;
			}
		}
		else{
			if(VERY_VERBOSE == 1){
				std::cout << "The polygon 1 and the intersection line do not intersect." << std::endl;
			}
			num_intersec_1 = 0;
		}
	}

	if(VERY_VERY_VERBOSE == 1){
		//std::cout << *s << std::endl;}
		std::cout << "(5) Re-check intersection between two planes." << std::endl;
		std::cout << std::setprecision(30) << Line_Int << std::endl;
	}

	//##########################
	// Ellipse 2

	if(VERY_VERBOSE == 1){
		std::cout << "###########################" << std::endl;
	}
	{
		unsigned int i = 0, end = polygon2.number_of_vertices();
		for( ; i < end; ++i) {
			vertex_descriptor vb_seg2(i);
			vertex_descriptor ve_seg2(i+1);
			vertex_descriptor ve_seg2_0(0);

			if (i==end-1)
				ve_seg2 = ve_seg2_0;

			//std::cout << "Check intersection on the line " << Line_Int << std::endl;

			Segment_3 mesh_2_seg = Segment_3(polygon2.point(vb_seg2), polygon2.point(ve_seg2));

			if(VERY_VERBOSE == 1){
				std::cout << "Check intersection between the line " << Line_Int << " and the segment "<< mesh_2_seg << std::endl;
			}
			result_type_ls result_sl2 = CGAL::intersection(mesh_2_seg, Line_Int);

			if (result_sl2) {
				if (const Segment_3* s = boost::get<Segment_3>(&*result_sl2)) {
					if(VERY_VERBOSE == 1){
						std::cout << "The segment lies on the plane-plane intersection." << std::endl;
						//std::cout << *s << std::endl;
					}
				}
				else
					if (const Point_3* p = boost::get<Point_3 >(&*result_sl2)) {
						//std::cout << "The segment and the plane-plane intersection intersect on a point." << std::endl;
						//std::cout << *p << std::endl;
						point_interst2 =true;
						if(VERY_VERBOSE == 1){
							std::cout << "The segment and the plane-plane intersection intersect on a point: "<<
									*p << " on the segment " << i << std::endl;
						}
						list_seg_interst2.push_back(i);
						intersection_points_2.push_back(*p);
					}
				no_ls_interst_vec2.push_back(false);
			}

			else {
				no_ls_interst_vec2.push_back(true); // No intersection between line-segment
				//std::cout << "The segment and the plane-plane intersection do not intersect." << std::endl;
			}
		}
		if(VERY_VERBOSE == 1){
			std::cout << "###########################" << std::endl;
		}
		if (point_interst2) {

			// Remove duplicate points

			// NTD 13/09/2016
			list_seg_index_unique<Point_3, NT>(intersection_points_2, list_seg_interst2);

			if(VERY_VERBOSE == 1){

				std::cout << "Size of intersection_points_2 : " << intersection_points_2.size() <<std::endl;
			}
			int it_(0);

			for (Pts_Iterator it=intersection_points_2.begin(); it!=intersection_points_2.end(); ++it){

				std::list<int>::iterator it_seg = list_seg_interst2.begin();
				std::advance(it_seg, it_);

				if(VERY_VERBOSE == 1){
					//std::cout << ' ' << *it << std::endl;
					std::cout << "The segment and the plane-plane intersection intersect on a point: "<<
							*it << " on the segment " << *it_seg << std::endl;
				}

				num_intersec_2 ++;
				it_++;
			}
			if(VERY_VERBOSE == 1){
				std::cout << "Number of intersection between the polygon 2 and the intersection line = " << num_intersec_2 << std::endl;

			}
		}
	}
}
} // end of namespace

#endif /* INCLUDE_CHECK_POLY_POLY_INTERSECTION_HH_ */
