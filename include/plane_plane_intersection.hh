/*
 * plane_plane_intersection.hh
 *
 *  Created on: 9 août 2016
 *      Author: ngotr
 */

#ifndef PLANE_PLANE_INTERSECTION_HH_
#define PLANE_PLANE_INTERSECTION_HH_

#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Kernel/global_functions.h>
#include <list>
#include "list_seg_index_unique.hh"
#include "modify_point_3.hh"

#include "include/Exact_predicates_exact_constructions_kernel_ntd.h"
#include "include/Exact_predicates_inexact_constructions_kernel_ntd_MP_Float.h"

//#include "include/Exact_predicates_inexact_constructions_kernel_ntd_float.h"

#define VERY_VERBOSE 0
#define VERY_VERY_VERBOSE 0

#ifndef USE_CONVERT_TO_EPECK
#define USE_CONVERT_TO_EPECK 1
#endif

namespace CGAL {

template <typename K, typename Output_IP>
Output_IP intersections_of_two_polys (CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon1, CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon2){

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

	Output_IP output_;

	Plane_3 pl_polygon1,pl_polygon2; // Planes which pass through the polygon1 and polygon2
	Point_3 pt_1_polygon1, pt_2_polygon1, pt_3_polygon1; 		// Choose 3 point to define the plane of polygon1
	Point_3 pt_1_polygon2, pt_2_polygon2, pt_3_polygon2; 		// Choose 3 point to define the plane of polygon2

	int num_intersec_1(0), num_intersec_2(0);
	std::list<Point_3> intersection_points_1, intersection_points_2;
	std::list<int> list_seg_interst1, list_seg_interst2;

	bool no_polygons_interst(false), point_interst(false), point_interst2(false);
	std::list<bool> no_ls_interst_vec1, no_ls_interst_vec2;

	{
		if (polygon1.number_of_vertices() ==3){
			unsigned int i = 0, end = polygon1.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				if (i==0){
					pt_1_polygon1 = polygon1.point(vh);
				}
				else if(i==1){
					pt_2_polygon1 = polygon1.point(vh);
				}
				else{
					pt_3_polygon1 = polygon1.point(vh);
				}
			}
		}
		else{
			unsigned int i = 0, end = polygon1.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				//if (i==polygon1.num_vertices/2)
				if (i==0)
					pt_1_polygon1 = polygon1.point(vh);
				else
					if(i==round(polygon1.num_vertices()/4))
						pt_2_polygon1 = polygon1.point(vh);
					else
						if (i==round(polygon1.num_vertices()/2))
							pt_3_polygon1 = polygon1.point(vh);
				//std::cout << polygon1.point(vh) << std::endl;
			}
		}
	}

	{
		if (polygon2.number_of_vertices() ==3){
			unsigned int i = 0, end = polygon2.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);

				if (i==0){
					pt_1_polygon2 = polygon2.point(vh);
				}
				else if(i==1){
					pt_2_polygon2 = polygon2.point(vh);
				}
				else{
					pt_3_polygon2 = polygon2.point(vh);
				}
			}
		}
		else{
			unsigned int i = 0, end = polygon2.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				if (i==0)
					pt_1_polygon2 = polygon2.point(vh);
				else if(i==round(polygon2.num_vertices()/4))
					pt_2_polygon2 = polygon2.point(vh);
				else if (i==(polygon2.num_vertices()/2))
					pt_3_polygon2 = polygon2.point(vh);
			}
		}
	}

	pl_polygon1 = Plane_3(pt_1_polygon1, pt_2_polygon1, pt_3_polygon1);
	pl_polygon2 = Plane_3(pt_1_polygon2, pt_2_polygon2, pt_3_polygon2);

	if(VERY_VERBOSE == 1){
		std::cout << "Plane 1 : " << pl_polygon1 << std::endl;
		std::cout << "Plane 2 : " << pl_polygon2 << std::endl;
	}
	// ###########################################################################
	// Check if two planes intersect
	// ###########################################################################

	result_type_pp result = CGAL::intersection(pl_polygon1, pl_polygon2);

	Line_3 Line_Int;

	if (result) {
		//if (const Line_3* li = boost::get<Line_3>(&*result)) {
		if (const Line_3* li = boost::get<Line_3>(&*result)) {
			if(VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "Two planes intersect in a line" << std::endl;
				std::cout << std::setprecision(20) << *li << std::endl;
			}
			Line_Int = *li;}
		else {
			const Plane_3* pl = boost::get<Plane_3 >(&*result);
			if(VERY_VERBOSE == 1){
				//std::cout << *p << std::endl;
				std::cout << "Two planes are coincident" << std::endl;}
		}

		// ###########################################################################
		// Check if plane-plane intersection intersect the polygons
		// ###########################################################################

		// Ellipse 1

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

				Segment_3 mesh_1_seg = Segment_3(polygon1.point(vb_seg1), polygon1.point(ve_seg1));
				result_type_ls result_sl = CGAL::intersection(mesh_1_seg, Line_Int);

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
								std::cout << "The segment "<< i <<" and the plane-plane intersection intersect on a point: "<<
										*p << std::endl;
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
						std::cout << "The segment " << *it_seg << " and the plane-plane intersection intersect on a point: "<<
								*it <<  std::endl;
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

				result_type_ls  result_sl2 = CGAL::intersection(mesh_2_seg, Line_Int);

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
								std::cout << "The segment "<< i <<" and the plane-plane intersection intersect on a point: "<<
										*p << std::endl;
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
						std::cout << "The segment " << *it_seg<< " and the plane-plane intersection intersect on a point: "<<
								*it << std::endl;
					}

					num_intersec_2 ++;
					it_++;
				}
				if(VERY_VERBOSE == 1){
					std::cout << "Number of intersection between the polygon 2 and the intersection line = " << num_intersec_2 << std::endl;

				}
			}
		}

		output_.num_intersec_1_ = num_intersec_1;
		output_.num_intersec_2_ = num_intersec_2;
		output_.intersection_points_1_ = intersection_points_1;
		output_.intersection_points_2_ = intersection_points_2;
		output_.list_seg_interst1_ = list_seg_interst1;
		output_.list_seg_interst2_ = list_seg_interst2;
	}
	else{
		if(VERY_VERBOSE == 1){
			std::cout << "Two planes are parallel" << std::endl;}
	}

	return output_;
}

template <typename K, typename Output_IPL>
Output_IPL intersections_of_two_polys_exp_line(CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon1, CGAL::Surface_mesh< CGAL::Point_3<K> > &polygon2) {

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

	typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
	typedef EPECK::Point_3                    			Point_3_epeck;
	typedef EPECK::Line_3                     			Line_3_epeck;
	typedef EPECK::Segment_3                  			Segment_3_epeck;
	typedef EPECK::Plane_3                  			Plane_3_epeck;
	typedef EPECK::Intersect_3 							Intersect_3_epeck;
	typedef CGAL::Surface_mesh<Point_3_epeck>       	Surface_mesh_epeck;
	typedef typename Surface_mesh_epeck::Vertex_index 	vertex_descriptor_epeck;

	typedef CGAL::Cartesian_converter<K,EPECK > Converter;

	struct Convert_sm{
		//mutable bool first_vertex;
		//Convert_sm():first_vertex(true) {}
		Convert_sm() {}

		Surface_mesh_epeck operator()(const Surface_mesh_&) const { return Surface_mesh_epeck(); }

		void operator()(const Surface_mesh_& src, Surface_mesh_epeck& tgt) const
		{
			unsigned int i = 0, end = src.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				src.point(vh);
				tgt.add_vertex( Converter()( src.point(vh) ));
			}

			Surface_mesh_epeck::Vertex_range current_m_range = tgt.vertices();
			tgt.add_face(current_m_range);
		}
	};

	Output_IPL output_;
	bool bool_intersecting(true);

	Plane_3 pl_polygon1,pl_polygon2; // Planes which pass through the polygon1 and polygon2
	Point_3 pt_1_polygon1, pt_2_polygon1, pt_3_polygon1; 		// Choose 3 point to define the plane of polygon1
	Point_3 pt_1_polygon2, pt_2_polygon2, pt_3_polygon2; 	// Choose 3 point to define the plane of polygon2

	int num_intersec_1(0), num_intersec_2(0);
	std::list<Point_3> intersection_points_1, intersection_points_2;
	std::list<int> list_seg_interst1, list_seg_interst2;

	bool no_polygons_interst(false), point_interst(false), point_interst2(false);
	std::list<bool> no_ls_interst_vec1, no_ls_interst_vec2;

	{
		if (polygon1.number_of_vertices() ==3){
			unsigned int i = 0, end = polygon1.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				if (i==0){
					pt_1_polygon1 = polygon1.point(vh);
				}
				else if(i==1){
					pt_2_polygon1 = polygon1.point(vh);
				}
				else{
					pt_3_polygon1 = polygon1.point(vh);
				}
			}
		}
		else{
			unsigned int i = 0, end = polygon1.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				//if (i==polygon1.num_vertices/2)
				if (i==0)
					pt_1_polygon1 = polygon1.point(vh);
				else
					if(i==round(polygon1.num_vertices()/4))
						pt_2_polygon1 = polygon1.point(vh);
					else
						if (i==round(polygon1.num_vertices()/2))
							pt_3_polygon1 = polygon1.point(vh);
				//std::cout << polygon1.point(vh) << std::endl;
			}
		}
	}

	{
		if (polygon2.number_of_vertices() ==3){
			unsigned int i = 0, end = polygon2.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);

				if (i==0){
					pt_1_polygon2 = polygon2.point(vh);
				}
				else if(i==1){
					pt_2_polygon2 = polygon2.point(vh);
				}
				else{
					pt_3_polygon2 = polygon2.point(vh);
				}
			}
		}
		else{
			unsigned int i = 0, end = polygon2.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor vh(i);
				if (i==0)
					pt_1_polygon2 = polygon2.point(vh);
				else if(i==round(polygon2.num_vertices()/4))
					pt_2_polygon2 = polygon2.point(vh);
				else if (i==(polygon2.num_vertices()/2))
					pt_3_polygon2 = polygon2.point(vh);
			}
		}
	}

	pl_polygon1 = Plane_3(pt_1_polygon1, pt_2_polygon1, pt_3_polygon1);
	pl_polygon2 = Plane_3(pt_1_polygon2, pt_2_polygon2, pt_3_polygon2);

	if(VERY_VERBOSE == 1){
		std::cout << "Plane 1 : " << pl_polygon1 << std::endl;
		std::cout << "Plane 2 : " << pl_polygon2 << std::endl;
	}
	// ###########################################################################
	// Check if two planes intersect
	// ###########################################################################

	result_type_pp result = CGAL::intersection(pl_polygon1, pl_polygon2);

	Line_3 Line_Int, Line_Int_Mod;
	Line_3 test_Line_Int;

	if (result) {
		//if (const Line_3* li = boost::get<Line_3>(&*result)) {
		if (const Line_3* li = boost::get<Line_3>(&*result)) {
			if(VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "Two planes intersect in a line" << std::endl;
				std::cout << std::setprecision(20) << *li << std::endl;
			}
			Line_Int = *li;

			if(VERY_VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "(1) Re-check intersection between two planes." << std::endl;
				std::cout << std::setprecision(20) << Line_Int << std::endl;
			}
		}
		else {
			const Plane_3* pl = boost::get<Plane_3 >(&*result);
			if(VERY_VERBOSE == 1){
				//std::cout << *p << std::endl;
				std::cout << "Two planes are coincident." << std::endl;
			}
			bool_intersecting = false;
		}
	}
	else{

		if(VERY_VERBOSE == 1){
			std::cout << "Two planes are parallel." << std::endl;
		}
		bool_intersecting = false;
	}

	// After calling CGAL::intersection(), Line_Int will be changed. This command return the original Lint_Int
	test_Line_Int = Line_Int;
	const Line_3 cst_Line_Int(Line_Int);

	if(VERY_VERY_VERBOSE == 1){
		//std::cout << *s << std::endl;}
		std::cout << "(2) Re-check intersection between two planes." << std::endl;
		std::cout << std::setprecision(20) << Line_Int << std::endl;
	}


	// ==============================================================================================
	// EPECK kernel

	Line_3_epeck Line_Int_Epeck;

	if (USE_CONVERT_TO_EPECK == 1 && bool_intersecting){
		Surface_mesh_epeck poly_1_epeck, poly_2_epeck;

		Convert_sm()(polygon1, poly_1_epeck);
		Convert_sm()(polygon2, poly_2_epeck);

		//std::cout << "Test EPECK." << std::endl;
		//std::cout << "poly_1_epeck.num_vetices() = " << poly_1_epeck.number_of_vertices()  << std::endl;

		Point_3_epeck pt_1_poly_1_epeck, pt_2_poly_1_epeck, pt_3_poly_1_epeck;
		Point_3_epeck pt_1_poly_2_epeck, pt_2_poly_2_epeck, pt_3_poly_2_epeck;

		{
			unsigned int i = 0, end = poly_1_epeck.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor_epeck vh(i);
				if (i==0)
					pt_1_poly_1_epeck = poly_1_epeck.point(vh);
				else if(i==round(poly_1_epeck.num_vertices()/4))
					pt_2_poly_1_epeck = poly_1_epeck.point(vh);
				else if (i==(poly_1_epeck.num_vertices()/2))
					pt_3_poly_1_epeck = poly_1_epeck.point(vh);
			}
		}

		{
			unsigned int i = 0, end = poly_2_epeck.number_of_vertices();
			for( ; i < end; ++i) {
				vertex_descriptor_epeck vh(i);
				//std::cout << "EPECK: poly_2_epeck.point(vh) = " << poly_2_epeck.point(vh)<< std::endl;
				if (i==0)
					pt_1_poly_2_epeck = poly_2_epeck.point(vh);
				else if(i==round(poly_2_epeck.num_vertices()/4))
					pt_2_poly_2_epeck = poly_2_epeck.point(vh);
				else if (i==(poly_2_epeck.num_vertices()/2))
					pt_3_poly_2_epeck = poly_2_epeck.point(vh);
			}
		}

		Plane_3_epeck test_plane1_epeck = CGAL::Plane_3<EPECK>(pt_1_poly_1_epeck, pt_2_poly_1_epeck, pt_3_poly_1_epeck);
		Plane_3_epeck test_plane2_epeck = CGAL::Plane_3<EPECK>(pt_1_poly_2_epeck, pt_2_poly_2_epeck, pt_3_poly_2_epeck);

		CGAL::cpp11::result_of<EPECK::Intersect_3(EPECK::Plane_3, EPECK::Plane_3)>::type
				result_epec_pp_1 = CGAL::intersection(test_plane1_epeck, test_plane2_epeck);

		if (const Line_3_epeck* l = boost::get<Line_3_epeck >(&*result_epec_pp_1)) {
			Line_Int_Epeck = *l;
			if (VERY_VERBOSE == 1){
				std::cout << "(1) EPECK: Found intersection between planes " << std::endl;
				std::cout << std::setprecision(30) << *l << std::endl;
				std::cout << "\n" << std::endl;
			}
		}
	}

	// ==============================================================================================

	// ###########################################################################
	// Check if plane-plane intersection intersect the polygons
	// ###########################################################################

	// Ellipse 1

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
				std::cout << std::setprecision(20) << Line_Int << std::endl;
			}

			Segment_3 mesh_1_seg = Segment_3(polygon1.point(vb_seg1), polygon1.point(ve_seg1));
			result_type_ls result_sl = CGAL::intersection(mesh_1_seg, test_Line_Int);

			if(VERY_VERY_VERBOSE == 1){
				//std::cout << *s << std::endl;}
				std::cout << "(4) Re-check intersection between two planes." << std::endl;
				std::cout << std::setprecision(20) << Line_Int << std::endl;

				std::cout << "(5) Re-check intersection between two planes." << std::endl;
				std::cout << std::setprecision(20) << cst_Line_Int << std::endl;

				//std::cout << std::setprecision(20) << *li << std::endl;
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
							std::cout << "The segment " << i <<" and the plane-plane intersection intersect on a point: "<<
									*p << std::endl;
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
		std::cout << std::setprecision(20) << Line_Int << std::endl;
	}

	//##########################
	// Ellipse 2

	test_Line_Int = Line_Int;

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
							std::cout << "The segment " << i <<" and the plane-plane intersection intersect on a point: "<<
									*p << std::endl;
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
					std::cout << "The segment " << *it_seg << " and the plane-plane intersection intersect on a point: "<<
							*it << std::endl;
				}

				num_intersec_2 ++;
				it_++;
			}
			if(VERY_VERBOSE == 1){
				std::cout << "Number of intersection between the polygon 2 and the intersection line = " << num_intersec_2 << std::endl;

			}
		}
	}

	// NTD 21/09/2016: Imprecision of intersection between lines consisting of close-to-zero coordinate point
	//NT eps(1E-30);
	//Point_3 pmod_1, pmod_2;
	//pmod_1 = Line_Int.point(1);
	//pmod_1 = CGAL::zero_p3<Point_3, NT>(pmod_1, eps);
	//pmod_2 = Line_Int.point(2);

	if(VERY_VERY_VERBOSE == 1){
		//std::cout << *s << std::endl;}
		std::cout << "(10) Re-check intersection between two planes." << std::endl;
		std::cout << std::setprecision(20) << Line_Int << std::endl;
	}


	// ==============================================================================================
	// EPECK kernel
	/*
	if (VERY_VERBOSE == 1){
		CGAL::cpp11::result_of<EPECK::Intersect_3(EPECK::Plane_3, EPECK::Plane_3)>::type
				result_epec_pp_1 = CGAL::intersection(test_plane1_epeck, test_plane2_epeck);

		if (const Line_3_epeck* l = boost::get<Line_3_epeck >(&*result_epec_pp_1)) {
			std::cout << "(2) EPECK: Found intersection between planes " << std::endl;
			std::cout << std::setprecision(30) << *l << std::endl;
			std::cout << "\n" << std::endl;
		}
	}
	 */
	// ==============================================================================================

	output_.num_intersec_1_ 		= num_intersec_1;
	output_.num_intersec_2_ 		= num_intersec_2;
	output_.intersection_points_1_ 	= intersection_points_1;
	output_.intersection_points_2_ 	= intersection_points_2;
	output_.list_seg_interst1_ 		= list_seg_interst1;
	output_.list_seg_interst2_ 		= list_seg_interst2;
	output_.intersection_line_ 		= Line_Int;
	output_.intersection_line_epeck_= Line_Int_Epeck;
	output_.bool_intersecting_ 		= bool_intersecting;

	return output_;
}
} // end namespace CGAL

#endif /* PLANE_PLANE_INTERSECTION_HH_ */
