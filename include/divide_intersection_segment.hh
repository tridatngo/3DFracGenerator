/*
 * divide_intersection_segment.hh
 *
 *  Created on: Sep 16, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_DIVIDE_INTERSECTION_SEGMENT_HH_
#define INCLUDE_DIVIDE_INTERSECTION_SEGMENT_HH_

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>
#include <list>
#include <include/myglobal_functions.hh>

#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Kernel/global_functions.h>
#include <list>
#include <CGAL/result_of.h>

#include "list_seg_index_unique.hh"
#include "modify_point_3.hh"
#include "my_line_3_intersection.hh"

#define VERY_VERBOSE 0
#define VERY_VERY_VERBOSE 0

#define USE_MY_INTERSECTION_L3L3 0

struct Output_ll_Int{
	bool is_coincident_, no_intersect_;
	Point_3 pts_;
};

struct Pts_Distance
{
	Point_3 pts_;
	double dis_;
};

bool compareDis(const Pts_Distance &lhs, const Pts_Distance &rhs)
{
	return lhs.dis_ < rhs.dis_;
}

namespace CGAL{


template <typename K, typename NT, typename Ell_>
std::list< CGAL::Point_3<K> > divide_intersection_segment_2poly(std::list< CGAL::Point_3<K> > &pl, Ell_ & Ell1, Ell_ & Ell2){

	typedef typename K::Point_3                    		Point_3;
	typedef typename K::Vector_3                   		Vector_3;
	typedef typename K::Line_3                     		Line_3;
	typedef typename K::Segment_3                  		Segment_3;
	typedef typename K::Plane_3                  		Plane_3;
	typedef typename K::Intersect_3 					Intersect_3;
	typedef CGAL::Surface_mesh<Point_3>        			Surface_mesh_;
	typedef typename Surface_mesh_::Vertex_index 		vertex_descriptor;

	//typedef CGAL::Exact_predicates_inexact_constructions_kernel_MP_Float Kernel;

	typedef typename std::list<Line_3> 					LineList;
	typedef typename std::list<Line_3>::iterator 		LineList_Iter;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Segment_3) >::type result_type;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Line_3) >::type result_type_ll;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Plane_3) >::type result_type_lp;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Plane_3, Plane_3) >::type result_type_pp;

	typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
	typedef EPECK::Point_3                    			Point_3_epeck;
	typedef EPECK::Line_3                     			Line_3_epeck;
	typedef EPECK::Segment_3                  			Segment_3_epeck;
	typedef EPECK::Plane_3                  			Plane_3_epeck;
	typedef EPECK::Intersect_3 							Intersect_3_epeck;
	typedef CGAL::Surface_mesh<Point_3_epeck>       	Surface_mesh_epeck;
	typedef typename Surface_mesh_epeck::Vertex_index 	vertex_descriptor_epeck;

	typedef CGAL::Cartesian_converter<K,EPECK > Converter;

	NT eps_(1E-16);

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

	int nameEp1 (Ell1.E_name), nameEp2 (Ell2.E_name);

	/*
	typedef typename CGAL::Point_3<K>                    	Point_3;
	typedef typename CGAL::Line_3<K>                      	Line_3;
	typedef typename CGAL::Segment_3<K>                   	Segment_3;

	typedef typename CGAL::Point_3<K>                    	Point_3;
	typedef typename CGAL::Vector_3<K>                     	Vector_3;
	typedef typename CGAL::Line_3<K>                     	Line_3;
	typedef typename CGAL::Segment_3<K>                 	Segment_3;
	typedef typename CGAL::Plane_3<K>                  		Plane_3;
	typedef typename K::Intersect_3   						Intersect_3;

	typedef typename std::list<Line_3> 					LineList;
	typedef typename std::list<Line_3>::iterator 		LineList_Iter;
	typedef typename CGAL::cpp11::result_of< Intersect_3(Line_3, Segment_3) >::type result_type;
	 */

	std::cout << "Em An xinh dep tuyet tran." << std::endl;

	LineList 	intlinelist;
	Line_3 		intline;

	std::list<int> 			enameList;

	// List of points
	std::list<Point_3> 		PointList, finalPointList;
	Point_3 pt_begin, pt_end;
	pt_begin = *pl.begin();
	Pts_Iterator pt_end_it = pl.begin();
	std::advance(pt_end_it,pl.size()-1);
	pt_end = *pt_end_it;

	PointList = pl;

	int inum(0);
	for (Pts_Iterator it = PointList.begin(); it != PointList.end(); it++){
		std::cout << "PointList [" << inum << "] = " << *it << std::endl;
		inum++;
	}

	double p0_0, p0_1, p0_2, p1_0, p1_1, p1_2, eps(CGAL::to_double(eps_));

	std::string fname;
	fname = "intersections/intersection_line_" + CGAL::int2str(nameEp1) + ".dat";
	std::ifstream    infile(fname.c_str());
	if (VERY_VERBOSE == 1){
		std::cout << "Check intersection with lines from intersections/intersection_line_" << CGAL::int2str(nameEp1) << ".dat" << std::endl;
	}

	std::string line;
	//std::istringstream ist(line);

	while (!infile.eof()){
		int num_ell;
		std::getline(infile, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell >> std::setprecision(30) >> p0_0 >> p0_1 >> p0_2 >> p1_0 >> p1_1 >> p1_2;

		if (CGAL::abs(p0_0) < eps){
			p0_0 = 0.E+0;
		}
		if (CGAL::abs(p0_1) < eps){
			p0_1 = 0.E+0;
		}
		if (CGAL::abs(p0_2) < eps){
			p0_2 = 0.E+0;
		}
		if (CGAL::abs(p1_0) < eps){
			p1_0 = 0.E+0;
		}
		if (CGAL::abs(p1_1) < eps){
			p1_1 = 0.E+0;
		}
		if (CGAL::abs(p1_2) < eps){
			p1_2 = 0.E+0;
		}

		NT p0_0_NT(p0_0), p0_1_NT(p0_1), p0_2_NT(p0_2), p1_0_NT(p1_0), p1_1_NT(p1_1), p1_2_NT(p1_2);

		if (num_ell != nameEp2){
			enameList.push_back(num_ell);
			//intlinelist.push_back( Line_3( Point_3(p0_0, p0_1, p0_2), Point_3(p1_0, p1_1, p1_2)) );
			intlinelist.push_back( Line_3( Point_3(p0_0_NT, p0_1_NT, p0_2_NT), Point_3(p1_0_NT, p1_1_NT, p1_2_NT)) );
		}
	}
	infile.close();

	fname = "intersections/intersection_line_" + CGAL::int2str(nameEp2) + ".dat";
	std::ifstream    infile_(fname.c_str());
	if (VERY_VERBOSE == 1){
		std::cout << "Check intersection with lines from intersections/intersection_line_" << CGAL::int2str(nameEp2) << ".dat" << std::endl;
	}

	while (!infile_.eof()){
		int num_ell;
		std::getline(infile_, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> num_ell >> std::setprecision(30) >> p0_0 >> p0_1 >> p0_2 >> p1_0 >> p1_1 >> p1_2;

		if (CGAL::abs(p0_0) < eps){
			p0_0 = 0.E+0;
		}
		if (CGAL::abs(p0_1) < eps){
			p0_1 = 0.E+0;
		}
		if (CGAL::abs(p0_2) < eps){
			p0_2 = 0.E+0;
		}
		if (CGAL::abs(p1_0) < eps){
			p1_0 = 0.E+0;
		}
		if (CGAL::abs(p1_1) < eps){
			p1_1 = 0.E+0;
		}
		if (CGAL::abs(p1_2) < eps){
			p1_2 = 0.E+0;
		}

		NT p0_0_NT(p0_0), p0_1_NT(p0_1), p0_2_NT(p0_2), p1_0_NT(p1_0), p1_1_NT(p1_1), p1_2_NT(p1_2);

		if (num_ell != nameEp1){
			enameList.push_back(num_ell);
			//intlinelist.push_back( Line_3( Point_3(p0_0, p0_1, p0_2), Point_3(p1_0, p1_1, p1_2)) );
			intlinelist.push_back( Line_3( Point_3(p0_0_NT, p0_1_NT, p0_2_NT), Point_3(p1_0_NT, p1_1_NT, p1_2_NT)) );
		}
	}
	infile_.close();

	// NTD 21/09/2016: Imprecision of intersection between line consisting of close-to-zero coordinate point
	Pts_Iterator itp1 = pl.begin();
	Pts_Iterator itp2 = pl.begin();
	std::advance(itp2, pl.size()-1);

	Point_3 p1(*itp1), p2(*itp2);

	Segment_3 test_seg = Segment_3(p1,p2);
	//Segment_3 test_seg = Segment_3(CGAL::zero_p3<Point_3,NT>(p1, eps), CGAL::zero_p3<Point_3,NT>(p2, eps));

	Point_3 pt_1_polygon1, pt_2_polygon1, pt_3_polygon1;
	Point_3 pt_1_polygon2, pt_2_polygon2, pt_3_polygon2;
	{
		unsigned int i = 0, end = Ell1.E_mesh.number_of_vertices();
		for( ; i < end; ++i) {
			vertex_descriptor vh(i);
			if (i==0)
				pt_1_polygon1 = Ell1.E_mesh.point(vh);
			else if(i==round(Ell1.E_mesh.num_vertices()/4))
				pt_2_polygon1 = Ell1.E_mesh.point(vh);
			else if (i==(Ell1.E_mesh.num_vertices()/2))
				pt_3_polygon1 = Ell1.E_mesh.point(vh);
		}
	}
	{
		unsigned int i = 0, end = Ell2.E_mesh.number_of_vertices();
		for( ; i < end; ++i) {
			vertex_descriptor vh(i);
			if (i==0)
				pt_1_polygon2 = Ell2.E_mesh.point(vh);
			else if(i==round(Ell2.E_mesh.num_vertices()/4))
				pt_2_polygon2 = Ell2.E_mesh.point(vh);
			else if (i==(Ell2.E_mesh.num_vertices()/2))
				pt_3_polygon2 = Ell2.E_mesh.point(vh);
		}
	}


	Plane_3 test_plane1 = Plane_3(pt_1_polygon1, pt_2_polygon1, pt_3_polygon1);
	Plane_3 test_plane2 = Plane_3(pt_1_polygon2, pt_2_polygon2, pt_3_polygon2);

	// ==============================================================================================
	// EPECK kernel

	//CGAL::Plane_3<EPECK> test_plane_epec1 = CGAL::Cartesian_converter<K, EPECK, NT_MP> (test_plane1);
	//CGAL::Plane_3<EPECK> test_plane_epec2(test_plane2);


	Surface_mesh_epeck poly_1_epeck, poly_2_epeck;

	Convert_sm()(Ell1.E_mesh, poly_1_epeck);
	Convert_sm()(Ell2.E_mesh, poly_2_epeck);

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


	if (VERY_VERBOSE == 1){
		CGAL::cpp11::result_of<EPECK::Intersect_3(EPECK::Plane_3, EPECK::Plane_3)>::type
				result_epec_pp_1 = CGAL::intersection(test_plane1_epeck, test_plane2_epeck);

		if (const Line_3_epeck* l = boost::get<Line_3_epeck >(&*result_epec_pp_1)) {
			std::cout << "EPECK: Found intersection between the plane " << std::setprecision(30) << test_plane1 << " and the plane " << std::setprecision(30) << test_plane2 << std::endl;
			std::cout << std::setprecision(30) << *l << std::endl;
			std::cout << "\n" << std::endl;
		}
	}

	// ==============================================================================================

	if (VERY_VERBOSE == 1){
		result_type_pp result_pp1 = CGAL::intersection(test_plane1, test_plane2);
		if (const Line_3* l = boost::get<Line_3 >(&*result_pp1)) {
			std::cout << "Found intersection between the plane " << std::setprecision(30) << test_plane1 << " and the plane " << std::setprecision(30) << test_plane2 << std::endl;
			std::cout << std::setprecision(30) << *l << std::endl;
			std::cout << "\n" << std::endl;
		}
	}



	inum = 0;
	for (LineList_Iter iter = intlinelist.begin(); iter != intlinelist.end(); iter++){
		Int_Iterator enameList_iter = enameList.begin();

		Line_3 testline = *iter;

		/*
		if (VERY_VERBOSE == 1){
			std::cout << "Checking intersection between the line " << std::setprecision(30) << testline << " and the plane " << std::setprecision(30) << test_plane1 << std::endl;
			result_type_lp result_lp1 = CGAL::intersection(testline, test_plane1);

			if (const Point_3* p = boost::get<Point_3 >(&*result_lp1)) {
				if (VERY_VERBOSE == 1){
					std::cout << "Intersection between the line " << testline << " and the plane " << test_plane1 << " found: "<< *p << std::endl;
				}

				std::cout << "P_1 : " << pt_begin << std::endl;
				std::cout << "P_2 : " << *p << std::endl;
				std::cout << "P_3 : " << pt_end << std::endl;
				std::cout << "Collinear : " << CGAL::collinear(pt_begin, pt_end, *p) << std::endl;


				if ( (CGAL::squared_distance(pt_begin, *p) < CGAL::squared_distance(pt_begin, pt_end)) &&
						(Vector_3(pt_begin, *p) * Vector_3(pt_begin, pt_end) >=0) ){
					PointList.push_back(*p);
				}

			}
		}
		 */

		if (USE_MY_INTERSECTION_L3L3 == 1){
			std::cout << "Checking intersection between the line " << std::setprecision(30) << testline << " and the plane " << std::setprecision(30) << test_plane1 << std::endl;
			Output_ll_Int output_ll_int;
			Line_3 test_seg_supp_line = test_seg.supporting_line();
			output_ll_int = my_line_3_intersection<K, NT, Output_ll_Int>(testline, test_seg_supp_line, eps_);

			std::cout << "output_ll_int.no_intersect_ =  "<< output_ll_int.no_intersect_ << std::endl;
			if (!output_ll_int.no_intersect_ && !output_ll_int.is_coincident_) {

				std::cout << "P_1 : " << pt_begin << std::endl;
				std::cout << "P_2 : " << output_ll_int.pts_ << std::endl;
				std::cout << "P_3 : " << pt_end << std::endl;
				std::cout << "Collinear : " << CGAL::collinear(pt_begin, pt_end, output_ll_int.pts_) << std::endl;

				if ( (CGAL::squared_distance(pt_begin, output_ll_int.pts_) < CGAL::squared_distance(pt_begin, pt_end)) &&
						(Vector_3(pt_begin, output_ll_int.pts_) * Vector_3(pt_begin, pt_end) >=0) ){
					PointList.push_back(output_ll_int.pts_);
				}
			}
		}

		if (VERY_VERBOSE == 1){
			std::cout << "Checking intersection between the segment " << std::setprecision(30) << test_seg << " and the line " << std::setprecision(30) << testline << std::endl;
		}

		result_type result_sl = CGAL::intersection(testline, test_seg);
		//CGAL::cpp11::result_of<K::Intersect_3(Line_3, Segment_3)>::type result_sl = CGAL::intersection(testline, test_seg);

		if (result_sl) {

			if (const Point_3* p = boost::get<Point_3 >(&*result_sl)) {
				if (VERY_VERBOSE == 1){
					std::cout << "Intersection between the segment " << test_seg << " and the line " << testline << " found: "<< *p << std::endl;
				}
				PointList.push_back(*p);
			}
		}

		inum++;
	}

	// Sorting the point by the distance from p1
	Pts_Distance pts_dis;
	std::list<Pts_Distance> list_pts_dis;

	inum = 0;
	for (Pts_Iterator it = PointList.begin(); it != PointList.end(); it++){
		if (VERY_VERBOSE == 1){
			std::cout << "PointList [" << inum << "] = " << *it << std::endl;
		}
		inum++;

		pts_dis.pts_ = *it;
		pts_dis.dis_ =  CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(*it, p1)));
		list_pts_dis.push_back(pts_dis);
	}


	list_pts_dis.sort(&compareDis);

	for (std::list<Pts_Distance>::iterator it = list_pts_dis.begin(); it !=list_pts_dis.end(); it++){
		finalPointList.push_back( (*it).pts_ );
	}

	/* persons list initialization */
	//persons.sort(boost::bind(&Person::age, _1) < boost::bind(&Person::age, _2));

	// Remove duplicate points
	finalPointList.unique();
	//sorted_list_points_unique <Point_3, NT_MP>( finalPointList, eps_);

	// Test
	if (VERY_VERY_VERBOSE == 1){
		/*
		NT testP1_x(-0.18863), testP1_y(-0.41137), testP1_z(-0.4);
		NT testP2_x(0.139037), testP2_y(-0.739037), testP2_z(-0.4);
		NT testP3_x(0.0), testP3_y(-0.), testP3_z(-0.4);
		NT testP4_x(0.766951), testP4_y(-15.339), testP4_z(-0.4);
		 */
		/*
		NT testP1_x(0.499608), testP1_y(-0.174507), testP1_z(-0.811194);
		NT testP2_x(0.861021), testP2_y(-0.182516), testP2_z(-1.5);
		//NT testP3_x(1.564790), testP3_y(0.9155960),	testP3_z(-0);
		//NT testP4_x(1.426760), testP4_y(0.7666970), testP4_z(-0.124605);
		NT testP3_x(-19.6173), testP3_y(13.9962), testP3_z(-0);
		NT testP4_x(-17.7925), testP4_y(12.7141), testP4_z(-0.0826208);
		 */

		NT testP1_x(0.499608), testP1_y(-0.174507), testP1_z(-0.811194);
		NT testP2_x(0.861021), testP2_y(-0.182516), testP2_z(-1.500008);
		NT testP3_x(1.56479188356251386338E+00), testP3_y(9.15595602576112277404E-01),	testP3_z(-0.00000000000000000000E+00);
		//NT testP4_x(1.42676123432908008581E+00), testP4_y(7.66697416639087636625E-01), testP4_z(-1.24604740194022789446E-01);
		NT testP4_x(0.1844854), testP4_y(- 0.5733863), testP4_z(- 1.2460474);


		//NT testP3_x(-19.), testP3_y(13.9962), testP3_z(-0);
		//NT testP4_x(-17.7925), testP4_y(12.7141), testP4_z(-0.0826208);


		Point_3 testP1( Point_3(testP1_x, testP1_y, testP1_z) );
		Point_3 testP2( Point_3(testP2_x, testP2_y, testP2_z) );
		Point_3 testP3( Point_3(testP3_x, testP3_y, testP3_z) );
		Point_3 testP4( Point_3(testP4_x, testP4_y, testP4_z) );

		Segment_3 test_seg (Segment_3(testP1, testP2));
		Line_3 testline (Line_3(testP3, testP4));
		std::cout << "Remark: Checking intersection between the segment " << test_seg << " and the line " << testline << std::endl;

		//CGAL::cpp11::result_of<Kernel::Intersect_3(Line_3, Segment_3)>::type result_sl = CGAL::intersection(testline, test_seg);
		result_type result_sl = CGAL::intersection(testline, test_seg);

		if (result_sl) {

			if (const Point_3* p = boost::get<Point_3 >(&*result_sl)) {
				std::cout << "Intersection check: Yes" << std::endl;
				std::cout << "Remark: Intersection between the segment " << test_seg << " and the line " << testline << " found: "<< *p << std::endl;
			}
		}


		float p1_x(CGAL::to_double(p1.x())), p1_y(CGAL::to_double(p1.y())), p1_z(CGAL::to_double(p1.z()));
		float p2_x(CGAL::to_double(p2.x())), p2_y(CGAL::to_double(p2.y())), p2_z(CGAL::to_double(p2.z()));

		/*
		NT p1_x(p1.x()), p1_y(p1.y()), p1_z(p1.z());
		NT p2_x(p2.x()), p2_y(p2.y()), p2_z(p2.z());
		 */
		/*
		Point_3 p1_ ( Point_3(p1.x(), p1.y(), p1.z()) );
		Point_3 p2_ ( Point_3(p2.x(), p2.y(), p2.z()) );
		 */
		Point_3 p1_ ( Point_3(p1_x, p1_y, p1_z) );
		Point_3 p2_ ( Point_3(p2_x, p2_y, p2_z) );

		std::cout << std::setprecision(12) << "testP1.x() = " << testP1.x() << " ; " << "testP1.y() = " << testP1.y() << " ; " << "testP1.z() = " << testP1.z() << std::endl;
		std::cout << std::setprecision(12) << "testP2.x() = " << testP2.x() << " ; " << "testP2.y() = " << testP1.y() << " ; " << "testP2.z() = " << testP1.z() << std::endl;
		std::cout << std::setprecision(12) << "p1_x = " << p1_x << " ; " << "p1_y = " << p1_y << " ; " << "p1_z = " << p1_z << std::endl;
		std::cout << std::setprecision(12) << "p2_x = " << p2_x << " ; " << "p2_y = " << p2_y << " ; " << "p2_z = " << p2_z << std::endl;

		Segment_3 test_seg_1 (Segment_3(p1_, p2_));
		std::cout << "Remark: Checking intersection between the segment " << std::setprecision(30) << test_seg_1 << " and the line " << std::setprecision(30) << testline << std::endl;

		//CGAL::cpp11::result_of<Kernel::Intersect_3(Line_3, Segment_3)>::type result_sl = CGAL::intersection(testline, test_seg);
		result_type result_sl_1 = CGAL::intersection(testline, test_seg_1);

		if (result_sl_1) {
			std::cout << "Intersection check: Yes" << std::endl;
			if (const Point_3* p = boost::get<Point_3 >(&*result_sl_1)) {
				std::cout << "Remark: Intersection between the segment " << test_seg_1 << " and the line " << testline << " found: "<< *p << std::endl;
			}
		}

		/*
		//Line_3 test_line_1 (Line_3(p1_, p2_));
		Line_3 test_line_1 (*intlinelist.begin());
		//Line_3 test_line_2 (*boost::next(boost::next(intlinelist.begin())));
		Line_3 test_line_3 (*boost::next(intlinelist.begin()));
		std::cout << "Remark: Checking intersection between the line " << std::setprecision(30) << test_line_1 << " and the line " << std::setprecision(30) << test_line_3 << std::endl;

		//result_type_ll result_ll_1 = CGAL::intersection(test_line_1, test_line_2);
		CGAL::cpp11::result_of<Kernel::Intersect_3(Kernel::Line_3, Kernel::Line_3)>::type result_ll_1 = CGAL::intersection(test_line_1, test_line_3);

		if (result_ll_1) {
			std::cout << "Intersection check: Yes" << std::endl;
			if (const Point_3* p = boost::get<Point_3 >(&*result_sl_1)) {
				std::cout << "Remark: Intersection between the segment " << test_seg_1 << " and the line " << test_line_2 << " found: "<< *p << std::endl;
			}
		}
		 */
	}

	return finalPointList;
}
} // end of namespace


#endif /* INCLUDE_DIVIDE_INTERSECTION_SEGMENT_HH_ */
