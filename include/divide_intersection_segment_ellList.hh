/*
 * divide_intersection_segment_ellList.hh
 *
 *  Created on: Sep 30, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_DIVIDE_INTERSECTION_SEGMENT_ELLLIST_HH_
#define INCLUDE_DIVIDE_INTERSECTION_SEGMENT_ELLLIST_HH_


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

template <typename K, typename NT, typename Ell_, typename Output_IPL>
std::list< CGAL::Point_3<K> > divide_intersection_segment_2poly_ellList(std::list< CGAL::Point_3<K> > &pl,
		std::list<Ell_> & ellList,
		Ell_ &Ell1,
		Ell_ &Ell2){

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

	// Get ellipses
	typedef typename std::list<Ell_>::const_iterator Ell_Iterator;

	/*
	Ell_ Ell1, Ell2;
	Ell_Iterator ell_1_Iter(ellList.begin()), ell_2_Iter(ellList.begin());
	std:advance(ell_1_Iter, nameEp1);
	std:advance(ell_2_Iter, nameEp2);
	Ell1 = *ell_1_Iter;
	Ell2 = *ell_1_Iter;
	 */

	//int nameEp1(Ell1.E_name), nameEp2(Ell2.E_name);
	int nameEp1(Ell1.Parent_E_name - 1), nameEp2(Ell2.Parent_E_name - 1);

	//std::cout << "Em An xinh dep tuyet tran." << std::endl;

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
	if (VERY_VERBOSE == 1){
		for (Pts_Iterator it = PointList.begin(); it != PointList.end(); it++){
			std::cout << "PointList [" << inum << "] = " << *it << std::endl;
			inum++;
		}
	}

	// =====================================================================================================
	//                                      Get the list of intersection line
	// =====================================================================================================

	// All intersection lines which lie on ellipse 1
	if (VERY_VERBOSE == 1){
		std::cout << "Checking intersection between the segment with all intersection lines of ellipse " <<  CGAL::int2str(nameEp1) << std::endl;
	}

	std::string fname;
	fname = "intersections/intersection_orig_" + CGAL::int2str(nameEp1) + ".dat";
	std::ifstream    infile(fname.c_str());
	std::string line;

	int num_ell(0);
	while (!infile.eof() && num_ell < ellList.size()){
		std::getline(infile, line);

		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}

		std::istringstream ist(line);
		ist.str(line);
		ist.clear();

		int check_intersection;
		ist >> check_intersection;
		if (VERY_VERBOSE == 1){
			std::cout << "num_ell = " << num_ell << ", check_intersection = " << check_intersection << std::endl;
		}

		if (num_ell != nameEp1 && num_ell != nameEp2 && check_intersection > 0){
			enameList.push_back(num_ell);

			Ell_Iterator ell_test_Iter(ellList.begin());
			std::advance(ell_test_Iter, num_ell);
			Ell_ ell_test = *ell_test_Iter;
			Output_IPL out_ipl;

			out_ipl = CGAL::intersections_of_two_polys_exp_line<K, Output_IPL>( Ell1.E_mesh, ell_test.E_mesh);

			Line_3 Int_Line_ = out_ipl.intersection_line_;
			if (VERY_VERBOSE == 1){
				std::cout << "Intersection line : " <<  Int_Line_<< std::endl;
			}
			//result_type result_ls = CGAL::intersection(test_seg, Int_Line_);

			intlinelist.push_back(Int_Line_);
		}
		num_ell ++;
	}
	infile.close();

	// --------------------------------------------------------------------------------------------------------
	// All intersection lines which lie on ellipse 2
	if (VERY_VERBOSE == 1){
		std::cout << "Checking intersection between the segment with all intersection lines of ellipse "<<  CGAL::int2str(nameEp2) << std::endl;
	}

	fname = "intersections/intersection_orig_" + CGAL::int2str(nameEp2) + ".dat";
	std::ifstream    infile_(fname.c_str());

	num_ell = 0;
	while (!infile_.eof() && num_ell < ellList.size()){
		std::getline(infile_, line);
		if (VERY_VERBOSE == 1){
			std::cout << "Line : " << line << std::endl;
		}

		int check_intersection;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> check_intersection;

		if (VERY_VERBOSE == 1){
			std::cout << "num_ell = " << num_ell << ", check_intersection = " << check_intersection << std::endl;
		}

		if (num_ell != nameEp1 && num_ell != nameEp2 && check_intersection > 0){
			enameList.push_back(num_ell);

			Ell_Iterator ell_test_Iter(ellList.begin());
			std::advance(ell_test_Iter, num_ell);
			Ell_ ell_test = *ell_test_Iter;
			Output_IPL out_ipl;

			out_ipl = CGAL::intersections_of_two_polys_exp_line<K, Output_IPL>( Ell1.E_mesh, ell_test.E_mesh);

			Line_3 Int_Line_ = out_ipl.intersection_line_;
			//result_type result_ls = CGAL::intersection(test_seg, Int_Line_);

			intlinelist.push_back(Int_Line_);
		}
		num_ell ++;
	}
	infile_.close();

	// =====================================================================================================
	//                                      Get the new point list
	// =====================================================================================================

	// NTD 21/09/2016: Imprecision of intersection between line consisting of close-to-zero coordinate point
	Pts_Iterator itp1 = pl.begin();
	Pts_Iterator itp2 = pl.begin();
	std::advance(itp2, pl.size()-1);

	Point_3 p1(*itp1), p2(*itp2);

	Segment_3 test_seg = Segment_3(p1,p2);
	//Segment_3 test_seg = Segment_3(CGAL::zero_p3<Point_3,NT>(p1, eps), CGAL::zero_p3<Point_3,NT>(p2, eps));

	inum = 0;
	for (LineList_Iter iter = intlinelist.begin(); iter != intlinelist.end(); iter++){
		Int_Iterator enameList_iter = enameList.begin();

		Line_3 testline = *iter;

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

	// Remove duplicate points
	finalPointList.unique();
	//sorted_list_points_unique <Point_3, NT_MP>( finalPointList, eps_);

	// NTD_01Fev2017
	//NT eps_filter(1E-4);
	//sorted_list_points_unique <Point_3, NT_MP>( finalPointList, eps_filter);
	//sorted_list_points <Point_3, NT_MP>( finalPointList, pl, eps_filter, eps_);
	// End NTD_01Fev2017

	return finalPointList;
}
} // end of namespace

#endif /* INCLUDE_DIVIDE_INTERSECTION_SEGMENT_ELLLIST_HH_ */
