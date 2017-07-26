/*
 * cutting_off_first_polygon_rmpts.hh
 *
 *  Created on: Oct 10, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CUTTING_OFF_FIRST_POLYGON_RMPTS_HH_
#define INCLUDE_CUTTING_OFF_FIRST_POLYGON_RMPTS_HH_


// Add some boost files
#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/include/deref.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/include/begin.hpp>

// Intersection
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>
#include <include/points_are_close.hh>

// Order of the intersection points
#include "include/order_intersectionpoints_3d.hh"

// Function to check if all elements of boolean list are true
#include "include/All_bool_true.hh"

// Minimum & maximum of an integer list
#include "include/Min_Int_List.hh"
#include "include/Max_Int_List.hh"

// Minimum & maximum of a list of numbers
#include "include/Min_Max_List.hh"

// Transpose a list
#include "include/transpose_a_list.hh"

#include "include/remove_close_points.hh"

//Intersection between two planes on which two polygons lie
#include "list_seg_index_unique.hh"
//#include "plane_plane_intersection_convert2MP.hh"
#include "plane_plane_intersection.hh"
//#include "Triangulation_MP_Float.hh"

//#include "divide_intersection_segment.hh"
#include "divide_intersection_segment_ellList.hh"

//#include <write_off_to_vtk_xml_file.h>
#include "Output_off_vtu.hh"

// check_conformity_poly_list
//#include "check_conformity_poly_list.hh"
#include "check_conformity_newpoly.hh"
#include "libmessages.hh"

#ifndef DEBUG_NTD
#define DEBUG_NTD 0
#endif

//#ifndef VERY_VERBOSE
#define VERY_VERBOSE 0
//#endif

#ifndef VERY_VERY_VERBOSE
#define VERY_VERY_VERBOSE 0
#endif

#ifndef PRINT_VTU
#define PRINT_VTU 0
#endif

namespace CGAL{

template <typename Kernel, typename Output_Child_Polys, typename Ell_, typename IP_Ell, typename Output_IPL, typename NT>
Output_Child_Polys cutting_off_first_polygon_rmpts(std::list<Ell_> &ell_list, Ell_ &cEll_1, Ell_ cEll_2, double &target_edge_length, NT &tol_length, bool bool_mergepts, double &mergeClosePointsRelCritLength){

	Output_Child_Polys output_;

	typedef typename std::list<IP_Ell> IPE_list;
	typedef typename CGAL::Point_3<Kernel>                    	Point_3;
	typedef typename CGAL::Vector_3<Kernel>                    	Vector_3;
	typedef typename CGAL::Line_3<Kernel>                      	Line_3;
	typedef typename CGAL::Segment_3<Kernel>                   	Segment_3;
	typedef typename CGAL::Plane_3<Kernel>                   	Plane_3;
	typedef typename CGAL::Surface_mesh<Point_3>				SurfaceMesh_;


	typedef typename SurfaceMesh_::Vertex_index 		   vertex_descriptor;
	typedef typename SurfaceMesh_::Vertex_iterator vertex_iterator;

	typedef typename IPE_list::const_iterator IPE_Iterator;
	typedef typename std::list<Point_3>::iterator Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator Vector_Iterator;
	typedef std::list<int>::iterator Int_Iterator;
	typedef std::list<char>::iterator Char_Iterator;

	IP_Ell IP1_Ell1, IP2_Ell1, IP1_Ell2, IP2_Ell2;

	SurfaceMesh_ 	polygon1(cEll_1.E_mesh), polygon2(cEll_2.E_mesh);
	int 			nameP1(cEll_1.Parent_E_name),  nameP2(cEll_2.Parent_E_name);

	std::string fname;

	int num_intersec_1(0), num_intersec_2(0);
	std::list<Point_3> intersection_points_1, intersection_points_2;
	std::list<Point_3> intersectionline;

	std::list<int> list_seg_interst1, list_seg_interst2;
	bool no_polygons_interst(false);
	Output_IPL output_ip;

	/* Intersection between two polygons */
	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IPL> (polygon1, polygon2);

	num_intersec_1 = output_ip.num_intersec_1_;
	num_intersec_2 = output_ip.num_intersec_2_;
	intersection_points_1 = output_ip.intersection_points_1_;
	intersection_points_2 = output_ip.intersection_points_2_;
	list_seg_interst1 = output_ip.list_seg_interst1_;
	list_seg_interst2 = output_ip.list_seg_interst2_;

	//NT tol_length = 0.1 * target_edge_length;
	/*
	const std::string off_file (filename_ + ".off");
	const char* filename = off_file.c_str();
	std::ifstream    infile(filename);
	 */

	// Cut off the polygons
	SurfaceMesh_ poly_dummy, poly_1_child_1, poly_1_child_2;

	if(VERY_VERBOSE == 1){
		for (Pts_Iterator it = intersection_points_1.begin(); it != intersection_points_1.end(); it++){
			std::cout << "intersection_points_1 : " << *it << std::endl;
		}
		std::cout << "intersection_points_2.size() : " << intersection_points_2.size() << std::endl;
		for (Pts_Iterator it = intersection_points_2.begin(); it != intersection_points_2.end(); it++){
			std::cout << "intersection_points_2 : " << *it << std::endl;

		}
		std::cout << "intersection_points_1.size() + intersection_points_2.size() : " << intersection_points_1.size() + intersection_points_2.size() << std::endl;

		std::cout << "num_intersec_1 = " << num_intersec_1 << std::endl;
		std::cout << "num_intersec_2 = " << num_intersec_1 << std::endl;
	}


	if (intersection_points_1.size() + intersection_points_2.size() < 4){
		if(VERY_VERBOSE == 1){
			std::cout << "Two polygons do not intersect: num_intersec_1 + num_intersec_2 < 4." << std::endl;
			std::cout << " " << std::endl;
		}
		no_polygons_interst = true;
	}
	else{

		// Be careful when intersection_points_1.size() + intersection_points_2.size() < 4
		Pts_Iterator ell1_interst_it_end = intersection_points_1.begin();
		std::advance(ell1_interst_it_end,1);

		Pts_Iterator ell2_interst_it_end = intersection_points_2.begin();
		std::advance(ell2_interst_it_end,1);

		IP1_Ell1.Ell_name_ = '1';
		IP2_Ell1.Ell_name_ = '1';
		IP1_Ell2.Ell_name_ = '2';
		IP2_Ell2.Ell_name_ = '2';


		IP1_Ell1.pts_ = *intersection_points_1.begin();
		IP2_Ell1.pts_ = *ell1_interst_it_end;
		IP1_Ell2.pts_ = *intersection_points_2.begin();
		IP2_Ell2.pts_ = *ell2_interst_it_end;

		std::list<IP_Ell> list_IP_Ell, list_IP_Ell_output;

		list_IP_Ell.push_back(IP1_Ell1);
		list_IP_Ell.push_back(IP2_Ell1);
		list_IP_Ell.push_back(IP1_Ell2);
		list_IP_Ell.push_back(IP2_Ell2);

		std::list<Point_3> AddPoints2NewPoly_1, AddPoints2NewPoly_2, AddPoints2NewPoly_1_orig, AddPoints2NewPoly_2_orig, AddPoints2IntersectionLine;
		std::list<Point_3> AddPoints2NewPoly_1_new,AddPoints2NewPoly_1_dummy;

		//std::cout << "bool_mergepts = " << bool_mergepts << std::endl;

		if (!bool_mergepts){
			list_IP_Ell_output = CGAL::Order_IntersectionPoins3D<Kernel, IP_Ell>(list_IP_Ell);
		}
		else{
			/* NTD 18/11/2016: Merge too close points : distance between two pts is less than 5% the target edge length*/
			//list_IP_Ell_output = CGAL::Order_IntersectionPoins3D_merge_close_pts<Kernel, IP_Ell>
			//				(list_IP_Ell, 0.05* CGAL::min(cEll_1.E_target_edge_length, cEll_2.E_target_edge_length));


			/* NTD 05/01/2017: Merge too close points : distance between two pts is less than 10% the target edge length*/
			list_IP_Ell_output = CGAL::Order_IntersectionPoins3D_merge_close_pts<Kernel, IP_Ell>
							(list_IP_Ell, mergeClosePointsRelCritLength * CGAL::min(cEll_1.E_target_edge_length, cEll_2.E_target_edge_length));
		}

		if ( (*list_IP_Ell_output.begin()).Ell_name_ == (*(boost::next(list_IP_Ell_output.begin()))).Ell_name_ ){
			if(VERY_VERBOSE == 1){
				std::cout << "Two polygons do not intersect: two intersection segments have no common point." << std::endl;
				std::cout << " " << std::endl;
			}
			no_polygons_interst = true;
		}
		else{
			switch ((*list_IP_Ell_output.begin()).Ell_name_){
			case '1':
			{
				switch ( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).Ell_name_){
				case '1':
				{
					AddPoints2NewPoly_1.push_back( (*list_IP_Ell_output.begin()).pts_); // First element
					AddPoints2NewPoly_1.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element

					AddPoints2NewPoly_2.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_); // Fourth element
				}
				break;

				case '2':
				{
					AddPoints2NewPoly_1.push_back( (*list_IP_Ell_output.begin()).pts_); // First element
					AddPoints2NewPoly_1.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_); // Fourth element

					AddPoints2NewPoly_2.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
				}
				break;
				}
			}
			break;

			case '2':
			{
				switch ( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).Ell_name_){

				case '1':
				{
					AddPoints2NewPoly_2.push_back( (*list_IP_Ell_output.begin()).pts_); // First element
					AddPoints2NewPoly_2.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_); // Fourth element

					AddPoints2NewPoly_1.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
				}
				break;

				case '2':
				{
					AddPoints2NewPoly_2.push_back( (*list_IP_Ell_output.begin()).pts_); // First element
					AddPoints2NewPoly_2.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_2.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element

					AddPoints2NewPoly_1.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_); // Second element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_); // Third element
					AddPoints2NewPoly_1.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_); // Fourth element
				}
				break;
				}
			}
			break;
			}

			if (VERY_VERBOSE == 1){
				for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
					std::cout << "AddPoints2NewPoly_1 = " << *iter << std::endl;
				}
				std::cout << " " << std::endl;

				for (Pts_Iterator iter=AddPoints2NewPoly_2.begin(); iter!=AddPoints2NewPoly_2.end(); iter++){
					std::cout << "AddPoints2NewPoly_2 = " << *iter << std::endl;
				}
				std::cout << " " << std::endl;

				std::cout << "###########################" << std::endl;
				std::cout << "Cut off the polygons" << std::endl;
				std::cout << "###########################" << std::endl;
			}

			if (VERY_VERY_VERBOSE ==1){
				std::cout<< "Min_List<int>(list_seg_interst1) = " << Min_List<int>(list_seg_interst1) <<std::endl;
				std::cout<< "Max_List<int>(list_seg_interst1) = "<< Max_List<int>(list_seg_interst1) <<std::endl;

				std::cout<< "Min_List<int>(list_seg_interst2) = "<< Min_List<int>(list_seg_interst2) <<std::endl;
				std::cout<< "Max_List<int>(list_seg_interst2) = "<< Max_List<int>(list_seg_interst2) <<std::endl;
			}

			fname = "intersections/intersection_line_" + CGAL::int2str(nameP2) + ".dat";
			/*========================================================================================================*/

			if (VERY_VERBOSE ==1){
				std::cout << "*list_IP_Ell_output.begin()).Ell_name_ = " << (*list_IP_Ell_output.begin()).Ell_name_ << std::endl;
				std::cout << "*next_list_IP_Ell_output.begin()).Ell_name_ = " << (*(boost::next(list_IP_Ell_output.begin()))).Ell_name_ << std::endl;
				std::cout << "*next_next_list_IP_Ell_output.begin()).Ell_name_ = " << (*(boost::next(boost::next(list_IP_Ell_output.begin())))).Ell_name_ << std::endl;
				std::cout << "*next_next_next_list_IP_Ell_output.begin()).Ell_name_ = " << (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).Ell_name_ << std::endl;
			}


			// Polygon 1
			bool reverse_lst_1(false);
			std::list< Point_3 > point_list;

			// Child polygon poly_1_child_1
			{
				// Part 1: Add points from the old "father" polygon_1 to the new "child" poly_1_child_1
				unsigned int begin = Min_List<int>(list_seg_interst1);
				unsigned int end = Max_List<int>(list_seg_interst1);

				for( ; begin < end; ++begin) {
					vertex_descriptor vb_seg1(begin);
					vertex_descriptor ve_seg1(begin+1);

					poly_1_child_1.add_vertex(polygon1.point(ve_seg1));
				}


				// Update the list AddPoints2NewPoly_1 by adding intersection points with intersection lines on the second ellipse
				if ((*list_IP_Ell_output.begin()).Ell_name_ == '1'){
					switch ( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).Ell_name_){
					case '1':
						// Case --- 1 --- 2 --- 1 --- 2
					{
						if(VERY_VERBOSE == 1){
							std::cout << "This is the case: --- 1 --- 2 --- 1 --- 2" << std::endl;
						}
						//std::cout << "(*list_IP_Ell_output.begin()).pts_ = " << (*list_IP_Ell_output.begin()).pts_ << std::endl;

						point_list.push_back( (*list_IP_Ell_output.begin()).pts_ );
						point_list.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_ );

						AddPoints2NewPoly_1_new.push_back( (*list_IP_Ell_output.begin()).pts_ );
						//AddPoints2NewPoly_1_dummy = CGAL::divide_intersection_segment_2poly<Kernel,NT,Ell_>(point_list, cEll_1, cEll_2);
						AddPoints2NewPoly_1_dummy = CGAL::divide_intersection_segment_2poly_ellList<Kernel,NT,Ell_,Output_IPL>(point_list, ell_list, cEll_1, cEll_2);

						point_list.clear();

						int iter_1(0);
						for (Pts_Iterator it = AddPoints2NewPoly_1_dummy.begin(); it != AddPoints2NewPoly_1_dummy.end(); it++){
							if(VERY_VERBOSE == 1){
								std::cout << "01Fev2017_AddPoints2NewPoly_1_dummy [" << iter_1 <<  "] = " << *it << std::endl;
							}
							iter_1++;
							AddPoints2NewPoly_1_new.push_back(*it);
						}

						AddPoints2NewPoly_1 = AddPoints2NewPoly_1_new;
						AddPoints2NewPoly_1_dummy.clear();
						AddPoints2NewPoly_1_new.clear();
					}
					break;

					case '2':
						// Case --- 1 --- 2 --- 2 --- 1
					{
						point_list.push_back( (*list_IP_Ell_output.begin()).pts_ );
						point_list.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_ );
						//AddPoints2NewPoly_1_dummy = CGAL::divide_intersection_segment_2poly<Kernel,NT,Ell_>(point_list, cEll_1, cEll_2);
						AddPoints2NewPoly_1_dummy =
								CGAL::divide_intersection_segment_2poly_ellList<Kernel,NT,Ell_,Output_IPL>(point_list, ell_list, cEll_1, cEll_2);

						point_list.clear();

						for (Pts_Iterator it = AddPoints2NewPoly_1_dummy.begin(); it != AddPoints2NewPoly_1_dummy.end(); it++){
							AddPoints2NewPoly_1_new.push_back(*it);
						}

						AddPoints2NewPoly_1 = AddPoints2NewPoly_1_new;
						AddPoints2NewPoly_1_dummy.clear();
						AddPoints2NewPoly_1_new.clear();
					}
					break;
					}
				}
				else
				{
					switch ( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).Ell_name_){
					case '1':
						// Case --- 2 --- 1 --- 1 --- 2
					{
						point_list.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_ );

						//AddPoints2NewPoly_1_dummy = CGAL::divide_intersection_segment_2poly<Kernel,NT,Ell_>(point_list, cEll_1, cEll_2);
						AddPoints2NewPoly_1_dummy =
								CGAL::divide_intersection_segment_2poly_ellList<Kernel,NT,Ell_,Output_IPL>(point_list, ell_list, cEll_1, cEll_2);

						point_list.clear();

						for (Pts_Iterator it = AddPoints2NewPoly_1_dummy.begin(); it != AddPoints2NewPoly_1_dummy.end(); it++){
							AddPoints2NewPoly_1_new.push_back(*it);
						}

						AddPoints2NewPoly_1 = AddPoints2NewPoly_1_new;
						AddPoints2NewPoly_1_dummy.clear();
						AddPoints2NewPoly_1_new.clear();
					}
					break;

					case '2':
						// Case --- 2 --- 1 --- 2 --- 1
					{

						point_list.push_back( (*(boost::next(list_IP_Ell_output.begin()))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(list_IP_Ell_output.begin())))).pts_ );
						point_list.push_back( (*(boost::next(boost::next(boost::next(list_IP_Ell_output.begin()))))).pts_ );

						//AddPoints2NewPoly_1_dummy = CGAL::divide_intersection_segment_2poly<Kernel,NT,Ell_>(point_list, cEll_1, cEll_2);
						AddPoints2NewPoly_1_dummy =
								CGAL::divide_intersection_segment_2poly_ellList<Kernel,NT,Ell_,Output_IPL>(point_list, ell_list, cEll_1, cEll_2);

						point_list.clear();

						for (Pts_Iterator it = AddPoints2NewPoly_1_dummy.begin(); it != AddPoints2NewPoly_1_dummy.end(); it++){
							AddPoints2NewPoly_1_new.push_back(*it);
						}

						AddPoints2NewPoly_1 = AddPoints2NewPoly_1_new;
						AddPoints2NewPoly_1_dummy.clear();
						AddPoints2NewPoly_1_new.clear();
					}
					break;
					}
				}

				if (VERY_VERBOSE == 1){
					std::cout << "=======================================" << std::endl;
					std::cout << "AddPoints2NewPoly_1 after updating" << std::endl;
					std::cout << "=======================================" << std::endl;
					for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
						std::cout << "AddPoints2NewPoly_1 = " << *iter << std::endl;
					}
				}

				// First of all, we must check conformity (see check_conformity_newpoly.hh for more details)
				CGAL::check_conformity_newpoly<Kernel, SurfaceMesh_, Point_3>(poly_1_child_1, AddPoints2NewPoly_1);
				reverse_lst_1 = !CGAL::bool_conformity_newpoly<Kernel, SurfaceMesh_, Point_3>(poly_1_child_1, AddPoints2NewPoly_1);

				// Update AddPoints2NewPoly_1 with some new intermediate points to ensure homogeneity of triangle edge length
				std::list<Point_3> AddPoints2NewPoly_1_updated;

				AddPoints2NewPoly_1_updated.push_back(*AddPoints2NewPoly_1.begin());
				unsigned int ip_ =0;
				for( ; ip_ < AddPoints2NewPoly_1.size()-1; ++ip_) {
					Pts_Iterator it_p1, it_p2;
					it_p1 = AddPoints2NewPoly_1.begin();
					std::advance(it_p1, ip_);
					it_p2 = AddPoints2NewPoly_1.begin();
					std::advance(it_p2, ip_+1);

					//int numAddedPoints = round(CGAL::abs(*it_p2 - *it_p1)/ target_edge_length);
					double distance_ =  CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(*it_p2, *it_p1)));
					int numAddedPoints(0);
					//numAddedPoints = ceil(distance_ / target_edge_length)-1;
					numAddedPoints = round(distance_ / target_edge_length);

					if (VERY_VERBOSE == 1){
						std::cout << "Number of added points (poly_1) = " << numAddedPoints << std::endl;
					}

					for (int i = 0; i!= (numAddedPoints+1); i++)
						AddPoints2NewPoly_1_updated.push_back(*it_p1 + (i+1) * (*it_p2 - *it_p1)/(numAddedPoints+1));

				}
				AddPoints2NewPoly_1 = AddPoints2NewPoly_1_updated;
				AddPoints2NewPoly_1_orig = AddPoints2NewPoly_1;

				if (VERY_VERBOSE == 1){
					for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
						std::cout << "Point to add 1 : " << *iter << std::endl;
					}
				}

				// Remove too close points _ NTD 14/09/2016 _ Deactivate
				//remove_close_points <SurfaceMesh_, Point_3, NT_MP> (poly_1_child_1, AddPoints2NewPoly_1, list_seg_interst1);

				if (VERY_VERBOSE == 1){
					for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
						std::cout << "Point to add 1 new : " << *iter << std::endl;
					}

				}

				if (VERY_VERY_VERBOSE == 1){
					std::cout << "(1) poly_1_child_1.number_of_vertices() = " << poly_1_child_1.number_of_vertices() << std::endl;
				}

				// Remove too close points _ NTD 10/10/2016 : when distance between two points is inferior than target_edge_length/5
				if (poly_1_child_1.number_of_vertices() > 3){
					poly_1_child_1 = remove_close_points_from_poly <SurfaceMesh_, Point_3, NT_MP>
										(poly_1_child_1, AddPoints2NewPoly_1, tol_length);
				}

				// Part 2: Add points from the intersection of the polygon_1 to the new "child" poly_1_child_1
				for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
					poly_1_child_1.add_vertex(*iter);

					if (VERY_VERBOSE == 1){
						std::cout << "Point to add: " << *iter << std::endl;
					}
				}

				if (VERY_VERY_VERBOSE == 1){
					std::cout << "(2) poly_1_child_1.number_of_vertices() = " << poly_1_child_1.number_of_vertices() << std::endl;
				}

				poly_1_child_1.add_face(poly_1_child_1.vertices());

				if (PRINT_VTU == 1){
					fname = "out_sm_child_1_1.off";
					CGAL::output_off_vtu <SurfaceMesh_>(fname, poly_1_child_1, false, true);
				}
			}

			// Child polygon poly_1_child_2
			{
				//check_conformity_poly_list(polygon1, AddPoints2NewPoly_1, list_seg_interst1);
				//remove_close_points <SurfaceMesh_, Point_3, NT_MP> (polygon1, AddPoints2NewPoly_1, list_seg_interst1);

				// Add points from the old "father" polygon_1 to the new "child" poly_1_child_2
				unsigned int begin = Min_List<int>(list_seg_interst1);
				unsigned int end = Max_List<int>(list_seg_interst1);


				for(int it=begin;  it>-1; --it) {

					vertex_descriptor vb_seg1(it+1);
					vertex_descriptor ve_seg1(it);

					poly_1_child_2.add_vertex(polygon1.point(ve_seg1));
				}

				begin = polygon1.num_vertices()-1;

				for( ; begin >end; --begin) {
					vertex_descriptor vb_seg1(begin);
					vertex_descriptor ve_seg1(begin-1);

					poly_1_child_2.add_vertex(polygon1.point(vb_seg1));
				}

				// Add points from the intersection of the polygon_1 to the new "child" poly_1_child_2
				AddPoints2NewPoly_1 = AddPoints2NewPoly_1_orig;

				// First of all, we must check conformity (see check_conformity_newpoly.hh for more details)
				CGAL::check_conformity_newpoly<Kernel, SurfaceMesh_, Point_3>(poly_1_child_2, AddPoints2NewPoly_1);

				if (VERY_VERY_VERBOSE == 1){
					std::cout << "(1) poly_1_child_2.number_of_vertices() = " << poly_1_child_2.number_of_vertices() << std::endl;
				}

				// Remove too close points _ NTD 10/10/2016 : when distance between two points is inferior than tol_length
				if (poly_1_child_2.number_of_vertices() > 3){
					poly_1_child_2 = remove_close_points_from_poly <SurfaceMesh_, Point_3, NT_MP>
									(poly_1_child_2, AddPoints2NewPoly_1, tol_length);
				}

				for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
					poly_1_child_2.add_vertex(*iter);
					if (VERY_VERY_VERBOSE == 1){
						std::cout << "Point to add: " << *iter << std::endl;
					}
				}

				if (VERY_VERY_VERBOSE == 1){
					for (vertex_iterator it = poly_1_child_2.vertices_begin(); it != poly_1_child_2.vertices_end(); ++it){
						std::cout << "Point of poly_1_child_2 : " << poly_1_child_2.point(*it) << std::endl;
					}
				}

				if (VERY_VERY_VERBOSE == 1){
					std::cout << "(2) poly_1_child_2.number_of_vertices() = " << poly_1_child_2.number_of_vertices() << std::endl;
				}

				// Add face
				poly_1_child_2.add_face(poly_1_child_2.vertices());

				if (PRINT_VTU == 1){
					fname = "out_sm_child_1_2.off";
					CGAL::output_off_vtu <SurfaceMesh_>(fname, poly_1_child_2, false, true);
				}
			}

			// NTD 12/10/2016
			AddPoints2IntersectionLine = AddPoints2NewPoly_1;
			// Add end points from AddPoints2NewPoly_2
			for (Pts_Iterator iter=AddPoints2NewPoly_2.begin(); iter!=AddPoints2NewPoly_2.end(); iter++){
				if (VERY_VERY_VERBOSE == 1){
					std::cout << "AddPoints2NewPoly_2[] = " << *iter << std::endl;
				}
				AddPoints2IntersectionLine.push_back(*iter);
			}
			//intersectionline = AddPoints2NewPoly_1;
			intersectionline = AddPoints2IntersectionLine;


			if (VERY_VERY_VERBOSE == 1){
				std::cout << "AddPoints2NewPoly_1.size() = " << AddPoints2NewPoly_1.size() << std::endl;
				std::cout << "AddPoints2IntersectionLine.size() = " << AddPoints2IntersectionLine.size() << std::endl;
			}
		}
	}
	output_.poly_1_child_1_ = poly_1_child_1;
	output_.poly_1_child_2_ = poly_1_child_2;
	output_.no_polygons_interst_ = no_polygons_interst;
	output_.intersectionline_ =  intersectionline;

	if (!no_polygons_interst){
		if (VERY_VERBOSE == 1){
			std::cout <<"Cutting off the polygon done."<< std::endl;
		}
	}
	return output_;
}

} // end of namespace



#endif /* INCLUDE_CUTTING_OFF_FIRST_POLYGON_RMPTS_HH_ */
