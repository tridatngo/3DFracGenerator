/*
 * cutting_off_polygons_bbox.hh
 *
 *  Created on: 17 août 2016
 *      Author: ngotr
 */

/*
 *  Cutting off the polygons along the interface line with the bounding box
 */

#ifndef CUTTING_OFF_POLYGONS_BBOX_HH_
#define CUTTING_OFF_POLYGONS_BBOX_HH_

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
#include <include/Min_Max_List.hh>

// Transpose a list
#include <include/transpose_a_list.hh>

#include <include/remove_close_points.hh>

//Intersection between two planes on which two polygons lie
#include "list_seg_index_unique.hh"
//#include "plane_plane_intersection_convert2MP.hh"
#include "plane_plane_intersection.hh"
//#include "Triangulation_MP_Float.hh"

#include "Output_off_vtu.hh"


// check_conformity_poly_list
//#include "check_conformity_poly_list.hh"
#include "check_conformity_newpoly.hh"

#ifndef VERY_VERBOSE
#define VERY_VERBOSE 0
#endif

#ifndef PRINT_VTU
#define PRINT_VTU 0
#endif

namespace CGAL{

template <typename Kernel, typename Output_Child_Polys, typename IP_Ell, typename Output_IP>
Output_Child_Polys cutting_off_polygons_bbox(CGAL::Surface_mesh< CGAL::Point_3<Kernel> >  &polygon1,
		CGAL::Surface_mesh< CGAL::Point_3<Kernel> >  &polygon2,
		double &target_edge_length){

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

	int num_intersec_1(0), num_intersec_2(0);
	std::list<Point_3> intersection_points_1, intersection_points_2;
	std::list<int> list_seg_interst1, list_seg_interst2;
	bool no_polygons_interst(false);
	Output_IP output_ip;

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (polygon1, polygon2);

	num_intersec_1 = output_ip.num_intersec_1_;
	num_intersec_2 = output_ip.num_intersec_2_;
	intersection_points_1 = output_ip.intersection_points_1_;
	intersection_points_2 = output_ip.intersection_points_2_;
	list_seg_interst1 = output_ip.list_seg_interst1_;
	list_seg_interst2 = output_ip.list_seg_interst2_;

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
		std::cout << "Two polygons do not intersect: num_intersec_1 + num_intersec_2 < 4." << std::endl;
		std::cout << " " << std::endl;
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

		//list_IP_Ell_output = CGAL::Order_IntersectionPoins3D<Point_3, Vector_3, IP_Ell, IPE_list>(list_IP_Ell);
		list_IP_Ell_output = CGAL::Order_IntersectionPoins3D<Kernel, IP_Ell>(list_IP_Ell);

		std::list<Point_3> AddPoints2NewPoly_1, AddPoints2NewPoly_2;

		if ((*list_IP_Ell_output.begin()).Ell_name_ == (*(boost::next(list_IP_Ell_output.begin()))).Ell_name_){
			std::cout << "Two polygons do not intersect: two intersection segments have no common point." << std::endl;
			std::cout << " " << std::endl;
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

			if (VERY_VERBOSE ==1){
				//std::cout<<  CGAL::min( *list_seg_interst1.begin(), *(boost::next(list_seg_interst1.begin()))) <<std::endl;

				//std::cout<< *list_seg_interst1.begin()<<std::endl;

				std::cout<< "Min_List<int>(list_seg_interst1) = " << Min_List<int>(list_seg_interst1) <<std::endl;
				std::cout<< "Max_List<int>(list_seg_interst1) = " << Max_List<int>(list_seg_interst1) <<std::endl;

				std::cout<< "Max_List<int>(list_seg_interst2) = " << Min_List<int>(list_seg_interst2) <<std::endl;
				std::cout<< "Max_List<int>(list_seg_interst2) = " << Max_List<int>(list_seg_interst2) <<std::endl;
			}

			// Child polygon poly_1_child_1
			{
				//check_conformity_poly_list<SurfaceMesh_, Point_3>(polygon1, AddPoints2NewPoly_1, list_seg_interst1);

				// Add points from the old "father" polygon_1 to the new "child" poly_1_child_1
				unsigned int begin = Min_List<int>(list_seg_interst1);
				unsigned int end = Max_List<int>(list_seg_interst1);

				for( ; begin < end; ++begin) {
					vertex_descriptor vb_seg1(begin);
					vertex_descriptor ve_seg1(begin+1);

					poly_1_child_1.add_vertex(polygon1.point(ve_seg1));
				}


				// Add points from the intersection of the polygon_1 to the new "child" poly_1_child_1

				// First of all, we must check conformity (see check_conformity_newpoly.hh for more details)
				CGAL::check_conformity_newpoly<Kernel, SurfaceMesh_, Point_3>(poly_1_child_1, AddPoints2NewPoly_1);

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

					double distance_ =  std::sqrt( CGAL::to_double(CGAL::squared_distance(*it_p2, *it_p1)));
					int numAddedPoints = ceil(distance_ / target_edge_length)-1;

					for (int i = 0; i!= (numAddedPoints+1); i++)
						AddPoints2NewPoly_1_updated.push_back(*it_p1 + (i+1) * (*it_p2 - *it_p1)/(numAddedPoints+1));

				}
				AddPoints2NewPoly_1 = AddPoints2NewPoly_1_updated;

				// Remove too close points _ NTD 14/09/2016 _ Deactivate
				//remove_close_points <SurfaceMesh_, Point_3, NT_MP> (poly_1_child_1, AddPoints2NewPoly_1, list_seg_interst1);

				for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
					poly_1_child_1.add_vertex(*iter);
				}

				if (VERY_VERBOSE ==1){
					for (vertex_iterator iter=poly_1_child_1.vertices_begin(); iter!=poly_1_child_1.vertices_end(); iter++){
						std::cout <<  "Vertex index of poly_1_child_1 : " << *iter << std::endl;
						std::cout <<  "Vertex point of poly_1_child_1 : " << poly_1_child_1.point(*iter) << std::endl;
					}
				}

				poly_1_child_1.add_face(poly_1_child_1.vertices());
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

				// First of all, we must check conformity (see check_conformity_newpoly.hh for more details)
				CGAL::check_conformity_newpoly<Kernel, SurfaceMesh_, Point_3>(poly_1_child_2, AddPoints2NewPoly_1);

				// Remove too close points _ NTD 14/09/2016 _ Deactivate
				//remove_close_points <SurfaceMesh_, Point_3, NT_MP> (poly_1_child_2, AddPoints2NewPoly_1, list_seg_interst1);

				for (Pts_Iterator iter=AddPoints2NewPoly_1.begin(); iter!=AddPoints2NewPoly_1.end(); iter++){
					poly_1_child_2.add_vertex(*iter);
				}

				if (VERY_VERBOSE ==1){
					for (vertex_iterator iter=poly_1_child_2.vertices_begin(); iter!=poly_1_child_2.vertices_end(); iter++){
						std::cout <<  "Vertex index of poly_1_child_2 : " << *iter << std::endl;
						std::cout <<  "Vertex point of poly_1_child_2 : " << poly_1_child_2.point(*iter) << std::endl;
					}
				}

				// Add face
				poly_1_child_2.add_face(poly_1_child_2.vertices());

			}
		}
	}

	output_.poly_1_child_1_ = poly_1_child_1;
	output_.poly_1_child_2_ = poly_1_child_2;
	output_.no_polygons_interst_ = no_polygons_interst;

	return output_;
}

}  // end of namespace CGAL

#endif /* CUTTING_OFF_POLYGONS_BBOX_HH_ */
