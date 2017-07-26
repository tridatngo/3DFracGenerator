/*
 * remove_out_of_bb.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 *
 *      Remove out-of-DFN_BB parts of polygons (ellipses)
 */

//#include "plane_plane_intersection_convert2MP.hh"
#include "plane_plane_intersection.hh"
#include "cutting_off_polygons_bbox.hh"
#include <list>

// Add some boost files
#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/include/deref.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/include/begin.hpp>

// Return the lexicographically smallest and largest coordinates of bb
#include "include/find_extrema_bb.hh"


#ifndef REMOVE_OUT_OF_BB_HH_
#define REMOVE_OUT_OF_BB_HH_

#define VERY_VERBO 0

namespace CGAL{

template <typename Kernel, typename Output_Child_Polys, typename IP_Ell, typename Output_IP, typename DFN_BB_>
CGAL::Surface_mesh< CGAL::Point_3<Kernel> > remove_out_of_bb_parts (DFN_BB_ &DFN_BB, CGAL::Surface_mesh< CGAL::Point_3<Kernel> > &poly, double &target_edge_length){

	typedef typename std::list<IP_Ell> IPE_list;

	typedef typename CGAL::Point_3<Kernel>                    	Point_3;
	typedef typename CGAL::Vector_3<Kernel>                    	Vector_3;
	typedef typename CGAL::Line_3<Kernel>                      	Line_3;
	typedef typename CGAL::Segment_3<Kernel>                   	Segment_3;
	typedef typename CGAL::Plane_3<Kernel>                   	Plane_3;
	typedef typename CGAL::Surface_mesh<Point_3>			SurfaceMesh_;

	typedef typename SurfaceMesh_::Vertex_iterator vertex_iterator;
	typedef typename std::list<Point_3> list_Point_3;
	typedef typename list_Point_3::iterator Pts_iterator;

	Point_3 p_ref;

	Output_IP output_ip;
	Output_Child_Polys child_polys;
	SurfaceMesh_ poly_new;
	poly_new = poly;
	Vector_3 vec_1, vec_2;

	int num_intersec_1(0);
	std::list<Point_3> intersection_points_1;
	std::list<int> list_seg_interst1;
	bool no_polygons_interst;
	Pts_iterator p_begin, p_end;

	bool dfn_bound_interst;

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (X_L)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and X_L ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly, DFN_BB.X_L);

	num_intersec_1 = output_ip.num_intersec_1_;

	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of X_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly, DFN_BB.X_L, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(1, 0, 0); // Normal vector of X_L

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (X_U)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and X_U ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly_new, DFN_BB.X_U);

	num_intersec_1 = output_ip.num_intersec_1_;

	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of X_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.X_U, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(-1, 0, 0); // Normal vector of X_U

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (Y_L)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and Y_L ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly_new, DFN_BB.Y_L);

	num_intersec_1 = output_ip.num_intersec_1_;

	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Y_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Y_L, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(0, 1, 0); // Normal vector of Y_L

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (Y_U)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and Y_U ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly_new, DFN_BB.Y_U);

	num_intersec_1 = output_ip.num_intersec_1_;
	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Y_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Y_U, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(0, -1, 0); // Normal vector of Y_U

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (Z_L)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and Z_L ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly_new, DFN_BB.Z_L);

	num_intersec_1 = output_ip.num_intersec_1_;

	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Z_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4 ){
			child_polys = CGAL::cutting_off_polygons_bbox<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Z_L, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(0, 0, 1); // Normal vector of Z_L

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	/*
	 *  #######################################################################################
	 * Checking of intersection between the polygon and the lower face of the DFN_BB (Z_U)
	 * #######################################################################################
	 */

	if (VERY_VERBO == 1){
		std::cout << "Checking of intersection between the polygon and Z_U ..." << std::endl;
	}

	output_ip = CGAL::intersections_of_two_polys<Kernel,Output_IP> (poly_new, DFN_BB.Z_U);

	num_intersec_1 = output_ip.num_intersec_1_;

	if (num_intersec_1 > 0){
		dfn_bound_interst = true;
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Z_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Z_U, target_edge_length);

			if (VERY_VERBO == 1){
				std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;
			}

			if (!child_polys.no_polygons_interst_){
				poly_new = child_polys.poly_1_child_1_;
				vec_1 = Vector_3(0, 0, -1); // Normal vector of Z_U

				for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
					if ( !CGAL::collinear(poly_new.point(*it), *intersection_points_1.begin(), *(boost::next(intersection_points_1.begin())) ) ){

						//std::cout << "Point on polys :"<< poly_new.point(*it) << std::endl;
						vec_2 = Vector_3(p_ref, poly_new.point(*it) );
						break;
					}
				}

				if (vec_1 * vec_2 < 0){
					poly_new = child_polys.poly_1_child_2_;
				}
			}
		}
	}

	if (VERY_VERBO == 1){
		std::cout << "=======================   Done   =======================" << std::endl;
		std::cout <<"\n";
	}

	std::cout << "dfn_bound_interst = " << dfn_bound_interst << std::endl;

	return poly_new;
}
}

#endif /* REMOVE_OUT_OF_BB_HH_ */
