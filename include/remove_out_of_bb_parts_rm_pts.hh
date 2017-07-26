/*
 * remove_out_of_bb_parts_rm_pts.hh
 *
 *  Created on: Oct 11, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_REMOVE_OUT_OF_BB_PARTS_RM_PTS_HH_
#define INCLUDE_REMOVE_OUT_OF_BB_PARTS_RM_PTS_HH_



//#include "plane_plane_intersection_convert2MP.hh"
#include "plane_plane_intersection.hh"
#include "cutting_off_polygons_bbox_rmpts.hh"
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

#define VERY_VERBO 0

namespace CGAL{

template <typename Kernel, typename Output_Inside_BB, typename Output_Child_Polys, typename IP_Ell, typename Output_IP, typename DFN_BB_, typename NT>
//CGAL::Surface_mesh< CGAL::Point_3<Kernel> > remove_out_of_bb_parts_rm_pts
Output_Inside_BB remove_out_of_bb_parts_rm_pts
		(DFN_BB_ &DFN_BB, CGAL::Surface_mesh< CGAL::Point_3<Kernel> > &poly,
				double &target_edge_length, NT &tol_length){

	CGAL::Point_3<Kernel> p_min, p_max;
	Output_Inside_BB output_inBB;

	typedef typename std::list<IP_Ell> IPE_list;

	typedef typename CGAL::Point_3<Kernel>                    	Point_3;
	typedef typename CGAL::Vector_3<Kernel>                    	Vector_3;
	typedef typename CGAL::Line_3<Kernel>                      	Line_3;
	typedef typename CGAL::Segment_3<Kernel>                   	Segment_3;
	typedef typename CGAL::Plane_3<Kernel>                   	Plane_3;
	typedef typename CGAL::Surface_mesh<Point_3>			SurfaceMesh_;

	typedef typename SurfaceMesh_::Vertex_iterator 	vertex_iterator;
	typedef typename SurfaceMesh_::Vertex_index 	vertex_descriptor;
	typedef typename std::list<Point_3> list_Point_3;
	typedef typename list_Point_3::iterator Pts_iterator;

	vertex_descriptor dfn_bb_vert_1(0), dfn_bb_vert_2(2);
	p_min =  DFN_BB.X_L.point(dfn_bb_vert_1);
	p_max =  DFN_BB.Z_U.point(dfn_bb_vert_2);

	/*
	std::cout << "p_min = " << p_min << std::endl;
	std::cout << "p_max = " << p_max << std::endl;
	 */

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

	bool dfn_bound_interst(false);

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

	//std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
	//std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of X_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly, DFN_BB.X_L, target_edge_length, tol_length);

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

	//std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
	//std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of X_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.X_U, target_edge_length, tol_length);

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

	//std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
	//std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Y_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Y_L, target_edge_length, tol_length);

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

	//std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
	//std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Y_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Y_U, target_edge_length, tol_length);

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

	if (VERY_VERBO == 1){
		std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
		std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;
	}

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Z_L

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4 ){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Z_L, target_edge_length, tol_length);

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

	//std::cout << "output_ip.num_intersec_1_ = " << output_ip.num_intersec_1_ << std::endl;
	//std::cout << "output_ip.num_intersec_2_ = " << output_ip.num_intersec_2_ << std::endl;

	num_intersec_1 = output_ip.num_intersec_1_;

	if ( output_ip.num_intersec_2_ > 0){
		dfn_bound_interst = true;
	}

	if (num_intersec_1 > 0){
		intersection_points_1 = output_ip.intersection_points_1_;
		list_seg_interst1 = output_ip.list_seg_interst1_;

		// Choose a reference point lying on the intersection line -> p_ref
		p_ref = *intersection_points_1.begin();

		// Check if poly_1_child_1_ and p_max are on the same side of Z_U

		if (output_ip.num_intersec_1_ + output_ip.num_intersec_2_ == 4){
			child_polys = CGAL::cutting_off_polygons_bbox_rmpts
					<Kernel, Output_Child_Polys, IP_Ell, Output_IP>
			(poly_new, DFN_BB.Z_U, target_edge_length, tol_length);

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

	bool inside_bb(true);
	if (dfn_bound_interst == false){
		for (vertex_iterator it = poly_new.vertices_begin(); it != poly_new.vertices_end(); ++it){
			//if ((*it). < p_min[0] || (*it)[0] > p_max[0] || (*it)[1] < p_min[1] || (*it)[1] > p_max[1] ||(*it)[2] < p_min[2] || (*it)[2] > p_max[2]){
			if (poly_new.point(*it)[0]  < p_min[0] || poly_new.point(*it)[0]  > p_max[0] ||
				poly_new.point(*it)[1]  < p_min[1] || poly_new.point(*it)[1]  > p_max[1] ||
				poly_new.point(*it)[2]  < p_min[2] || poly_new.point(*it)[2]  > p_max[2]){

				/*
				std::cout << std::endl << "poly_new.point(*it) = "
						<< poly_new.point(*it)[0] << ", "
						<< poly_new.point(*it)[1] << ", "
						<< poly_new.point(*it)[2] << std::endl<< std::endl;
				 */

				inside_bb = false;
				break;
			}
		}
	}

	//std::cout << "inside_bb = " << inside_bb	<< std::endl;
	//std::cout << "dfn_bound_interst = " << dfn_bound_interst << std::endl;

	output_inBB.poly_inside_bb 	= poly_new;
	output_inBB.inside_bb 		= inside_bb;

	return output_inBB;
}
}

#endif /* INCLUDE_REMOVE_OUT_OF_BB_PARTS_RM_PTS_HH_ */
