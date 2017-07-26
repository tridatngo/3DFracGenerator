/*
 * order_intersectionpoints_3d.hh
 *
 *  Created on: 1 août 2016
 *  Modified 19 août 2016: V_19_08_2016
 *      Author: ngotr
 */

/*
 * The functions in this header file are in order to sort the intersection points
 */

//#include "Exact_predicates_inexact_constructions_kernel_ntd.h"

#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Kernel/global_functions.h>


#include <algorithm>
#include <iterator>     // std::next & std::prev
#include <list>
#include <vector>

// Intersection
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>

// To check all boolean list elements are true
#include "All_bool_true.hh"

#include "points_are_close.hh"
#include "transpose_a_list.hh"
#include "list_seg_index_unique.hh"

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel_MP_Float Kernel;

#ifndef USE_EPECK
#define USE_EPECK 0
#endif

#if (USE_EPECK)
typedef CGAL::Exact_predicates_exact_constructions_kernel 				Kernel;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel_MP_Float 	Kernel;
#endif

typedef Kernel::Point_3                    Point_3;
typedef Kernel::Vector_3                   Vector_3;
typedef Kernel::Line_3                     Line_3;
typedef Kernel::Segment_3                  Segment_3;

//typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT_MP;
typedef CGAL::Inverse_index<std::list<Point_3>::iterator> Inverse_Index;

typedef std::list<Point_3>::iterator Pts_Iterator;
typedef std::list<int>::iterator Int_Iterator;
typedef std::list<char>::iterator Char_Iterator;
typedef std::list<Vector_3>::iterator Vector_Iterator;


#ifndef CGAL_order_intersectionpoints_3d_hh
#define CGAL_order_intersectionpoints_3d_hh

/* [IMP] Prerequisite : All points must be collinear. */
#define VERY_VERBOSE 0

namespace CGAL {

template <typename Point_3>
void Order_IntersectionPoins3D_3pts (const Point_3 &pt_1, const Point_3 &pt_2, const Point_3 &pt_3){

	if (CGAL::collinear_are_strictly_ordered_along_line(pt_1, pt_2, pt_3)){
		std::cout << "The point (" << pt_2 << ") lies between two points (" << pt_1 <<") and (" << pt_3 <<")" << std::endl;
	}
	else
		if (CGAL::collinear_are_strictly_ordered_along_line(pt_2, pt_1 , pt_3)){
			std::cout << "The point (" << pt_1 << ") lies between two points (" << pt_2 <<") and (" << pt_3 <<")" << std::endl;
		}
		else{
			std::cout << "The point (" << pt_3 << ") lies between two points (" << pt_1 <<") and (" << pt_2 <<")" << std::endl;
		}
}

template <typename Point_3>
void Order_IntersectionPoins3D_3pts_reorder (const Point_3 &pt_1, const Point_3 &pt_2, const Point_3 &pt_3,
		const Point_3 &pt_1_new, const Point_3 &pt_2_new, const Point_3 &pt_3_new){

	if (CGAL::collinear_are_strictly_ordered_along_line(pt_1, pt_2, pt_3)){
		std::cout << "The point (" << pt_2 << ") lies between two points (" << pt_1 <<") and (" << pt_3 <<")" << std::endl;
		pt_1_new = pt_1; pt_2_new = pt_2; pt_3_new = pt_3;
	}
	else
		if (CGAL::collinear_are_strictly_ordered_along_line(pt_2, pt_1 , pt_3)){
			std::cout << "The point (" << pt_1 << ") lies between two points (" << pt_2 <<") and (" << pt_3 <<")" << std::endl;
			pt_1_new = pt_2; pt_2_new = pt_1; pt_3_new = pt_3;
		}
		else{
			std::cout << "The point (" << pt_3 << ") lies between two points (" << pt_1 <<") and (" << pt_2 <<")" << std::endl;
			pt_1_new = pt_3; pt_2_new = pt_1; pt_3_new = pt_2;
		}
}

//template <class Point_3, class Vector_3>
//template <typename Point_3>
template <typename Point_3, typename Vector_3>
void Order_IntersectionPoins3D_4pts (const Point_3 &pt_1, const Point_3 &pt_2, const Point_3 &pt_3, const Point_3 &pt_4){
	/*
	std::cout << " " << std::endl;
	std::cout << "############     Order of points    ############" << std::endl;
	std::cout << "-----   1   ---   2   ---   3   ---   4   ------" << std::endl;
	std::cout << "################################################" << std::endl;
	std::cout << " " << std::endl;
	 */

	// Search the segment ends
	std::list<Point_3> list_pts_, list_pts_new_, endpoints_, list_pts_reordered;
	std::list<Vector_3> list_vect_;
	std::list< std::list<Vector_3> > all_list_vect_;
	std::list<int> index_reordered;

	list_pts_.push_back(pt_1);
	list_pts_.push_back(pt_2);
	list_pts_.push_back(pt_3);
	list_pts_.push_back(pt_4);

	bool test_endpts(true); // Test endpoints
	int idx = 1;

	for (Pts_Iterator it1=list_pts_.begin(); it1!=list_pts_.end(); ++it1){
		test_endpts = true;

		for (Pts_Iterator it2=list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			list_vect_.push_back(Vector_3(*it1, *it2));
		}

		for (Vector_Iterator it3=list_vect_.begin(); it3!=list_vect_.end(); ++it3){
			for (Vector_Iterator it4=list_vect_.begin(); it4!=list_vect_.end(); ++it4){
				//std::cout << " test_endpts = " << test_endpts << std::endl;
				//std::cout << " (*it3) * (*it4) " << (*it3) * (*it4) << std::endl;
				if ((*it3) * (*it4) < 0){
					test_endpts = false;
					break;
				}
			}
			if(!test_endpts)
				break;
		}

		list_vect_.clear();


		if (test_endpts){
			std::cout << "The point " << idx << " is an endpoint." << std::endl;
			endpoints_.push_back(*it1);
		}

		idx++;
	}

	// Remove one endpoint and order the others
	Pts_Iterator i_ep_begin_ = endpoints_.begin();
	Pts_Iterator i_ep_end_ = i_ep_begin_;
	std::advance(i_ep_end_, 1);

	list_pts_new_ = list_pts_;
	for (Pts_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( (*it1 == *i_ep_begin_) || (*it1 == *i_ep_end_)){
			list_pts_new_.erase(it1);
			break;
		}
	}

	std::list<NT_MP> squ_distance_;

	list_pts_reordered.push_back(*i_ep_begin_);

	for (int i=0;i!=1;i++){
		Pts_Iterator it1=list_pts_new_.begin();
		Pts_Iterator it2 = it1;
		std::advance(it2, 1);
		if (CGAL::squared_distance(*i_ep_begin_, *it1) < CGAL::squared_distance(*i_ep_begin_, *it2) ){
			list_pts_reordered.push_back(*it1);
			list_pts_reordered.push_back(*it2);
		}
		else{
			list_pts_reordered.push_back(*it2);
			list_pts_reordered.push_back(*it1);
		}
	}

	list_pts_reordered.push_back(*i_ep_end_);

	// Print out the initial list of points
	for (Pts_Iterator it2 = list_pts_.begin(); it2!=list_pts_.end(); ++it2){
		std::cout << "Point " << *it2 << std::endl;}

	// Reorder the list of point
	for (Pts_Iterator it1 = list_pts_reordered.begin(); it1!=list_pts_reordered.end(); ++it1){
		int idx_= 0;
		for (Pts_Iterator it2 = list_pts_.begin(); it2!=list_pts_.end(); ++it2){

			if (*it1 == *it2){
				index_reordered.push_back(idx_);}
			idx_++;
		}
	}

	Int_Iterator idx_f_1 = index_reordered.begin();
	Int_Iterator idx_f_2 = index_reordered.begin();
	Int_Iterator idx_f_3 = index_reordered.begin();
	Int_Iterator idx_f_4 = index_reordered.begin();

	std::advance(idx_f_2, 1);
	std::advance(idx_f_3, 2);
	std::advance(idx_f_4, 3);

	if (VERY_VERBOSE == 1){
		std::cout << " " << std::endl;
		std::cout << "############     Order of points    ############" << std::endl;
		std::cout << "-----   "<< *idx_f_1 <<"   ---   "<< *idx_f_2 <<"   ---   "
				<< *idx_f_3 <<"   ---   "<< *idx_f_4 <<"   ------" << std::endl;
		std::cout << "################################################" << std::endl;
		std::cout << " " << std::endl;
	}
}

// Merge too close intersection points - NTD 18/11
template <typename Kernel, typename IP_Ell>
std::list<IP_Ell> Order_IntersectionPoins3D_merge_close_pts (const std::list<IP_Ell> &list_IP_Ell_in, const double crit_length){

	typedef typename std::list<IP_Ell> 					IPE_list;
	typedef typename std::list<IP_Ell>::const_iterator 	IPE_Iterator;
	typedef typename CGAL::Point_3<Kernel>              Point_3;
	typedef typename CGAL::Vector_3<Kernel>             Vector_3;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator 		Vector_Iterator;
	typedef typename std::list<IP_Ell>::const_iterator 	IPE_Iterator;

	IPE_list list_IP_Ell_out;
	IP_Ell IPE_new_1, IPE_new_2, IPE_new_3, IPE_new_4;


	std::list<Point_3> list_pts_, list_pts_new_, endpoints_, list_pts_reordered;
	std::list<Vector_3> list_vect_;
	std::list< std::list<Vector_3> > all_list_vect_;
	std::list<char> list_chars_, list_chars_reordered;

	if (VERY_VERBOSE == 1){
		std::cout << "############     Points from input    ############" << std::endl;
	}

	for (IPE_Iterator it=list_IP_Ell_in.begin(); it!=list_IP_Ell_in.end(); it++){
		list_pts_.push_back((*it).pts_);

		/*
		if (VERY_VERBOSE == 1){
			std::cout << "-----   "<< (*it).pts_ <<"   ---   "<<std::endl;
		}
		 */

		if (VERY_VERBOSE == 1){
			std::cout << "-----   "<< (*it).pts_ <<"   ---   " <<" from polygon " << (*it).Ell_name_ << "   ---   "<<std::endl;
		}

		list_chars_.push_back((*it).Ell_name_);
	}


	bool test_endpts(true); // Test endpoints
	int idx = 1;

	/*
	 * Search the endpoints :
	 * A point is a end-point if and only if all scalar product between two vectors
	 *  created by this point and other one are positive
	 */

	for (Pts_Iterator it1=list_pts_.begin(); it1!=list_pts_.end(); ++it1){
		test_endpts = true;

		for (Pts_Iterator it2=list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			list_vect_.push_back(Vector_3(*it1, *it2));
		}

		for (Vector_Iterator it3=list_vect_.begin(); it3!=list_vect_.end(); ++it3){
			for (Vector_Iterator it4=list_vect_.begin(); it4!=list_vect_.end(); ++it4){
				if ((*it3) * (*it4) < 0){
					test_endpts = false;
					break;
				}
			}
			if(!test_endpts)
				break;
		}

		list_vect_.clear();


		if (test_endpts){
			if (VERY_VERBOSE == 1){
				//std::cout << "The point " << idx << " is an endpoint." << std::endl;
				std::cout << "The point " << *it1 << " is an endpoint." << std::endl;
			}
			endpoints_.push_back(*it1);
		}

		idx++;
	}

	Pts_Iterator ipe_it_pt_1, ipe_it_pt_2, ipe_it_pt_3, ipe_it_pt_4;
	Char_Iterator ipe_it_ch_1, ipe_it_ch_2, ipe_it_ch_3, ipe_it_ch_4;

	//list_points_unique<Point_3, NT_MP>(endpoints_);
	// V_19_08_2016
	endpoints_.sort();
	endpoints_.unique();

	if (VERY_VERBOSE == 1){
		std::cout << "############     List end points    ############" << std::endl;
		for (Pts_Iterator it=endpoints_.begin(); it!=endpoints_.end(); it++){
			std::cout <<  (*it) <<std::endl;
		}
	}

	// Remove the two end-points and re-order the others
	Pts_Iterator i_ep_begin_ = endpoints_.begin();
	Pts_Iterator i_ep_end_ = i_ep_begin_;
	std::advance(i_ep_end_, 1);

	list_pts_new_ = list_pts_;
	for (Pts_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( (*it1 == *i_ep_begin_)){
			list_pts_new_.erase(it1);
			break;
		}
	}

	for (Pts_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( (*it1 == *i_ep_end_)){
			list_pts_new_.erase(it1);
			break;
		}
	}
	if (VERY_VERBOSE == 1){
		std::cout << "############     List new    ############" << std::endl;
		for (Pts_Iterator it=list_pts_new_.begin(); it!=list_pts_new_.end(); it++){
			std::cout <<  (*it) <<std::endl;
		}
	}

	std::list<NT_MP> squ_distance_;

	list_pts_reordered.push_back(*i_ep_begin_);

	for (int i=0;i!=1;i++){
		Pts_Iterator it1=list_pts_new_.begin();
		Pts_Iterator it2 = it1;
		std::advance(it2, 1);
		if (CGAL::squared_distance(*i_ep_begin_, *it1) < CGAL::squared_distance(*i_ep_begin_, *it2) ){
			list_pts_reordered.push_back(*it1);
			list_pts_reordered.push_back(*it2);
		}
		else{
			list_pts_reordered.push_back(*it2);
			list_pts_reordered.push_back(*it1);
		}
	}

	list_pts_reordered.push_back(*i_ep_end_);


	// [1st time] Reorder the list of point

	ipe_it_pt_1 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_1, 0);
	IPE_new_1.pts_ = *ipe_it_pt_1;

	ipe_it_pt_2 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_2, 1);
	IPE_new_2.pts_ = *ipe_it_pt_2;

	ipe_it_pt_3 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_3, 2);
	IPE_new_3.pts_ = *ipe_it_pt_3;

	ipe_it_pt_4 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_4, 3);
	IPE_new_4.pts_ = *ipe_it_pt_4;


	// Be careful  NTD 14/09/2016
	Char_Iterator ipe_it_ch = list_chars_.begin();

	for (Pts_Iterator it1 = list_pts_reordered.begin(); it1!=list_pts_reordered.end(); ++it1){
		int idx_= 0;
		for (Pts_Iterator it2 = list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			if (*it1 == *it2){
				ipe_it_ch = list_chars_.begin();
				std::advance(ipe_it_ch, idx_);
				list_chars_reordered.push_back(*ipe_it_ch);}
			idx_++;
		}
	}

	if (VERY_VERBOSE == 1){
		std::cout << "############     list_chars_reordered    ############" << std::endl;
		std::cout << "-----   ";
		for (Char_Iterator it = list_chars_reordered.begin(); it != list_chars_reordered.end(); it++){
			std::cout<< * it <<"   -----   ";
		}
		std::cout <<"\n";
	}

	Pts_Iterator pts_ib(list_pts_reordered.begin());
	Char_Iterator char_ib(list_chars_reordered.begin());
	if (list_chars_reordered.size() > list_pts_reordered.size()){

		if ( *pts_ib == *(boost::next(pts_ib))){
			if ( *(boost::next(boost::next(pts_ib))) != *(boost::next(boost::next(boost::next(pts_ib)))) ){

				// Case: A    A    B    C  --->  1 --- 2 --- 1 --- 2 --- 2 --- 1
				list_chars_reordered.erase(char_ib);
				char_ib = list_chars_reordered.begin();
				list_chars_reordered.erase(char_ib);
			}
			else
			{
				// Case: A    A    B    B  --->  1 --- 2 --- 1 --- 2 --- 1 --- 2 --- 1 --- 2
				list_chars_reordered.erase(char_ib);
				for (int it =0; it!=3; it++){
					char_ib = list_chars_reordered.begin();
					list_chars_reordered.erase(char_ib);
				}
			}
		}

		else{
			if ( *(boost::next(pts_ib)) == *(boost::next(boost::next(pts_ib))) ){
				if (*char_ib == *(boost::next(char_ib))){
					// Case: A    B    B    C  --->  1 --- 1 --- 2 --- 1 --- 2 --- 2
					std::advance(char_ib,1);
					char_ib = list_chars_reordered.begin();
					std::advance(char_ib,3);
					list_chars_reordered.erase(char_ib);
				}
				else{
					// Case: C    B    B    A  --->  2 --- 1 --- 2 --- 1 --- 2 --- 1
					std::advance(char_ib,2);
					char_ib = list_chars_reordered.begin();
					std::advance(char_ib,2);
					list_chars_reordered.erase(char_ib);

				}
			}

		}
	}

	if (VERY_VERBOSE == 1){
		std::cout << "############     list_chars_reordered    ############" << std::endl;
		std::cout << "-----   ";
		for (Char_Iterator it = list_chars_reordered.begin(); it != list_chars_reordered.end(); it++){
			std::cout<< * it <<"   -----   ";
		}
		std::cout <<"\n";
	}
	// Reorder the list of char

	ipe_it_ch_1 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_1, 0);
	IPE_new_1.Ell_name_ = *ipe_it_ch_1;

	ipe_it_ch_2 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_2, 1);
	IPE_new_2.Ell_name_ = *ipe_it_ch_2;

	ipe_it_ch_3 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_3, 2);
	IPE_new_3.Ell_name_ = *ipe_it_ch_3;

	ipe_it_ch_4 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_4, 3);
	IPE_new_4.Ell_name_ = *ipe_it_ch_4;

	/*=================================================================================== */
	// [2nd time] Reorder the list of point - 03/01/2017
	/*
	ipe_it_pt_1 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_1, 0);
	IPE_new_1.pts_ = *ipe_it_pt_1;

	ipe_it_pt_2 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_2, 1);
	if (CGAL::to_double(CGAL::squared_distance(*ipe_it_pt_1, *ipe_it_pt_2)) < crit_length * crit_length){
		IPE_new_2.pts_ = *ipe_it_pt_1;
	}
	else{
		IPE_new_2.pts_ = *ipe_it_pt_2;
	}
	*/

	/*=================================================================================== */
	// [2nd time] Reorder the list of point - 04/01/2017
	ipe_it_pt_1 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_1, 0);

	ipe_it_pt_2 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_2, 1);

	if (CGAL::to_double(CGAL::squared_distance(*ipe_it_pt_1, *ipe_it_pt_2)) < crit_length * crit_length){
		IPE_new_1.pts_ = *ipe_it_pt_2;
		IPE_new_2.pts_ = *ipe_it_pt_2;
	}
	else{
		IPE_new_1.pts_ = *ipe_it_pt_1;
		IPE_new_2.pts_ = *ipe_it_pt_2;
	}
	/*=================================================================================== */

	ipe_it_pt_3 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_3, 2);
	if (CGAL::to_double(CGAL::squared_distance(*ipe_it_pt_2, *ipe_it_pt_3)) < crit_length * crit_length){
		IPE_new_3.pts_ = *ipe_it_pt_2;
	}
	else{
		IPE_new_3.pts_ = *ipe_it_pt_3;
	}

	ipe_it_pt_4 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_4, 3);
	if (CGAL::to_double(CGAL::squared_distance(*ipe_it_pt_3, *ipe_it_pt_4)) < crit_length * crit_length){
		IPE_new_4.pts_ = *ipe_it_pt_3;
	}
	else{
		IPE_new_4.pts_ = *ipe_it_pt_4;
	}

	list_IP_Ell_out.push_back(IPE_new_1);
	list_IP_Ell_out.push_back(IPE_new_2);
	list_IP_Ell_out.push_back(IPE_new_3);
	list_IP_Ell_out.push_back(IPE_new_4);

	if (VERY_VERBOSE == 1){
		int type_print_out = 1;

		switch (type_print_out){
		case 1:
		{
			std::cout << " " << std::endl;
			std::cout << "crit_length = " << crit_length << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_pt_1 <<"   ---   "<< *ipe_it_pt_2 <<"   ---   "
					<< *ipe_it_pt_3 <<"   ---   "<< *ipe_it_pt_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;
			std::cout << "############     On the ellipse    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_ch_1 <<"   ---   "<< *ipe_it_ch_2 <<"   ---   "
					<< *ipe_it_ch_3 <<"   ---   "<< *ipe_it_ch_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;

		}
		break;

		case 2:
		{
			std::cout << " " << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_pt_1 <<"   ---   "<< *ipe_it_pt_2 <<"   ---   "
					<< *ipe_it_pt_3 <<"   ---   "<< *ipe_it_pt_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;

		}
		break;
		}
	}

	return list_IP_Ell_out;
}

// New version 14/09/2016
template <typename Kernel, typename IP_Ell>
std::list<IP_Ell> Order_IntersectionPoins3D_new (const std::list<IP_Ell> &list_pts_){

	typedef typename std::list<IP_Ell> 					IPE_list;
	typedef typename std::list<IP_Ell>::const_iterator 	IPE_Iterator;
	typedef typename CGAL::Point_3<Kernel>              Point_3;
	typedef typename CGAL::Vector_3<Kernel>             Vector_3;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator 		Vector_Iterator;

	IPE_list listendpoints_, list_pts_new_, list_pts_reordered;
	IP_Ell IPE_new_1, IPE_new_2, IPE_new_3, IPE_new_4;

	std::list<Vector_3> list_vect_;
	//std::list< std::list<Vector_3> > all_list_vect_;
	//std::list<char> list_chars_, list_chars_reordered;

	if (VERY_VERBOSE == 1){
		std::cout << "############     Points from input    ############" << std::endl;
	}

	for (IPE_Iterator it=list_pts_.begin(); it!=list_pts_.end(); it++){
		if (VERY_VERBOSE == 1){
			std::cout << "-----   "<< (*it).pts_ <<"   ---   " <<" from polygon " << (*it).Ell_name_ << "   ---   "<<std::endl;
		}
	}

	bool test_endpts(true); // Test endpoints
	int idx = 1;

	/*
	 * Search the endpoints :
	 * A point is a end-point if and only if all scalar product between two vectors
	 *  created by this point and other one are positive
	 */

	for (IPE_Iterator it1=list_pts_.begin(); it1!=list_pts_.end(); ++it1){
		test_endpts = true;

		for (IPE_Iterator it2=list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			list_vect_.push_back(Vector_3((*it1).pts_, (*it2).pts_));
		}

		for (Vector_Iterator it3=list_vect_.begin(); it3!=list_vect_.end(); ++it3){
			for (Vector_Iterator it4=list_vect_.begin(); it4!=list_vect_.end(); ++it4){
				if ((*it3) * (*it4) < 0){
					test_endpts = false;
					break;
				}
			}
			if(!test_endpts)
				break;
		}

		list_vect_.clear();


		if (test_endpts){
			if (VERY_VERBOSE == 1){
				std::cout << "The point " << (*it1).pts_ << " is an endpoint." << std::endl;
			}
			listendpoints_.push_back(*it1);
		}

		idx++;
	}

	//Pts_Iterator ipe_it_pt_1, ipe_it_pt_2, ipe_it_pt_3, ipe_it_pt_4;
	//Char_Iterator ipe_it_ch_1, ipe_it_ch_2, ipe_it_ch_3, ipe_it_ch_4;

	//list_points_unique<Point_3, NT_MP>(listendpoints_);
	// V_19_08_2016
	//listendpoints_.sort();
	//listendpoints_.unique();

	/*
	 * TODO sort and remove duplicate points
	 */

	if (VERY_VERBOSE == 1){
		std::cout << "############     List end points    ############" << std::endl;
		for (IPE_Iterator it=listendpoints_.begin(); it!=listendpoints_.end(); it++){
			std::cout <<  (*it).pts_ <<std::endl;
		}
	}

	// Remove the two end-points and re-order the others
	IPE_Iterator i_ep_begin_ = listendpoints_.begin();
	IPE_Iterator i_ep_end_ = i_ep_begin_;
	std::advance(i_ep_end_, 1);

	list_pts_new_ = list_pts_;
	for (IPE_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( ( (*it1).pts_ == (*i_ep_begin_).pts_)){
			list_pts_new_.erase(it1);
			break;
		}
	}

	for (IPE_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( ((*it1).pts_ == (*i_ep_end_).pts_)){
			list_pts_new_.erase(it1);
			break;
		}
	}

	if (VERY_VERBOSE == 1){
		std::cout << "############     List new    ############" << std::endl;
		for (IPE_Iterator it=list_pts_new_.begin(); it!=list_pts_new_.end(); it++){
			std::cout <<  (*it).pts_ <<std::endl;
		}
	}

	std::list<NT_MP> squ_distance_;

	list_pts_reordered.push_back(*i_ep_begin_);

	for (int i=0;i!=1;i++){
		IPE_Iterator it1=list_pts_new_.begin();
		IPE_Iterator it2 = it1;
		std::advance(it2, 1);
		if (CGAL::squared_distance((*i_ep_begin_).pts_, (*it1).pts_) < CGAL::squared_distance((*i_ep_begin_).pts_, (*it2).pts_) ){
			list_pts_reordered.push_back(*it1);
			list_pts_reordered.push_back(*it2);
		}
		else{
			list_pts_reordered.push_back(*it2);
			list_pts_reordered.push_back(*it1);
		}
	}

	list_pts_reordered.push_back(*i_ep_end_);

	// Reorder the list of point
	IPE_Iterator ipe_it_pt_1, ipe_it_pt_2, ipe_it_pt_3, ipe_it_pt_4;

	ipe_it_pt_1 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_1, 0);

	ipe_it_pt_2 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_2, 1);

	ipe_it_pt_3 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_3, 2);

	ipe_it_pt_4 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_4, 3);

	if (VERY_VERBOSE == 1){
		int type_print_out = 1;

		switch (type_print_out){
		case 1:
		{
			std::cout << " " << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< (*ipe_it_pt_1).pts_ <<"   ---   "<< (*ipe_it_pt_2).pts_ <<"   ---   "
					<< (*ipe_it_pt_3).pts_ <<"   ---   "<< (*ipe_it_pt_4).pts_ <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;
			std::cout << "############     On the ellipse    ############" << std::endl;
			std::cout << "-----   "<< (*ipe_it_pt_1).Ell_name_ <<"   ---   "<< (*ipe_it_pt_2).Ell_name_ <<"   ---   "
					<< (*ipe_it_pt_3).Ell_name_ <<"   ---   "<< (*ipe_it_pt_4).Ell_name_ <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;

		}
		break;

		case 2:
		{
			std::cout << " " << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< (*ipe_it_pt_1).pts_ <<"   ---   "<< (*ipe_it_pt_2).pts_ <<"   ---   "
					<< (*ipe_it_pt_3).pts_ <<"   ---   "<< (*ipe_it_pt_4).pts_ <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;
		}
		break;
		}
	}

	return list_pts_reordered;
}


template <typename Kernel, typename IP_Ell>
std::list<IP_Ell> Order_IntersectionPoins3D (const std::list<IP_Ell> &list_IP_Ell_in){

	typedef typename std::list<IP_Ell> 					IPE_list;
	typedef typename std::list<IP_Ell>::const_iterator 	IPE_Iterator;
	typedef typename CGAL::Point_3<Kernel>              Point_3;
	typedef typename CGAL::Vector_3<Kernel>             Vector_3;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename std::list<Vector_3>::iterator 		Vector_Iterator;
	typedef typename std::list<IP_Ell>::const_iterator 	IPE_Iterator;

	IPE_list list_IP_Ell_out;
	IP_Ell IPE_new_1, IPE_new_2, IPE_new_3, IPE_new_4;


	std::list<Point_3> list_pts_, list_pts_new_, endpoints_, list_pts_reordered;
	std::list<Vector_3> list_vect_;
	std::list< std::list<Vector_3> > all_list_vect_;
	std::list<char> list_chars_, list_chars_reordered;

	if (VERY_VERBOSE == 1){
		std::cout << "############     Points from input    ############" << std::endl;
	}

	for (IPE_Iterator it=list_IP_Ell_in.begin(); it!=list_IP_Ell_in.end(); it++){
		list_pts_.push_back((*it).pts_);

		/*
		if (VERY_VERBOSE == 1){
			std::cout << "-----   "<< (*it).pts_ <<"   ---   "<<std::endl;
		}
		 */

		if (VERY_VERBOSE == 1){
			std::cout << "-----   "<< (*it).pts_ <<"   ---   " <<" from polygon " << (*it).Ell_name_ << "   ---   "<<std::endl;
		}

		list_chars_.push_back((*it).Ell_name_);
	}


	bool test_endpts(true); // Test endpoints
	int idx = 1;

	/*
	 * Search the endpoints :
	 * A point is a end-point if and only if all scalar product between two vectors
	 *  created by this point and other one are positive
	 */

	for (Pts_Iterator it1=list_pts_.begin(); it1!=list_pts_.end(); ++it1){
		test_endpts = true;

		for (Pts_Iterator it2=list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			list_vect_.push_back(Vector_3(*it1, *it2));
		}

		for (Vector_Iterator it3=list_vect_.begin(); it3!=list_vect_.end(); ++it3){
			for (Vector_Iterator it4=list_vect_.begin(); it4!=list_vect_.end(); ++it4){
				if ((*it3) * (*it4) < 0){
					test_endpts = false;
					break;
				}
			}
			if(!test_endpts)
				break;
		}

		list_vect_.clear();


		if (test_endpts){
			if (VERY_VERBOSE == 1){
				//std::cout << "The point " << idx << " is an endpoint." << std::endl;
				std::cout << "The point " << *it1 << " is an endpoint." << std::endl;
			}
			endpoints_.push_back(*it1);
		}

		idx++;
	}

	Pts_Iterator ipe_it_pt_1, ipe_it_pt_2, ipe_it_pt_3, ipe_it_pt_4;
	Char_Iterator ipe_it_ch_1, ipe_it_ch_2, ipe_it_ch_3, ipe_it_ch_4;

	//list_points_unique<Point_3, NT_MP>(endpoints_);
	// V_19_08_2016
	endpoints_.sort();
	endpoints_.unique();

	if (VERY_VERBOSE == 1){
		std::cout << "############     List end points    ############" << std::endl;
		for (Pts_Iterator it=endpoints_.begin(); it!=endpoints_.end(); it++){
			std::cout <<  (*it) <<std::endl;
		}
	}

	// Remove the two end-points and re-order the others
	Pts_Iterator i_ep_begin_ = endpoints_.begin();
	Pts_Iterator i_ep_end_ = i_ep_begin_;
	std::advance(i_ep_end_, 1);

	list_pts_new_ = list_pts_;
	for (Pts_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( (*it1 == *i_ep_begin_)){
			list_pts_new_.erase(it1);
			break;
		}
	}

	for (Pts_Iterator it1=list_pts_new_.begin(); it1!=list_pts_new_.end(); ++it1){
		if ( (*it1 == *i_ep_end_)){
			list_pts_new_.erase(it1);
			break;
		}
	}
	if (VERY_VERBOSE == 1){
		std::cout << "############     List new    ############" << std::endl;
		for (Pts_Iterator it=list_pts_new_.begin(); it!=list_pts_new_.end(); it++){
			std::cout <<  (*it) <<std::endl;
		}
	}

	std::list<NT_MP> squ_distance_;

	list_pts_reordered.push_back(*i_ep_begin_);

	for (int i=0;i!=1;i++){
		Pts_Iterator it1=list_pts_new_.begin();
		Pts_Iterator it2 = it1;
		std::advance(it2, 1);
		if (CGAL::squared_distance(*i_ep_begin_, *it1) < CGAL::squared_distance(*i_ep_begin_, *it2) ){
			list_pts_reordered.push_back(*it1);
			list_pts_reordered.push_back(*it2);
		}
		else{
			list_pts_reordered.push_back(*it2);
			list_pts_reordered.push_back(*it1);
		}
	}

	list_pts_reordered.push_back(*i_ep_end_);

	// Reorder the list of point

	ipe_it_pt_1 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_1, 0);
	IPE_new_1.pts_ = *ipe_it_pt_1;

	ipe_it_pt_2 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_2, 1);
	IPE_new_2.pts_ = *ipe_it_pt_2;

	ipe_it_pt_3 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_3, 2);
	IPE_new_3.pts_ = *ipe_it_pt_3;

	ipe_it_pt_4 = list_pts_reordered.begin();
	std::advance(ipe_it_pt_4, 3);
	IPE_new_4.pts_ = *ipe_it_pt_4;


	// Be careful  NTD 14/09/2016
	Char_Iterator ipe_it_ch = list_chars_.begin();

	for (Pts_Iterator it1 = list_pts_reordered.begin(); it1!=list_pts_reordered.end(); ++it1){
		int idx_= 0;
		for (Pts_Iterator it2 = list_pts_.begin(); it2!=list_pts_.end(); ++it2){
			if (*it1 == *it2){
				ipe_it_ch = list_chars_.begin();
				std::advance(ipe_it_ch, idx_);
				list_chars_reordered.push_back(*ipe_it_ch);}
			idx_++;
		}
	}

	if (VERY_VERBOSE == 1){
		std::cout << "############     list_chars_reordered    ############" << std::endl;
		std::cout << "-----   ";
		for (Char_Iterator it = list_chars_reordered.begin(); it != list_chars_reordered.end(); it++){
			std::cout<< * it <<"   -----   ";
		}
		std::cout <<"\n";
	}

	Pts_Iterator pts_ib(list_pts_reordered.begin());
	Char_Iterator char_ib(list_chars_reordered.begin());
	if (list_chars_reordered.size() > list_pts_reordered.size()){

		if ( *pts_ib == *(boost::next(pts_ib))){
			if ( *(boost::next(boost::next(pts_ib))) != *(boost::next(boost::next(boost::next(pts_ib)))) ){

				// Case: A    A    B    C  --->  1 --- 2 --- 1 --- 2 --- 2 --- 1
				list_chars_reordered.erase(char_ib);
				char_ib = list_chars_reordered.begin();
				list_chars_reordered.erase(char_ib);
			}
			else
			{
				// Case: A    A    B    B  --->  1 --- 2 --- 1 --- 2 --- 1 --- 2 --- 1 --- 2
				list_chars_reordered.erase(char_ib);
				for (int it =0; it!=3; it++){
					char_ib = list_chars_reordered.begin();
					list_chars_reordered.erase(char_ib);
				}
			}
		}

		else{
			if ( *(boost::next(pts_ib)) == *(boost::next(boost::next(pts_ib))) ){
				if (*char_ib == *(boost::next(char_ib))){
					// Case: A    B    B    C  --->  1 --- 1 --- 2 --- 1 --- 2 --- 2
					std::advance(char_ib,1);
					char_ib = list_chars_reordered.begin();
					std::advance(char_ib,3);
					list_chars_reordered.erase(char_ib);
				}
				else{
					// Case: C    B    B    A  --->  2 --- 1 --- 2 --- 1 --- 2 --- 1
					std::advance(char_ib,2);
					char_ib = list_chars_reordered.begin();
					std::advance(char_ib,2);
					list_chars_reordered.erase(char_ib);

				}
			}

		}
	}

	if (VERY_VERBOSE == 1){
		std::cout << "############     list_chars_reordered    ############" << std::endl;
		std::cout << "-----   ";
		for (Char_Iterator it = list_chars_reordered.begin(); it != list_chars_reordered.end(); it++){
			std::cout<< * it <<"   -----   ";
		}
		std::cout <<"\n";
	}
	// Reorder the list of char

	ipe_it_ch_1 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_1, 0);
	IPE_new_1.Ell_name_ = *ipe_it_ch_1;

	ipe_it_ch_2 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_2, 1);
	IPE_new_2.Ell_name_ = *ipe_it_ch_2;

	ipe_it_ch_3 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_3, 2);
	IPE_new_3.Ell_name_ = *ipe_it_ch_3;

	ipe_it_ch_4 = list_chars_reordered.begin();
	std::advance(ipe_it_ch_4, 3);
	IPE_new_4.Ell_name_ = *ipe_it_ch_4;

	list_IP_Ell_out.push_back(IPE_new_1);
	list_IP_Ell_out.push_back(IPE_new_2);
	list_IP_Ell_out.push_back(IPE_new_3);
	list_IP_Ell_out.push_back(IPE_new_4);

	if (VERY_VERBOSE == 1){
		int type_print_out = 1;

		switch (type_print_out){
		case 1:
		{
			std::cout << " " << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_pt_1 <<"   ---   "<< *ipe_it_pt_2 <<"   ---   "
					<< *ipe_it_pt_3 <<"   ---   "<< *ipe_it_pt_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;
			std::cout << "############     On the ellipse    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_ch_1 <<"   ---   "<< *ipe_it_ch_2 <<"   ---   "
					<< *ipe_it_ch_3 <<"   ---   "<< *ipe_it_ch_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;

		}
		break;

		case 2:
		{
			std::cout << " " << std::endl;
			std::cout << "############     Order of points    ############" << std::endl;
			std::cout << "-----   "<< *ipe_it_pt_1 <<"   ---   "<< *ipe_it_pt_2 <<"   ---   "
					<< *ipe_it_pt_3 <<"   ---   "<< *ipe_it_pt_4 <<"   ------" << std::endl;
			std::cout << "################################################" << std::endl;
			std::cout << " " << std::endl;

		}
		break;
		}
	}

	return list_IP_Ell_out;
}
} // end of namespace CGAL
#endif
