/*
 * remove_close_points.hh
 *
 *  Created on: 18 août 2016
 *      Author: ngotr
 */

#ifndef REMOVE_CLOSE_POINTS_HH_
#define REMOVE_CLOSE_POINTS_HH_

#include <include/points_are_close.hh>

#define VERBOSE_NTD 0

template<typename SM_, typename PT_, typename NT_>
void remove_close_points (SM_ &poly, std::list<PT_> &list, std::list<int> &seg_of_endpoints){

	unsigned int begin = 0;
	unsigned int end = poly.num_vertices() -1;
	typedef typename SM_::Vertex_index 		   vertex_descriptor;

	// Begin_segment
	vertex_descriptor vb_seg_b(begin);
	vertex_descriptor ve_seg_b(begin+1);
	vertex_descriptor ve_seg_b_0(0);

	if (VERBOSE_NTD == 1){
		std::cout << "*list.begin() = " << *list.begin() << std::endl;
		std::cout << "poly.point(vb_seg_b) = " << poly.point(vb_seg_b) << std::endl;
		std::cout << "poly.point(ve_seg_b) = " << poly.point(ve_seg_b) << std::endl;
	}

	if (begin==poly.num_vertices()-1)
		ve_seg_b = ve_seg_b_0;

	if ( points_are_close< PT_, NT_>( *list.begin(), poly.point(vb_seg_b) ) || points_are_close< PT_, NT_>(*list.begin(), poly.point(ve_seg_b)) ){
		list.erase(list.begin());
		//*list.begin() = poly.point(vb_seg_b);
	}

	typename std::list<PT_>::iterator prev_end = list.begin();
	std::advance(prev_end, list.size()-1);

	if ( points_are_close< PT_, NT_>( *prev_end, poly.point(vb_seg_b) ) || points_are_close< PT_, NT_>(*prev_end, poly.point(ve_seg_b)) ){
		list.erase(prev_end);
		//*prev_end = poly.point(vb_seg_b);
	}

	// End_segment

	vertex_descriptor vb_seg_e(end);
	vertex_descriptor ve_seg_e(end+1);
	vertex_descriptor ve_seg_e_0(0);

	if (end==poly.num_vertices()-1)
		ve_seg_e = ve_seg_e_0;

	if ( points_are_close< PT_, NT_>( *list.begin(), poly.point(vb_seg_e) ) || points_are_close< PT_, NT_>(*list.begin(), poly.point(ve_seg_e)) ){
		list.erase(list.begin());
		//*list.begin() = poly.point(vb_seg_e);
	}

	prev_end = list.begin();
	std::advance(prev_end, list.size()-1);

	if ( points_are_close< PT_, NT_>( *prev_end, poly.point(vb_seg_e) ) || points_are_close< PT_, NT_>(*prev_end, poly.point(ve_seg_e)) ){
		list.erase(prev_end);
		//*prev_end = poly.point(vb_seg_e);
	}
}

template<typename SM_, typename PT_, typename NT_>
SM_ remove_close_points_from_poly (SM_ &poly, std::list<PT_> &list, NT_ &tol_length){

	unsigned int begin = 0;
	unsigned int end = poly.num_vertices() -1;
	SM_ poly_new;
	typedef typename SM_::Vertex_index 		   vertex_descriptor;

	bool no_polybegin(false), no_polyend(false);

	// Begin_segment
	vertex_descriptor vb_seg_b(begin);
	vertex_descriptor ve_seg_b(begin+1);
	vertex_descriptor ve_seg_b_0(0);

	if (begin==poly.num_vertices()-1)
		ve_seg_b = ve_seg_b_0;

	typename std::list<PT_>::iterator prev_end = list.begin();
	std::advance(prev_end, list.size()-1);

	if (points_are_close_custom< PT_, NT_>( *list.begin(), poly.point(vb_seg_b), tol_length ) ||
			points_are_close_custom< PT_, NT_>( *prev_end, poly.point(vb_seg_b), tol_length ) ){
		no_polybegin = true;
	}

	/*
	if (points_are_close_custom< PT_, NT_>(*list.begin(), poly.point(ve_seg_b), tol_length)  ||
			points_are_close_custom< PT_, NT_>(*prev_end, poly.point(ve_seg_b), tol_length) ){
		no_polyend = true;
	}
	 */

	// End_segment

	vertex_descriptor vb_seg_e(end);
	vertex_descriptor ve_seg_e(end+1);
	vertex_descriptor ve_seg_e_0(0);

	if (end==poly.num_vertices()-1)
		ve_seg_e = ve_seg_e_0;

	prev_end = list.begin();
	std::advance(prev_end, list.size()-1);

	if (VERBOSE_NTD == 1){

		std::cout << "*list.begin() = " << *list.begin() << std::endl;
		std::cout << "poly.point(vb_seg_b) = " << poly.point(vb_seg_b) << std::endl;
		std::cout << "poly.point(ve_seg_b) = " << poly.point(ve_seg_b) << std::endl;

		std::cout << "*list.prev_end() = " << *prev_end << std::endl;
		std::cout << "poly.point(vb_seg_e) = " << poly.point(vb_seg_e) << std::endl;
		std::cout << "poly.point(ve_seg_e) = " << poly.point(ve_seg_e) << std::endl;
	}
	if ( points_are_close_custom< PT_, NT_>(*list.begin(), poly.point(vb_seg_e), tol_length) ||
			points_are_close_custom< PT_, NT_>( *prev_end, poly.point(vb_seg_e), tol_length) ){
		no_polyend = true;
	}

	/*
	if ( points_are_close_custom< PT_, NT_>(*list.begin(), poly.point(ve_seg_e), tol_length) ||
			points_are_close_custom< PT_, NT_>( *prev_end, poly.point(ve_seg_e), tol_length) ){
		no_polyend = true;
	}
	*/

	vertex_descriptor vb_poly(0);
	vertex_descriptor ve_poly(poly.num_vertices() -1);

	if (!no_polybegin){
		poly_new.add_vertex(poly.point(vb_poly));
	}

	unsigned ibegin = 0;
	for( ; ibegin < end-1; ++ibegin) {
		vertex_descriptor vb_seg1(ibegin);
		vertex_descriptor ve_seg1(ibegin+1);

		poly_new.add_vertex(poly.point(ve_seg1));
	}

	if (!no_polyend){
		poly_new.add_vertex(poly.point(ve_poly));
	}

	return poly_new;
}

#endif /* REMOVE_CLOSE_POINTS_HH_ */
