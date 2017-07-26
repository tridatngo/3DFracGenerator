/*
 * export_inp_files_from_ell_list.hh
 *
 *  Created on: Sep 8, 2016
 *      Author: ngotr
 */

#ifndef EXPORT_INP_FILES_FROM_ELL_LIST_HH_
#define EXPORT_INP_FILES_FROM_ELL_LIST_HH_

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>
#include <include/myglobal_functions.hh>
#include <include/reserve_orientation_of_clockwise_polygon.hh>

namespace CGAL{

template <typename Ell_, typename Surface_mesh>
void export_inp_files_from_ell_list(std::list<Ell_> & ell_list, std::string fname_, const int &digits){

	Ell_			cEll_;			// Current ellipse
	typedef	typename std::list<Ell_>::iterator Ell_Iterator;
	typedef typename Surface_mesh::Vertex_iterator 			vertex_iterator;

	for (int it = 0; it !=ell_list.size(); it++){
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		// Convert all polygon into counterclockwise-oriented polygon
		//cEll_ = CGAL::reverse_polygon_point_order(cEll_);

		std::string fname;
		fname = fname_ + CGAL::int2str_setw(it,digits) + ".inp";
		std::ofstream cE_name_ss(fname.c_str()) ;

		cE_name_ss << cEll_.E_mesh.num_vertices() << " 0  0 0 0 \n";

		int pindex(1);
		for (vertex_iterator it = cEll_.E_mesh.vertices_begin(); it != cEll_.E_mesh.vertices_end(); it++ ){
			cE_name_ss << pindex << "  " << std::uppercase << std::scientific << cEll_.E_mesh.point(*it) <<std::endl;
			pindex++;
		}
		cE_name_ss.close();
	}


}

} // end of namespace

#endif /* EXPORT_INP_FILES_FROM_ELL_LIST_HH_ */
