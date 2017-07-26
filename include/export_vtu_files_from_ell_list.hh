/*
 * export_vtu_files_from_ell_list.hh
 *
 *  Created on: Sep 8, 2016
 *      Author: ngotr
 */

#ifndef EXPORT_VTU_FILES_FROM_ELL_LIST_HH_
#define EXPORT_VTU_FILES_FROM_ELL_LIST_HH_

#include "myglobal_functions.hh"

namespace CGAL{

template <typename Ell_, typename Surface_mesh>
void export_vtu_files_from_ell_list(std::list<Ell_> & ell_list, std::string fname_){

	typedef	typename std::list<Ell_>::iterator Ell_Iterator;

	for (int it = 0; it !=ell_list.size(); it++){
		Ell_			cEll_;			// Current ellipse
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		std::string fname;
		//std::ostringstream cE_name_ss;
		//cE_name_ss << cEll_.E_name;
		fname = fname_ + CGAL::int2str(it);
		CGAL::output_off_vtu <Surface_mesh>(fname, cEll_.E_mesh, false, true);
	}
}
} // end of namespace

#endif /* EXPORT_VTU_FILES_FROM_ELL_LIST_HH_ */
