/*
 * remove_polygons.hh
 *
 *  Created on: Sep 6, 2016
 *      Author: ngotr
 *
 * Remove disconnected polygons / polygons which do not lie in DFN_BB
 */

#ifndef INCLUDE_REMOVE_POLYGONS_HH_
#define INCLUDE_REMOVE_POLYGONS_HH_

#include <algorithm>
#include "update_polygon_numbering.hh"

using namespace std;

namespace CGAL{

template < typename K, typename Ell_>
std::list<Ell_> remove_polygons(std::list<Ell_> &ell_list, std::list<int> & list_) {

	typedef typename std::list<Ell_>::iterator Ell_iterator;
	std::list<Ell_> out_list;

	for (Ell_iterator iter = ell_list.begin(); iter != ell_list.end(); iter++){
		if (list_.size() != 0 ){

			bool found = (std::find(list_.begin(), list_.end(), (*iter).E_name) != list_.end());
			if (!found){
				out_list.push_back(*iter);
			}
		}
		else{
			out_list.push_back(*iter);
		}
	}

	out_list = CGAL::update_polygon_numbering <K, Ell_>(out_list);

	return out_list;
}

} // end of namespace

#endif /* INCLUDE_REMOVE_POLYGONS_HH_ */
