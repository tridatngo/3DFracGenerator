/*
 * update_polygon_numbering.hh
 *
 *  Created on: Sep 6, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_UPDATE_POLYGON_NUMBERING_HH_
#define INCLUDE_UPDATE_POLYGON_NUMBERING_HH_

namespace CGAL{

template < typename K, typename Ell_>
std::list<Ell_> update_polygon_numbering (std::list<Ell_> &ell_list) {

	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	std::list<Ell_> out_list;
	Ell_ ell_;
	int inum(0);

	for (Ell_Iterator iter = ell_list.begin(); iter != ell_list.end(); iter++){
		ell_ = *iter;
		ell_.E_name = inum;
		out_list.push_back(ell_);
		inum++;
	}
	/*
	std::cout << "OK : update_polygon_numbering." << std::endl;
	std::cout << "ell_list.size() = " << ell_list.size()  << std::endl;
	std::cout << "out_list.size() = " << out_list.size()  << std::endl;
	*/
	return out_list;
}

template < typename K, typename Ell_>
std::list<Ell_> update_parent_polygon_numbering (std::list<Ell_> &ell_list) {

	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	std::list<Ell_> out_list;
	Ell_ ell_;
	int inum(0);

	for (Ell_Iterator iter = ell_list.begin(); iter != ell_list.end(); iter++){
		ell_ = *iter;
		ell_.Parent_E_name = inum + 1;
		out_list.push_back(ell_);
		inum++;
	}
	return out_list;
}

} // end of namespace

#endif /* INCLUDE_UPDATE_POLYGON_NUMBERING_HH_ */
