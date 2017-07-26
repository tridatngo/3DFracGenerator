/*
 * vectors_are_concident.hh
 *
 *  Created on: Sep 6, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_VECTORS_ARE_CONCIDENT_HH_
#define INCLUDE_VECTORS_ARE_CONCIDENT_HH_

namespace CGAL{

template <typename K, typename NT_>
bool vectors_are_concident (CGAL::Vector_3<K> first, CGAL::Vector_3<K> second)
{
	static NT_ eps_(1E-16);

	return (CGAL::compare(eps_,
			CGAL::max(
					CGAL::max(CGAL::abs(first[0]-second[0]),CGAL::abs(first[1]-second[1])),
					CGAL::abs(first[2]-second[2]))) > 0);
}

} // end of namespace

#endif /* INCLUDE_VECTORS_ARE_CONCIDENT_HH_ */
