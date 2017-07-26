/*
 * points_are_close.hh
 *
 *  Created on: 9 août 2016
 *      Author: ngotr
 */

#ifndef POINTS_ARE_NEAR_HH_
#define POINTS_ARE_NEAR_HH_

template <typename PT_, typename NT_>
bool points_are_close (PT_ first, PT_ second)
{
	NT_ eps_(1E-10);

	return (CGAL::is_negative(CGAL::max(
					CGAL::max(CGAL::abs(first[0]-second[0]),CGAL::abs(first[1]-second[1])),
					CGAL::abs(first[2]-second[2]))- eps_) );
}

template <typename PT_, typename NT_>
bool points_are_close_custom (PT_ first, PT_ second, NT_ &tol)
{
	NT_ eps_(tol);

	return (CGAL::is_negative(CGAL::max(
					CGAL::max(CGAL::abs(first[0]-second[0]),CGAL::abs(first[1]-second[1])),
					CGAL::abs(first[2]-second[2]))- eps_) );
}

#endif /* POINTS_ARE_NEAR_HH_ */
