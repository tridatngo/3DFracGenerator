/*
 * list_seg_index_unique.hh
 *
 *  Created on: 9 août 2016
 *      Author: ngotr
 */

#ifndef LIST_SEG_INDEX_UNIQUE_HH_
#define LIST_SEG_INDEX_UNIQUE_HH_

#include "points_are_close.hh"
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>

/*
template <typename PT_, typename NT_>
bool points_are_close (PT_ first, PT_ second)
{
	NT_ eps_(1E-3);

	return ( CGAL::is_positive(eps_ - CGAL::max(
					CGAL::max(CGAL::abs(first[0]-second[0]),CGAL::abs(first[1]-second[1])),
					CGAL::abs(first[2]-second[2]))) );
}
 */

template <typename PT_, typename NT_>
// Function to check if two points are near or coincident
void list_points_unique (std::list<PT_> &input)
{
	PT_ point_ref, point_com;
	typedef typename std::list<PT_>::iterator Iterator_;

	int it_(0);
	int test_num(0);

	for (Iterator_ iter = input.begin(); iter !=  input.end(); iter++ ){
		PT_ num_ref = *iter;
		int n = 1;
		Iterator_ it_begin = iter;
		std::advance(it_begin, 1);

		while( it_begin != input.end() && n < input.size()){

			Iterator_ it_begin = iter;
			std::advance(it_begin, n);
			if(points_are_close<PT_, NT_>(*it_begin, num_ref) ==1){
				input.erase(it_begin);
			}
			else
				n++;
		}
	}
}

template <typename PT_, typename NT_>
bool in_points_list (PT_ &point_, std::list<PT_> &input, NT_ &eps_closepts){
	bool bool_in_pl(false);
	typedef typename std::list<PT_>::iterator Iterator_;

	unsigned int i = 1, end = input.size();
	for( ; i < end; ++i) {
		Iterator_ it_pt = input.begin();
		std::advance(it_pt, i);
		if (points_are_close_custom<PT_, NT_>(point_, *it_pt, eps_closepts)){
			bool_in_pl = true;
			break;
		}
	}

	return bool_in_pl;
}

template <typename PT_, typename NT_>
void sorted_list_points_unique (std::list<PT_> &input, NT_ &eps_)
{
	std::list<PT_> output;
	PT_ point_ref, point_com;
	typedef typename std::list<PT_>::iterator Iterator_;
	output.push_back( *input.begin() );

	unsigned int i = 1, end = input.size();
	for( ; i < end; ++i) {
		Iterator_ it_to_add = input.begin();
		std::advance(it_to_add, i);

		if (!points_are_close_custom<PT_, NT_>(*it_to_add, *boost::prior(it_to_add), eps_)){
			output.push_back(*it_to_add);
		}
	}

	input = output;
}

template <typename PT_, typename NT_>
void sorted_list_points (std::list<PT_> &input, std::list<PT_> &point_list, NT_ &eps_, NT_ &eps_closepts)
{
	std::list<PT_> output;
	PT_ point_ref, point_com;
	typedef typename std::list<PT_>::iterator Iterator_;
	output.push_back( *input.begin() );

	unsigned int i = 1, end = input.size();
	for( ; i < end; ++i) {
		Iterator_ it_to_add = input.begin();
		std::advance(it_to_add, i);

		if (!points_are_close_custom<PT_, NT_>(*it_to_add, *boost::prior(it_to_add), eps_)){
			if (i == end -1){
				// The last point --> to add
				output.push_back(*it_to_add);
			}
			else{
				if (!points_are_close_custom<PT_, NT_>(*it_to_add, *boost::next(it_to_add), eps_)){
					output.push_back(*it_to_add);
				}
				else{
					if (in_points_list<PT_, NT_>(*it_to_add, point_list, eps_closepts)){
						output.push_back(*it_to_add);
					}
				}
			}
		}
		else{
			if (in_points_list<PT_, NT_>(*it_to_add, point_list, eps_closepts)){
				output.push_back(*it_to_add);
			}
		}
	}

	input = output;
}

template <typename PT_, typename NT_>
void list_seg_index_unique (std::list<PT_> &input, std::list<int> &ind_input)
{
	PT_ point_ref, point_com;
	typedef typename std::list<PT_>::iterator Iterator_;

	int it_(0);
	int test_num(0);

	CGAL::Inverse_index<Iterator_ > I_idx;

	for (Iterator_ iter = input.begin(); iter !=  input.end(); iter++ ){
		PT_ num_ref = *iter;
		int n = 1;

		Iterator_ it_begin = iter;
		std::advance(it_begin, 1);

		std::list<int>::iterator it_ind_begin = ind_input.begin();
		std::advance(it_ind_begin, 1 + I_idx.operator[](iter));


		// TODO Remove the duplicated points (0,0,0) to avoid error when comparing this point with input.end() which are normally very small

		while( it_begin != input.end() && n < input.size()){

			Iterator_ it_begin = iter;
			std::advance(it_begin, n);

			std::list<int>::iterator it_ind_begin = ind_input.begin();
			std::advance(it_ind_begin, n + I_idx.operator[](iter));

			if(points_are_close<PT_, NT_>(*it_begin, num_ref) ==1){
				input.erase(it_begin);
				ind_input.erase(it_ind_begin);
			}
			else
				n++;
		}
	}
}

#endif /* LIST_SEG_INDEX_UNIQUE_HH_ */
