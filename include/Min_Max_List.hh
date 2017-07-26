/*
 * Min_list.hh
 *
 *  Created on: 8 août 2016
 *      Author: ngotr
 *      Return Min and max of a list
 */

#include <list>
#include <iterator>

#ifndef MIN_MAX_LIST_HH_
#define MIN_MAX_LIST_HH_

template<typename NumType_>
NumType_ Min_List (std::list<NumType_> &first){

	typedef typename std::list<NumType_>::iterator Iterator_;

	NumType_ min_ = *first.begin();
	for (Iterator_ iter=first.begin(); iter!=first.end(); iter++){
		if (*iter <= min_){
			min_ = *iter;
		}
	}

	return min_;
}


template<typename NumType_>
NumType_ Max_List (std::list<NumType_> &first){

	typedef typename std::list<NumType_>::iterator Iterator_;

	NumType_ max_ = *first.begin();
	for (Iterator_ iter=first.begin(); iter!=first.end(); iter++){
		if (*iter >= max_){
			max_ = *iter;
		}
	}
	return max_;
}

#endif /* MIN_MAX_LIST_HH_ */
