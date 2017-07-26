/*
 * Min_Int_List.hh
 *
 *  Created on: 8 août 2016
 *      Author: ngotr
 *
 * Return the minimum of a integer list
 */
#include <list>
#include <iterator>

#ifndef MIN_INT_LIST_HH_
#define MIN_INT_LIST_HH_

int Min_Int_List (std::list<int> &first){
	int min_ = *first.begin();
	for (std::list<int>::iterator iter=first.begin(); iter!=first.end(); iter++){
		if (*iter <= min_){
			min_ = *iter;
		}
	}

	return min_;
}

#endif /* MIN_INT_LIST_HH_ */
