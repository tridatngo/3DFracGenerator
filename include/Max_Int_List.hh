/*
 * Max_Int_List.hh
 *
 *  Created on: 8 août 2016
 *      Author: ngotr
 */
#include <list>
#include <iterator>

#ifndef MAX_INT_LIST_HH_
#define MAX_INT_LIST_HH_

int Max_Int_List (std::list<int> &first){
	int max_ = *first.begin();
	for (std::list<int>::iterator iter=first.begin(); iter!=first.end(); iter++){
		if (*iter >= max_){
			max_ = *iter;
		}
	}

	return max_;
}

#endif /* MAX_INT_LIST_HH_ */
