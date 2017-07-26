/*
 * Max_List.hh
 *
 *  Created on: 8 août 2016
 *      Author: ngotr
 */

#ifndef MAX_LIST_HH_
#define MAX_LIST_HH_

template<typename NumType_>
NumType_ Max_Int_List (std::list<NumType_> &first){

	typedef typename std::list<NumType_>::iterator Iterator_;

	NumType_ max_ = *first.begin();
	for (Iterator_ iter=first.begin(); iter!=first.end(); iter++){
		if (*iter >= max_){
			max_ = *iter;
		}
	}
	return max_;
}

#endif /* MAX_LIST_HH_ */
