/*
 * transpose_a_list.hh
 *
 * To transpose a list
 *
 *  Created on: 9 août 2016
 *      Author: ngotr
 */

#ifndef TRANSPOSE_A_LIST_HH_
#define TRANSPOSE_A_LIST_HH_

template<typename ListType_>

ListType_ transpose_a_list(ListType_ &list){
	typedef typename ListType_::iterator Iterator;
	ListType_ list_new;

	int idx = 0;

	for (Iterator it = list.begin(); it != list.end(); it++){

		Iterator iter = list.begin();
		std::advance(iter, list.size()-idx-1);

		list_new.push_back(*iter);
		idx++;
	}

	return list_new;
}

#endif /* TRANSPOSE_A_LIST_HH_ */
