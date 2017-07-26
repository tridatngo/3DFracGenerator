/*
 * All_bool_true.hh
 *
 *  Created on: 5 août 2016
 *      Author: ngotr
 */

#ifndef ALL_BOOL_TRUE_HH
#define ALL_BOOL_TRUE_HH

#include <list>
#include <iterator>

// Function to check if all elements of boolean list are true
bool All_bool_true (std::list<bool> first){
	bool AllTrue=true;
	for(std::list<bool>::iterator it=first.begin(); it!=first.end();++it){
		if(*it==false)
		{
			AllTrue=false;
			break;
		}
	}
	return AllTrue;
}

#endif /* ALL_BOOL_TRUE_HH */
