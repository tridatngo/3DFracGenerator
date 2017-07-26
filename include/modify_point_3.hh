/*
 * modify_point_3.hh
 *
 *  Created on: Sep 21, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_MODIFY_POINT_3_HH_
#define INCLUDE_MODIFY_POINT_3_HH_

#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>
#include <list>
#include <include/myglobal_functions.hh>

#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Kernel/global_functions.h>
#include <list>

namespace CGAL {

template<typename PT, typename NT>
PT zero_p3(PT &p, NT &eps){

	NT pmod_x,pmod_y, pmod_z;
	if (CGAL::abs(p.x()) < eps){
		pmod_x = 0.E+00;
	}
	else{
		pmod_x = p.x();
	}
	if (CGAL::abs(p.y()) < eps){
		pmod_y = 0.E+00;
	}
	else{
		pmod_y = p.y();
	}
	if (CGAL::abs(p.z()) < eps){
		pmod_z = 0.E+00;
	}
	else{
		pmod_z = p.z();
	}

	return PT(pmod_x, pmod_y, pmod_z);
}
}

#endif /* INCLUDE_MODIFY_POINT_3_HH_ */
