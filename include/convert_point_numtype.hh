/*
 * convert_point_numtype.hh
 *
 *  Created on: 10 août 2016
 *      Author: ngotr
 */

#ifndef CONVERT_POINT_NUMTYPE_HH_
#define CONVERT_POINT_NUMTYPE_HH_

namespace CGAL {

template<typename NT_, typename MP_NT_>
MP_NT_ NT_to_NT_MP(NT_ &cp){

	return MP_NT_ (cp.hx(), cp.hy(), cp.hz(), cp.hw());
}

template<typename MP_NT_, typename NT_>
NT_ NT_MP_to_NT(MP_NT_ &cp){

	return NT_ (cp.hx(), cp.hy(), cp.hz(), cp.hw());
}

}

#endif /* CONVERT_POINT_NUMTYPE_HH_ */
