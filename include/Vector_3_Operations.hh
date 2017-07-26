/*
 * Vector_3_Operations.hh
 *
 *  Created on: Sep 5, 2016
 *      Author: ngotr
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Norm[v]	give the norm of v
 * Normalize[v]	give a unit vector in the direction of v
 * Standardize[v]	shift v to have zero mean and unit sample variance
 * Standardize[v,f1]	shift v by f1[v] and scale to have unit sample variance
 * UnitVector[n,k] gives the n-dimensional unit vector in the k^(th) direction.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef INCLUDE_VECTOR_3_OPERATIONS_HH_
#define INCLUDE_VECTOR_3_OPERATIONS_HH_

namespace CGAL{

template < typename K, typename NT>
NT Norm (const CGAL::Vector_3<K> & vect)
{
	float square_ = CGAL::to_double(CGAL::square(vect[0]) + CGAL::square(vect[1]) + CGAL::square(vect[2]));
	return CGAL::sqrt( square_);
}

template < typename K, typename NT>
CGAL::Vector_3<K> Normalize (const CGAL::Vector_3<K> & vect)
{
	typedef typename CGAL::Vector_3<K> Vector_3;
	Vector_3 vect_norm(0,0,0);
	/*
	 if (CGAL::Norm<K, NT>(vect) == 0) {
		std::cout << "Vector " << vect << " is a zero vector. Return the unit vector (0,0,1)." <<std::endl;
	}
	else{
	 */
	if (CGAL::Norm<K, NT>(vect) != 0) {
	vect_norm = CGAL::Vector_3<K>(	vect[0] / CGAL::Norm<K, NT>(vect),
									vect[1] / CGAL::Norm<K, NT>(vect),
									vect[2] / CGAL::Norm<K, NT>(vect));
	}
	return vect_norm;
}

template < typename K>
CGAL::Vector_3<K>  UnitVector_3 (const int & n)
{
	typedef typename CGAL::Vector_3<K> Vector_3;
	Vector_3 unit_vect_3;
	switch (n){
		case 1:
		{
			unit_vect_3 = CGAL::Vector_3<K> (1,0,0);
		}
		break;

		case 2:
		{
			unit_vect_3 = CGAL::Vector_3<K> (0,1,0);
		}
		break;

		case 3:
		{
			unit_vect_3 = CGAL::Vector_3<K> (0,0,1);
		}
		break;
	}

	return unit_vect_3;
}

template < typename K>
CGAL::Vector_3<K> MyCross (const CGAL::Vector_3<K> & v1_, const CGAL::Vector_3<K> & v2_)
{
	return CGAL::Vector_3<K> (	v1_[1]*v2_[2] - v1_[2]*v2_[1],
								v1_[2]*v2_[0] - v1_[0]*v2_[2],
								v1_[0]*v2_[1] - v1_[1]*v2_[0]);
}

template < typename K, typename NT>
NT MyAngle (const CGAL::Vector_3<K>  & v1_, const CGAL::Vector_3<K> & v2_)
{
	/*
	NT angle_ = acos (CGAL::to_double(CGAL::Normalize<K, NT>(v1_) *  CGAL::Normalize<K, NT>(v2_)));
	return angle_;
	*/
	// NTD 10/11/2016
	if ( CGAL::Norm<K,NT>( CGAL::MyCross<K>(v1_, v2_) ) == 0 ){
		return 0;
	}
	else{
		return acos (CGAL::to_double(CGAL::Normalize<K, NT>(v1_) *  CGAL::Normalize<K, NT>(v2_)));
	}
}

} // end of namespace

#endif /* INCLUDE_VECTOR_3_OPERATIONS_HH_ */
