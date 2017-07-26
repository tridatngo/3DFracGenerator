/*
 * MyRotation.hh
 *
 *  Created on: Sep 5, 2016
 *      Author: ngotr
 *
 * MyRotationMatrix[\Theta, V]: Gives the 3D rotation matrix for a counterclockwise rotation around the 3D vector V.
 * MyRotation: Return the result of 3D rotation of a point.
 */

#ifndef INCLUDE_MYROTATION_HH_
#define INCLUDE_MYROTATION_HH_

// Transformation of the polygon
#include <CGAL/Aff_transformation_3.h>
#include "Vector_3_Operations.hh"
#include "vectors_are_concident.hh"
namespace CGAL{

template <typename K, typename NT>
CGAL::Aff_transformation_3<K> MyRotationMatrix (CGAL::Vector_3<K> &normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis
	Vector_3 rotvect = CGAL::MyCross<K>(normal_vect, e3);
	Vector_3 u = CGAL::Normalize<K,NT>(rotvect);

	// NTD - 25/10/2016
	//NT theta = CGAL::MyAngle<Kernel, NT_MP>(normal_vect,e3);
	NT theta = CGAL::MyAngle<K, NT>(normal_vect,e3);
	NT c = cos(CGAL::to_double(theta));
	NT s = sin(CGAL::to_double(theta));

	Aff_transformation_3 rotMat(
			u[0]*u[0]*(1-c) +      c,  u[1]*u[0]*(1-c) - u[2]*s,  u[2]*u[0]*(1-c) + u[1]*s,  transvect[0],
			u[0]*u[1]*(1-c) + u[2]*s,  u[1]*u[1]*(1-c) +      c,  u[2]*u[1]*(1-c) - u[0]*s,  transvect[1],
			u[0]*u[2]*(1-c) - u[1]*s,  u[1]*u[2]*(1-c) + u[0]*s,  u[2]*u[2]*(1-c) +      c,  transvect[2],
			1);
	return rotMat;
}

template <typename K, typename NT>
CGAL::Aff_transformation_3<K> MyInverseRotationMatrix (CGAL::Vector_3<K> &normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis
	Vector_3 rotvect = CGAL::MyCross<K>(e3, normal_vect);
	Vector_3 u = CGAL::Normalize<K,NT>(rotvect);

	//NT theta = CGAL::MyAngle<Kernel, NT_MP>(e3, normal_vect);

	// NTD - 25/10/2016
	//NT theta = CGAL::MyAngle<Kernel, NT_MP>(normal_vect,e3);
	NT theta = CGAL::MyAngle<K, NT>(normal_vect,e3);
	NT c = cos(CGAL::to_double(theta));
	NT s = sin(CGAL::to_double(theta));

	Aff_transformation_3 rotMat(
			u[0]*u[0]*(1-c) +      c,  u[1]*u[0]*(1-c) - u[2]*s,  u[2]*u[0]*(1-c) + u[1]*s,  transvect[0],
			u[0]*u[1]*(1-c) + u[2]*s,  u[1]*u[1]*(1-c) +      c,  u[2]*u[1]*(1-c) - u[0]*s,  transvect[1],
			u[0]*u[2]*(1-c) - u[1]*s,  u[1]*u[2]*(1-c) + u[0]*s,  u[2]*u[2]*(1-c) +      c,  transvect[2],
			1);
	return rotMat;
}

template <typename K, typename NT>
CGAL::Point_3<K> PointRotation_Eps (CGAL::Point_3<K> &p, NT & eps_, CGAL::Vector_3<K> & normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;
	typedef typename CGAL::Point_3<K>			   	Point_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis
	Vector_3 normal_vect_norm = CGAL::Normalize<K,NT>(normal_vect);

	CGAL::Aff_transformation_3<K> rotMat;
	Point_3 prot, pout;
	NT pout_x, pout_y, pout_z;

	if (CGAL::vectors_are_concident<K,NT>(normal_vect_norm,e3)){
		pout = p;
	}
	else
	{
		Vector_3 rotvect = CGAL::MyCross<K>(normal_vect, e3);
		Vector_3 u = CGAL::Normalize<K,NT>(rotvect);
		//if ( normal_vect != e3 ){
		rotMat = CGAL::MyRotationMatrix<K,NT> (normal_vect, transvect);
		prot = rotMat(p);

		// If the rotation causes Points to be eps_ close of zero, set those Points to 0.
		if ( CGAL::abs(prot.x()) < eps_){
			pout_x = 0;
		}
		else
			pout_x = prot.x();

		if ( CGAL::abs(prot.y()) < eps_)
			pout_y = 0;
		else
			pout_y = prot.y();

		if ( CGAL::abs(prot.z()) < eps_)
			pout_z = 0;
		else
			pout_z = prot.z();

		pout = Point_3(pout_x, pout_y, pout_z);
	}

	return pout;
}


template <typename K, typename NT>
CGAL::Point_3<K> PointRotation (CGAL::Point_3<K> &p, CGAL::Vector_3<K> & normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;
	typedef typename CGAL::Point_3<K>			   	Point_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along thz z-axis
	Vector_3 normal_vect_norm = CGAL::Normalize<K,NT>(normal_vect);

	CGAL::Aff_transformation_3<K> rotMat;
	Point_3 prot, pout;
	NT pout_x, pout_y, pout_z;

	if (CGAL::vectors_are_concident<K,NT>(normal_vect_norm,e3)){
		pout = p;
	}
	else
	{
		rotMat = CGAL::MyRotationMatrix<K,NT> (normal_vect, transvect);
		prot = rotMat(p);

		pout = prot;
	}

	return pout;
}

template <typename K, typename NT>
CGAL::Point_3<K> PointInverseRotation (CGAL::Point_3<K> &p, CGAL::Vector_3<K> & normal_vect, CGAL::Vector_3<K> & transvect){

	typedef typename CGAL::Aff_transformation_3<K> 	Aff_transformation_3;
	typedef typename CGAL::Vector_3<K>			   	Vector_3;
	typedef typename CGAL::Point_3<K>			   	Point_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis
	Vector_3 normal_vect_norm = CGAL::Normalize<K,NT>(normal_vect);

	CGAL::Aff_transformation_3<K> rotMat;
	Point_3 prot, pout;
	NT pout_x, pout_y, pout_z;

	if (CGAL::vectors_are_concident<K,NT>(normal_vect_norm,e3)){
		pout = p + transvect;
	}
	else
	{
		rotMat = CGAL::MyInverseRotationMatrix<K,NT> (normal_vect, transvect);
		prot = rotMat(p);

		pout = prot;
	}

	return pout;
}

} // end of namespace

#endif /* INCLUDE_MYROTATION_HH_ */
