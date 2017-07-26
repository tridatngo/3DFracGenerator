/*
 * Create_DFN_BB.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 *
 *      Create the bounding box of the DFN
 */

#include <list>

#ifndef CREATE_DFN_BB_HH_
#define CREATE_DFN_BB_HH_

namespace CGAL{

template <typename SurfaceMesh_, typename Point_3, typename Output>
Output create_dfn_bb (Point_3 &p_min, Point_3 &p_max){

	// The rectangle of the Lower and upper of the DFN bounding box
	SurfaceMesh_ X_L, X_U, Y_L, Y_U, Z_L, Z_U;
	Output output;

	// X_L
	X_L.add_vertex( Point_3(p_min.x(),p_min.y(),p_min.z()) );
	X_L.add_vertex( Point_3(p_min.x(),p_max.y(),p_min.z()) );
	X_L.add_vertex( Point_3(p_min.x(),p_max.y(),p_max.z()) );
	X_L.add_vertex( Point_3(p_min.x(),p_min.y(),p_max.z()) );

	X_L.add_face(X_L.vertices());

	// X_U
	X_U.add_vertex( Point_3(p_max.x(),p_min.y(),p_min.z()) );
	X_U.add_vertex( Point_3(p_max.x(),p_max.y(),p_min.z()) );
	X_U.add_vertex( Point_3(p_max.x(),p_max.y(),p_max.z()) );
	X_U.add_vertex( Point_3(p_max.x(),p_min.y(),p_max.z()) );

	X_U.add_face(X_U.vertices());

	// Y_L
	Y_L.add_vertex( Point_3(p_min.x(),p_min.y(),p_min.z()) );
	Y_L.add_vertex( Point_3(p_max.x(),p_min.y(),p_min.z()) );
	Y_L.add_vertex( Point_3(p_max.x(),p_min.y(),p_max.z()) );
	Y_L.add_vertex( Point_3(p_min.x(),p_min.y(),p_max.z()) );

	Y_L.add_face(Y_L.vertices());

	// Y_U
	Y_U.add_vertex( Point_3(p_min.x(),p_max.y(),p_min.z()) );
	Y_U.add_vertex( Point_3(p_max.x(),p_max.y(),p_min.z()) );
	Y_U.add_vertex( Point_3(p_max.x(),p_max.y(),p_max.z()) );
	Y_U.add_vertex( Point_3(p_min.x(),p_max.y(),p_max.z()) );

	Y_U.add_face(Y_U.vertices());

	// Z_L
	Z_L.add_vertex( Point_3(p_min.x(),p_min.y(),p_min.z()) );
	Z_L.add_vertex( Point_3(p_max.x(),p_min.y(),p_min.z()) );
	Z_L.add_vertex( Point_3(p_max.x(),p_max.y(),p_min.z()) );
	Z_L.add_vertex( Point_3(p_min.x(),p_max.y(),p_min.z()) );

	Z_L.add_face(Z_L.vertices());

	// Z_U
	Z_U.add_vertex( Point_3(p_min.x(),p_min.y(),p_max.z()) );
	Z_U.add_vertex( Point_3(p_max.x(),p_min.y(),p_max.z()) );
	Z_U.add_vertex( Point_3(p_max.x(),p_max.y(),p_max.z()) );
	Z_U.add_vertex( Point_3(p_min.x(),p_max.y(),p_max.z()) );

	Z_U.add_face(Z_L.vertices());

	output.X_L = X_L;
	output.X_U = X_U;
	output.Y_L = Y_L;
	output.Y_U = Y_U;
	output.Z_L = Z_L;
	output.Z_U = Z_U;

	return output;
}
}

#endif /* CREATE_DFN_BB_HH_ */
