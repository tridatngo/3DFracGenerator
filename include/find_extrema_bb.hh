/*
 * find_extrema_bb.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 */

/*
 * Return the lexicographically smallest and largest coordinates
 * of the bounding box
 */

#ifndef FIND_EXTREMA_BB_HH_
#define FIND_EXTREMA_BB_HH_

#define VERBOSE_NTD 0

namespace CGAL{

template <typename Point_3, typename SurfaceMesh_, typename DFN_BB_>
std::list<Point_3> find_extrema_bb (DFN_BB_ &bbox){

	SurfaceMesh_ X_L(bbox.X_L);
	SurfaceMesh_ X_U(bbox.X_U);
	SurfaceMesh_ Y_L(bbox.Y_L);
	SurfaceMesh_ Y_U(bbox.Y_U);
	SurfaceMesh_ Z_L(bbox.Z_L);
	SurfaceMesh_ Z_U(bbox.Z_U);

	std::list<Point_3> output;

	Point_3 p_min, p_max, p_xl, p_yl, p_zl, p_xu, p_yu, p_zu;

	p_xl = X_L.point(*(X_L.vertices_begin()));
	p_yl = Y_L.point(*(Y_L.vertices_begin()));
	p_zl = Z_L.point(*(Z_L.vertices_begin()));

	p_xu = X_U.point(*(X_U.vertices_begin()));
	p_yu = Y_U.point(*(Y_U.vertices_begin()));
	p_zu = Z_U.point(*(Z_U.vertices_begin()));

	p_min = Point_3(p_xl.x(), p_yl.y(), p_zl.z());
	p_max = Point_3(p_xu.x(), p_yu.y(), p_zu.z());

	output.push_back(p_min);
	output.push_back(p_max);

	if (VERBOSE_NTD == 1){
		std::cout << "P_min of the BB :"<< p_min << std::endl;
		std::cout << "P_max of the BB :"<< p_max << std::endl;
	}

	return output;
}
}

#endif /* FIND_EXTREMA_BB_HH_ */
