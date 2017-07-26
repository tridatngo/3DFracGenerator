/*
 * print_dfn_bb.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 *
 *      print out to *.vtu file the dfn bouding box
 */


#ifndef PRINT_DFN_BB_HH_
#define PRINT_DFN_BB_HH_

namespace CGAL{

template <typename SurfaceMesh_, typename DFN_BB_>
void print_dfn_bb (DFN_BB_ &DFN_BB){

	std::string fname;

	fname = "output/DFN_BB_X_L";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.X_L, false, true);

	fname = "output/DFN_BB_X_U";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.X_U, false, true);

	fname = "output/DFN_BB_Y_L";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.Y_L, false, true);

	fname = "output/DFN_BB_Y_U";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.Y_U, false, true);

	fname = "output/DFN_BB_Z_L";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.Z_L, false, true);

	fname = "output/DFN_BB_Z_U";
	CGAL::output_off_vtu <SurfaceMesh_>(fname, DFN_BB.Z_U, false, true);

}
}

#endif /* PRINT_DFN_BB_HH_ */
