/*
 * isotropic_triangulation_sm.hh
 *
 *  Created on: 16 août 2016
 *      Author: ngotr
 */

/*
 * This function is in order to do isotropic triangulation of surface meshes
 */
#ifndef ISOTROPIC_TRIANGULATION_SM_HH_
#define ISOTROPIC_TRIANGULATION_SM_HH_

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
//#include <CGAL/Polygon_mesh_processing/remesh.h>
#include "Polygon_mesh_processing/remesh_ntd.h"

#include <CGAL/Polygon_mesh_processing/border.h>

namespace CGAL{

template <typename SurfaceMesh_>
void isotropic_triangulation_sm (SurfaceMesh_ &poly, double &target_edge_length, unsigned int &nb_iter){

	// Triangulating the surface mesh poly
	CGAL::Polygon_mesh_processing::triangulate_faces(poly);

	// Refinement
	CGAL::Polygon_mesh_processing::isotropic_remeshing(
			faces(poly),
			target_edge_length,
			poly,
			CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
	.protect_constraints(true) //i.e. protect border, here
	);
}

} // end of namespace CGAL

#endif /* ISOTROPIC_TRIANGULATION_SM_HH_ */
