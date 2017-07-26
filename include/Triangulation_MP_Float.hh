/*
 * Triangulation_MP_Float.hh
 *
 *  Created on: 10 août 2016
 *      Author: ngotr
 */

#ifndef TRIANGULATION_MP_FLOAT_HH_
#define TRIANGULATION_MP_FLOAT_HH_

namespace CGAL {



#include "Exact_predicates_inexact_constructions_kernel_ntd_float.h"
//#include "Exact_predicates_inexact_constructions_kernel_ntd_MP_Float.h"

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>
#include <CGAL/Cartesian_converter.h>

//#include "convert_point_numtype.hh"

//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
//#include <CGAL/Polygon_mesh_processing/refine.h>
//#include <CGAL/Polygon_mesh_processing/fair.h>
//#include <CGAL/Polygon_mesh_processing/remesh.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Float_Kernel;

template <typename MP_FK_, typename FK_, typename Surface_mesh>
Surface_mesh triangulation_mp_float (Surface_mesh &poly_MP_in, double target_edge_length, int nb_iter){

	//typedef CGAL::Exact_predicates_inexact_constructions_kernel Float_Kernel;

	typedef typename FK_::Point_3                    		Point_3_fl;
	typedef CGAL::Surface_mesh<Point_3_fl>        	 		Surface_mesh_fl;
	typedef typename Surface_mesh_fl::Halfedge_iterator		halfedge_iterator_fl;

	typedef typename Surface_mesh::Vertex_index 		   	vertex_descriptor_mp;
	typedef typename Surface_mesh_fl::Vertex_index 		   	vertex_descriptor_fl;
	typedef typename MP_FK_::Point_3		   				Point_3_mp_fl;

	typedef Simple_cartesian<double>  						DK_;

	typedef CGAL::Cartesian_converter<MP_FK_,DK_>                         MP_FK_to_DK;
	typedef CGAL::Cartesian_converter<DK_,FK_>                        	  DK_to_FK;
	typedef CGAL::Cartesian_converter<FK_,DK_>                            FK_to_DK;
	typedef CGAL::Cartesian_converter<DK_,MP_FK_>                         DK_to_MP_FK;

	MP_FK_to_DK to_inexact_1;
	DK_to_FK to_inexact_2;
	FK_to_DK to_exact_1;
	DK_to_MP_FK to_exact_2;

	/*
	typedef CGAL::Cartesian_converter<FK_,MP_FK_>                         FK_to_MP_FK;
	typedef CGAL::Cartesian_converter<MP_FK_,FK_>                         MP_FK_to_FK;

	FK_to_MP_FK to_exact;
	MP_FK_to_FK to_inexact;
	*/

	Surface_mesh poly_MP_out;
	Surface_mesh_fl poly_fl;

	BOOST_FOREACH(vertex_descriptor_mp vd, poly_MP_in.vertices()){
		//poly_fl.add_vertex( to_inexact( poly_MP_in.point(vd)) );
		poly_fl.add_vertex( to_inexact_2(to_inexact_1( poly_MP_in.point(vd)) ));
		std::cout << "Point of poly_2_child_2_converted : " << poly_fl.point(vd) << std::endl;

	}

	poly_fl.add_face(poly_fl.vertices());

	std::string output0("test0.off");
	std::ofstream os0(output0.c_str()) ; os0 << poly_fl ;
	os0.close();

	CGAL::read_off_write_vtk_xml_file(output0);


	// Triangulate the polygon
	CGAL::Polygon_mesh_processing::triangulate_faces(poly_fl);

	std::string output1("test1.off");
	std::ofstream os1(output1.c_str()) ; os1 << poly_fl ;
	os1.close();

	// Refinement
	CGAL::Polygon_mesh_processing::isotropic_remeshing(
			faces(poly_fl),
			target_edge_length,
			poly_fl,
			CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
	.protect_constraints(true)//i.e. protect border, here
	);

	std::cout << "Print to check: OK_2" << std::endl;

	std::string output2("test2.off");
	std::ofstream os2(output2.c_str()) ; os2 << poly_fl ;
	os2.close();

	BOOST_FOREACH(vertex_descriptor_fl vd, poly_fl.vertices()){
		//poly_MP_out.add_vertex( to_exact( poly_fl.point(vd)) );
		poly_MP_out.add_vertex( to_exact_2( to_exact_1( poly_fl.point(vd)) ));
	}

	/*
	for (halfedge_iterator_fl ed_it = poly_fl.halfedges_begin(); ed_it !=poly_fl.halfedges_end(); ed_it++){
		poly_MP_out.add_edge(*ed_it);
	}
	*/

	poly_MP_out.add_face(poly_MP_out.vertices());

	return poly_MP_out;
}

} // end namespace CGAL

#endif /* TRIANGULATION_MP_FLOAT_HH_ */
