/*
 * My_Bounding_Box.hh
 *
 *  Created on: Sep 6, 2016
 *      Author: ngotr
 *
 * Create a axis-aligned bounding box around a surface mesh
 * If the surface mesh plane is also axis-aligned, a positive width will be set.
 */

#ifndef INCLUDE_MY_BOUNDING_BOX_HH_
#define INCLUDE_MY_BOUNDING_BOX_HH_

#include <CGAL/Iso_cuboid_3.h>

namespace CGAL{

//typedef CGAL::Surface_mesh<Point_3>
template < typename K, typename NT_>
CGAL::Iso_cuboid_3<K> My_Bounding_Box(CGAL::Surface_mesh< CGAL::Point_3<K> > & sm){

	typedef typename CGAL::Point_3<K> Point_3;
	typedef typename CGAL::Surface_mesh<Point_3> Surface_mesh;
	typedef typename CGAL::Iso_cuboid_3<K> Iso_cuboid_3;
	typedef typename Surface_mesh::Vertex_index vertex_descriptor;

	Iso_cuboid_3 isoC;
	vertex_descriptor vh(0);
	NT_ xmin(sm.point(vh)[0]), xmax(sm.point(vh)[0]), ymin(sm.point(vh)[1]), ymax(sm.point(vh)[1]),
		zmin(sm.point(vh)[2]), zmax(sm.point(vh)[2]);

	unsigned int i = 0, end = sm.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vh(i);

		if (xmin >= sm.point(vh)[0]){
			xmin = sm.point(vh)[0];}

		if (xmax <= sm.point(vh)[0]){
			xmax = sm.point(vh)[0];}

		if (ymin >= sm.point(vh)[1]){
			ymin = sm.point(vh)[1];}

		if (ymax <= sm.point(vh)[1]){
			ymax = sm.point(vh)[1];}

		if (zmin >= sm.point(vh)[2]){
			zmin = sm.point(vh)[2];}

		if (zmax <= sm.point(vh)[2]){
			zmax = sm.point(vh)[2];}
	}

	isoC = Iso_cuboid_3(xmin, ymin, zmin, xmax, ymax, zmax);

	return isoC;
}
}

#endif /* INCLUDE_MY_BOUNDING_BOX_HH_ */
