/*
 * export_polygon_points.hh
 *
 *  Created on: 26 août 2016
 *      Author: ngotr
 */

#ifndef EXPORT_POLYGON_POINTS_HH_
#define EXPORT_POLYGON_POINTS_HH_


#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <sstream>

namespace CGAL{

template <typename K, typename SM_>
void export_polygon_points(SM_ sm, std::string outputName){

	typedef typename K::Point_3                   	Point_3;
	typedef typename SM_::Vertex_index 		   		vertex_descriptor;
	typedef typename SM_::Vertex_iterator 			vertex_iterator;

	/*
	std::fstream    outfile;
	FILE * fp;
	fp = fopen (outputName.c_str(), "w");

	// header
	fprintf (fp,"%i  0  0 0 0 \n", sm.num_vertices());

	// print points

	int pindex(1);
	for (vertex_iterator it = sm.vertices_begin(); it != sm.vertices_end(); it++ ){
		fprintf (fp, "%i7", pindex);
		char buffer_0 [13], buffer_1 [13], buffer_2 [13];
		sprintf (buffer_0, "%11.6E", sm.point(*it)[0]);
		sprintf (buffer_1, "%11.6E", sm.point(*it)[1]);
		sprintf (buffer_2, "%11.6E", sm.point(*it)[2]);

		fprintf (fp, "%13.13s %13.13s %13.13s", buffer_0, buffer_1, buffer_2);
		fprintf (fp,"\n");
		pindex++;
	}

	fclose(fp);
	*/

	std::ofstream os(outputName.c_str());
	os << sm.num_vertices() << " 0  0 0 0 \n";

	int pindex(1);
		for (vertex_iterator it = sm.vertices_begin(); it != sm.vertices_end(); it++ ){
			os << pindex << "  " << std::uppercase << std::scientific << sm.point(*it) <<std::endl;
			pindex++;
		}
	os.close();
}

} // end of namespace
#endif /* EXPORT_POLYGON_POINTS_HH_ */
