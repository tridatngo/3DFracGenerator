#include "include/Exact_predicates_inexact_constructions_kernel_ntd.h"
#include <CGAL/Surface_mesh.h>

#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Kernel/global_functions.h>

#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>     // std::next & std::prev

// Intersection
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                    Point_3;
typedef Kernel::Line_3                     Line_3;
typedef Kernel::Segment_3                  Segment_3;

typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
typedef CGAL::Inverse_index<std::list<Point_3>::iterator> Inverse_Index;

#ifndef CGAL_ORDER_INTERSECTIONPOINTS_3D_HH
#define CGAL_ORDER_INTERSECTIONPOINTS_3D_HH

namespace CGAL {

void Order_IntersectionPoins3D_3pt (const Point_3 &first, const Point_3 &second, const Point_3 &third){
	if (CGAL::collinear_are_strictly_ordered_along_line(first, second, third)){
		std::cout << "The point (" << first << ") lies between two points (" << second <<") and (" << third <<")" << std::endl;}
	else
		if (CGAL::collinear_are_strictly_ordered_along_line(second, first , third)){
			std::cout << "The point (" << second << ") lies between two points (" << first <<") and (" << third <<")" << std::endl;}
		else
			std::cout << "The point (" << third << ") lies between two points (" << first <<") and (" << second <<")" << std::endl;

} // end of namespace
