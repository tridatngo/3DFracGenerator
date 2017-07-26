/*
 * get_target_edge_length.hh
 *
 *  Created on: 12 sept. 2016
 *      Author: ngo
 */

#ifndef INCLUDE_GET_TARGET_EDGE_LENGTH_HH_
#define INCLUDE_GET_TARGET_EDGE_LENGTH_HH_

#define HARMONIC_MEAN 1

namespace CGAL{

template <typename K>
double get_target_edge_length (CGAL::Surface_mesh< CGAL::Point_3 <K> > & poly){

	typedef typename CGAL::Point_3<K> 					Point_3;
	typedef typename CGAL::Segment_3<K> 				Segment_3;
	typedef typename CGAL::Surface_mesh< Point_3> 		Surface_mesh;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename Surface_mesh::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 		vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 		vertex_range;

	double dL0(0);
	unsigned int i = 0, end = poly.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vb_seg1(i);
		vertex_descriptor ve_seg1(i+1);
		vertex_descriptor ve_seg1_0(0);

		if (i==end-1)
			ve_seg1 = ve_seg1_0;

		Segment_3 seg_ = Segment_3( poly.point(vb_seg1), poly.point(ve_seg1));
		if (seg_.squared_length () >= dL0){
			dL0 = CGAL::to_double(seg_.squared_length ());
		}
	}

	dL0 = CGAL::sqrt(CGAL::to_double(dL0));

	return dL0;
}

template < typename K, typename Ell_>
std::list<Ell_> update_target_edge_length (std::list<Ell_> &ell_list) {

	typedef typename std::list<Ell_>::iterator Ell_iterator;
	std::list<Ell_> out_list;
	Ell_ ell_;

	for (Ell_iterator iter = ell_list.begin(); iter != ell_list.end(); iter++){
		ell_ = *iter;
		ell_.E_target_edge_length = CGAL::get_target_edge_length<K>(ell_.E_mesh);
		out_list.push_back(ell_);
	}
	return out_list;
}

template <typename K>
double get_target_edge_length_rm_pts (CGAL::Surface_mesh< CGAL::Point_3 <K> > & poly, double &ratio){

	typedef typename CGAL::Point_3<K> 					Point_3;
	typedef typename CGAL::Segment_3<K> 				Segment_3;
	typedef typename CGAL::Surface_mesh< Point_3> 		Surface_mesh;
	typedef typename std::list<Point_3>::iterator 		Pts_Iterator;
	typedef typename Surface_mesh::Vertex_index 		vertex_descriptor;
	typedef typename Surface_mesh::Vertex_iterator 		vertex_iterator;
	typedef typename Surface_mesh::Vertex_range 		vertex_range;

	double dL0(0), Lbar(0);
	unsigned int i = 0, end = poly.number_of_vertices();
	for( ; i < end; ++i) {
		vertex_descriptor vb_seg1(i);
		vertex_descriptor ve_seg1(i+1);
		vertex_descriptor ve_seg1_0(0);

		if (i==end-1)
			ve_seg1 = ve_seg1_0;

		Segment_3 seg_ = Segment_3( poly.point(vb_seg1), poly.point(ve_seg1));
		if (seg_.squared_length () >= dL0){
			dL0 = CGAL::to_double(seg_.squared_length ());
		}


		if (HARMONIC_MEAN){
			// Harmonic mean
			if (CGAL::sqrt(CGAL::to_double(seg_.squared_length ())) > 1E-15){
				Lbar = Lbar +  1/CGAL::sqrt(CGAL::to_double(seg_.squared_length ()));
			}
		}
		else{
			// Arithmetic mean
			Lbar = Lbar +  1/CGAL::sqrt(CGAL::to_double(seg_.squared_length ()));
		}
	}

	//std::cout << "Lbar = " << Lbar << std::endl;

	if (HARMONIC_MEAN){
		// Harmonic mean
		Lbar = poly.number_of_vertices()/Lbar;
	}
	else {
		// Arithmetic mean
		Lbar = Lbar/poly.number_of_vertices();
	}

	dL0 = CGAL::sqrt(CGAL::to_double(dL0));
	//std::cout << "poly.number_of_vertices() = " << poly.number_of_vertices() << ", Lbar = " << Lbar << ", dL0 = " << dL0 <<std::endl;

	if  (HARMONIC_MEAN){
		return Lbar;
	}
	else{
		if (dL0 > Lbar * ratio){
			//dL0 = Lbar;
			dL0 = (dL0+Lbar)/2;
		}

		return dL0;
	}
}

template < typename K, typename Ell_>
std::list<Ell_> update_target_edge_length_rm_pts (std::list<Ell_> &ell_list, double ratio) {

	typedef typename std::list<Ell_>::iterator Ell_iterator;
	std::list<Ell_> out_list;
	Ell_ ell_;

	for (Ell_iterator iter = ell_list.begin(); iter != ell_list.end(); iter++){
		ell_ = *iter;
		ell_.E_target_edge_length = CGAL::get_target_edge_length_rm_pts<K>(ell_.E_mesh, ratio);
		out_list.push_back(ell_);
	}
	return out_list;
}

template <typename K>
double get_poly_perimeter (std::list<CGAL::Point_3<K> > & pt_list){

	typedef typename CGAL::Point_3<K> 							Point_3;
	typedef typename std::list<CGAL::Point_3<K> >::iterator   	Pts_Iterator;

	double perim(0);

	unsigned int ip_ = 0;
	for( ; ip_ < pt_list.size(); ++ip_) {
		Pts_Iterator it_p1, it_p2;
		if (ip_ == pt_list.size()-1){
			it_p1 = pt_list.begin();
			std::advance(it_p1, ip_);

			it_p2 = pt_list.begin();
		}
		else{
			it_p1 = pt_list.begin();
			std::advance(it_p1, ip_);
			it_p2 = pt_list.begin();
			std::advance(it_p2, ip_+1);
		}
		double distance_ =  CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(*it_p2, *it_p1)));

		perim = perim + distance_;
	}

	return perim;
}

}

#endif /* INCLUDE_GET_TARGET_EDGE_LENGTH_HH_ */
