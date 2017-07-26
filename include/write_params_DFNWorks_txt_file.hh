/*
 * write_params_DFNWorks_txt_file.hh
 *
 *  Created on: Oct 25, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_WRITE_PARAMS_DFNWORKS_TXT_FILE_HH_
#define INCLUDE_WRITE_PARAMS_DFNWORKS_TXT_FILE_HH_

#ifndef VERY_VERBOSE
#define VERY_VERBOSE 0
#endif

namespace CGAL{

template < typename K, typename Ell_, typename NT>
void write_params_DFNWorks_txt_file (
		std::list<Ell_> &ell_list,
		CGAL::Point_3<K> p_min,
		CGAL::Point_3<K> p_max,
		const std::string &filename_) {

	typedef typename std::list<Ell_>::iterator Ell_Iterator;
	typedef typename CGAL::Vector_3<K> Vector_3;
	typedef typename CGAL::Point_3<K> Point_3;

	Vector_3 e3 =  CGAL::UnitVector_3<K>(3); // unit vector along the z-axis

	const char* filename = filename_.c_str();
	std::ofstream    outfile(filename);
	Vector_3 vect_;
	NT x1,y1,z1,x2,y2,z2;

	//std::ofstream    outfile;
	//outfile.open(outfile_.c_str());

	std::list<Ell_> ell_list_ = ell_list;
	double target_edge_length(0.0);

	// Line 1 is the number of polygons
	outfile << ell_list_.size() << std::endl;

	target_edge_length =(*(ell_list_).begin()).E_target_edge_length;

	for (int i = 1; i !=ell_list_.size(); i++){
		Ell_			cEll_1;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list_.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;
		if (target_edge_length < cEll_1.E_target_edge_length){
			target_edge_length = cEll_1.E_target_edge_length;
		}
	}

	//Line 2 is the h scale
	outfile << target_edge_length << "\n";

	// Line 3 numberOfPoints in a Polygon
	outfile << "5\n";

	// Line 4 slope is the rate at which the mesh will get coarser in refinement
	//slope = float(fin.readline())
	outfile << "1\n";

	// Line 5 refine_dist is the distance from the intersections which refinement will be performed
	//refine_dist = float(fin.readline())
	outfile << "0.000001\n";

	// Line 6 is the visualization mode: '1' is on, '0' is off.
	outfile << "0\n";

	// Line 7 is the filename of the avs for the polygon
	outfile << "polys.inp\n";

	// Line 8 is the filename of the avs for the intersections
	outfile << "intersections.inp\n";

	int inum(1);
	for (Ell_Iterator it = ell_list_.begin(); it != ell_list_.end(); it++){
		vect_ = (*it).E_normalvect;
		vect_ = CGAL::Normalize<K, NT>(vect_);
		/*  Find Angle of Rotation*/
		if(VERY_VERBOSE){
			std::cout << "vect_ = " << vect_ << std::endl;
		}

		NT theta = CGAL::MyAngle<K, NT_MP>(vect_,e3);
		NT c = cos(CGAL::to_double(theta));
		NT s = sin(CGAL::to_double(theta));

		/* Convert angle to degrees*/
		theta = CGAL::Rad2Deg(theta);

		if(VERY_VERBOSE){
			std::cout << "theta = " << theta << std::endl;
		}

		/* Find angle to Rotate into xy plane */
		Vector_3 rotvect = CGAL::MyCross<K>(vect_, e3);
		Vector_3 u = CGAL::Normalize<K,NT>(rotvect);
		/*
		NT dk(CGAL::max(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-12){
			dk = CGAL::min(p_min[0], p_max[0]);
		}
		 */

		NT dk(CGAL::min(p_min[0], p_max[0]));

		if ( CGAL::abs(dk) < 1E-12){
			dk = CGAL::max(p_min[0], p_max[0]);
		}

		x1 = - dk * u[0]; x2 = dk * u[0];
		y1 = - dk * u[1]; y2 = dk * u[1];
		z1 = - dk * u[2]; z2 = dk * u[2];

		//outfile << inum << " " << std::setprecision(12) << -theta << " " << (*it).E_target_edge_length << " " << x1 << " " << y1 << " " << z1
		//		<< " " << x2 << " " << y2 << " " << z2 << "\n";

		//outfile << inum << " " << vect_[0] << " " <<  vect_[1] << " " << vect_[2] << std::endl;

		outfile << inum << " " << std::setprecision(12) << -theta << " " << x1 << " " << y1 << " " << z1
					<< " " << x2 << " " << y2 << " " << z2 << "\n";

		inum++;
	}

	outfile.close();
}
}

#endif /* INCLUDE_WRITE_PARAMS_DFNWORKS_TXT_FILE_HH_ */
