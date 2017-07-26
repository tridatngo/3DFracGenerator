/*
 * create_lagrit_scripts.hh
 *
 *  Created on: Sep 9, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_CREATE_LAGRIT_SCRIPTS_HH_
#define INCLUDE_CREATE_LAGRIT_SCRIPTS_HH_


#include<list>
#include<string>
#include<iostream>
#include<sstream>
#include "myglobal_functions.hh"

#ifndef VERY_VERBOSE
#define VERY_VERBOSE 0
#endif

namespace CGAL{

template<typename Ell_>
void create_lagrit_scripts(	const std::string &paramsfile,
		const std::string &linefile,
		const std::string &polyfile,
		const std::string &outfile_,
		Ell_ &ell){
	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	outfile << "# LaGriT Script \n";

	outfile << "# Read the parameter input file.\n";
	outfile << "infile / " << paramsfile << std::endl;

	outfile << "# Name the input file that contain the intersection lines.\n";

	outfile << "# define / LINE_FILE / LINEFILE\n";
	outfile << "define / LINE_FILE / " << linefile  << std::endl;

	outfile << "# Name the input file that contain the polygons.\n";

	outfile << "# define / POLY_FILE / POLYFILE\n";
	outfile << "define / POLY_FILE / " << polyfile  << std::endl;

	outfile << "# Define parameters:\n";

	// input ID, mo_poly_work, mo_line_work
	outfile << "# Read the first polygon\n";
	outfile << "cmo / create / mo_poly_work / / / tri\n";

	outfile << "# Read in line and polygon files:\n";
	outfile << "read / avs / POLY_FILE / mo_poly_work\n";
	outfile << "read / avs / LINE_FILE / mo_line_work\n";

	outfile << "cmo / create / mo_poly_pts / / / triplane \n";
	outfile << "cmo / select / mo_poly_pts \n";
	outfile << "copypts / mo_poly_pts / mo_poly_work \n";
	outfile << "cmo / select / mo_poly_work \n";

	/*
	cmo / create / cmotri1 / / / tri
	read / avs / poly_1_1.inp / cmotri1
	define / ICOLOR / 1
	 */

	//outfile << "compute / distance_field / mo_final / mo_line_work / dfield\n";

	outfile << "#--# Infile: Triangulate a set of point and set imt, itetclr, itp\n";

	outfile << "triangulate/counterclockwise \n";
	outfile << "cmo / setatt / mo_poly_work / imt / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "cmo / setatt / mo_poly_work / itetclr / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "resetpts / itp\n";
	outfile << "#\n";
	outfile << "# Create parent/child points at material interfaces\n";
	outfile << "#\n";
	outfile << "settets\n";
	outfile << "#\n";
	outfile << "# Reset the itetclr values of nodes to be the same as the element imt\n";
	outfile << "#\n";
	outfile << "resetpts / imt \n";

	outfile << "quality\n";

	outfile << "#\n";
	outfile << "# Create isotropic mesh.\n";
	outfile << "#\n";

	outfile << "# This set of arguments will remove degenerate elements from a mesh by merging nodes that have\n"
			"# the same coordinate values ( within 1.e-9). Mesh edges longer than ell.E_target_edge_length \n"
			"# will be bisected. \n";

	// here 1.05 is a factor of safety which ensure that all boundary segment will not be bisected.
	//outfile << "massage / "<< ell.E_target_edge_length * 1.05 << " / 1E-9 / 1.E-9 / nosmooth / strictmergelenth\n";

	// Version 15_09_2016 _
	// toldamage is not specified, no node annihilation will take place.
	// See http://lagrit.lanl.gov/docs/commands/MASSAGE.html
	//outfile << "massage / "<< ell.E_target_edge_length * 1.05 << " / 1E-9 \n";

	// Version 20_09_2016
	outfile << "massage / "<< ell.E_target_edge_length << " / 1E-9 / 1.E-9 / nosmooth / strictmergelenth\n";


	outfile << "# associate the pset name allpts with all points\n";
	outfile << "pset/allpts/seq/1,0,0/\n";

	outfile << "# associate the name pboundary with the points whose type field(itp1)\n";
	outfile << "# has value greater than or equal to 10 (these would be boundary nodes)\n";

	outfile << "pset/pboundary/attribute/itp/1,0,0/10/ge\n";

	outfile << "# associate the pset name psmooth with all point except boundary nodes\n";
	outfile << "pset / psmooth / not / allpts pboundary\n";

	outfile << "# massage / H_SCALE / 1.e-5 / 1.e-5 / pset get pref / &\n";
	outfile << "#nosmooth / strictmergelenth\n";

	outfile << "massage / "<< ell.E_target_edge_length << " / 1E-9 / 1E-9 / nosmooth / strictmergelenth\n";

	outfile << "massage / 1E+20 / " << ell.E_target_edge_length * 0.5 << " / " << ell.E_target_edge_length*0.5 << " / &\n";
	outfile << "nosmooth / strictmergelenth\n";

	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";

	outfile << "\n";
	outfile << "# cmo poly final \n";
	outfile << "cmo / create / mo_poly_final / / / triplane \n";
	outfile << "cmo / select / mo_poly_final \n";
	outfile << "copypts / mo_poly_final / mo_poly_work / \n";
	outfile << "pset/allpts/seq/1,0,0/ \n";
	outfile << "pset/pboundary/attribute/itp/1,0,0/10/ge \n";
	outfile << "pset / psmooth / not / allpts pboundary \n";
	outfile << "rmpoint / pset get pboundary \n";

	outfile << "rmpoint / compress \n";

	outfile << "copypts / mo_poly_final / mo_poly_pts \n";
	outfile << "cmo / select / mo_poly_final \n";
	outfile << "resetpts / itp \n";

	outfile << "cmo / setatt / mo_poly_final / imt / 1 0 0 / ID \n";
	outfile << "cmo / setatt / mo_poly_final / itp / 1 0 0 / 0 \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";


	// NTD - Modification 02/11/2016
	/*
	outfile << "connect \n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	*/
	// NTD - End of modification 02/11/2016

	outfile << "trans/ 1 0 0 / zero / xyz \n";

	outfile << "cmo / setatt / mo_poly_final / zic / 1 0 0 / 0.0 \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";
	outfile << "connect \n";

	outfile << "trans / 1 0 0 / original / xyz \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";

	/*
	outfile << "dump / gmv / output_tri_0.gmv / mo_poly_work\n";
	outfile << "dump / avs / output_tri_0.inp / mo_poly_work\n";
	 */

	outfile << "# Rotate\n";
	outfile << "rotateln / 1 0 0 / nocopy / X1, Y1, Z1 / X2, Y2, Z2 / THETA / 0.,0.,0.,/\n";

	// NTD 22/09/16
	outfile << "compute / distance_field / mo_poly_final / mo_line_work / dfield\n";
	// End

	outfile << "cmo / delete / mo_poly_pts \n";
	outfile << "cmo / delete / mo_poly_work \n";
	outfile << "cmo / delete / mo_line_work \n";
	outfile << "\n";

	//outfile << "rotateln / 1 0 0 / nocopy / X1, Y1, Z1 / X2, Y2, Z2 / THETA / X1 Y1 Z1 /\n";

	// NTD 22/09/2016: Do not need translation
	//outfile << "# Translation\n";

	//outfile << "trans/ 1 0 0 / X1, Y1, Z1 / X2, Y2, Z2 / \n";
	//outfile << "trans/ 1 0 0 / X1, Y1, Z1 / 0.,0.,0. / \n";
	//outfile << "trans/ 1 0 0 / 0.,0.,0. / X2, Y2, Z2 / \n";

	//outfile << "cmo / printatt / mo_poly_work / -xyz- / minmax\n";

	//outfile << "recon 1\n";

	//outfile << "resetpts / itp\n";

	/*
	outfile << "cmo / addatt / mo_poly_work / unit_area_normal / xyz / vnorm\n";
	outfile << "cmo / addatt / mo_poly_work / scalar / xnorm ynorm znorm / vnorm\n";
	outfile << "cmo / DELATT / mo_poly_work / vnorm\n";
	 */

	outfile << "cmo / setatt / mo_poly_final / imt / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "cmo / setatt / mo_poly_final / itetclr / 1 0 0 / " << ell.Parent_E_name << "\n";

	outfile << "dump / gmv / OUTFILE_GMV / mo_poly_final\n";
	outfile << "dump / avs / OUTFILE_AVS / mo_poly_final\n";
	outfile << "dump / stor / tri_poly_fracture_" << ell.E_name + 1 << " / mo_poly_final / ascii\n";

	outfile << "quality\n";
	//outfile << "cmo / delete / mo_poly_final\n";
	outfile << "cmo / status / brief\n";

	outfile << "finish\n";
	outfile.close();

	/*
	 *
	 * Merge all the triangulated polygons into a single MO
	 *
	addmesh / merge / cmo_all / cmotri1 / cmotri2

	dump / gmv / output_all_0.gmv / cmo_all
	dump / avs / output_all_0.inp / cmo_all

	define / EPS_FILTER / 1.e-4
	filter / 1 0 0 / EPS_FILTER
	rmpoint / compress

	 * SORT can affect a_b attribute

	sort / cmo_all / index / ascending / ikey / imt xic yic zic
	reorder / cmo_all / ikey
	cmo / DELATT / cmo_all / ikey
	resetpts / itp
	boundary_components

	dump / gmv / output_all.gmv / cmo_all
	dump / avs / output_all.inp / cmo_all
	 */
}

template<typename Ell_>
void create_lagrit_scripts_onefrac(	const std::string &paramsfile,
		const std::string &linefile,
		const std::string &polyfile,
		const std::string &outfile_,
		Ell_ &ell){
	/* Go through the list and write out parameter file for each polygon to be an input file for LaGriT */

	std::ofstream    outfile;
	outfile.open(outfile_.c_str());

	outfile << "# LaGriT Script \n";

	outfile << "# Read the parameter input file.\n";
	outfile << "infile / " << paramsfile << std::endl;

	outfile << "# Name the input file that contain the intersection lines.\n";

	outfile << "# define / LINE_FILE / LINEFILE\n";
	outfile << "define / LINE_FILE / " << linefile  << std::endl;

	outfile << "# Name the input file that contain the polygons.\n";

	outfile << "# define / POLY_FILE / POLYFILE\n";
	outfile << "define / POLY_FILE / " << polyfile  << std::endl;

	outfile << "# Define parameters:\n";

	// input ID, mo_poly_work, mo_line_work
	outfile << "# Read the first polygon\n";
	outfile << "cmo / create / mo_poly_work / / / tri\n";

	outfile << "# Read in line and polygon files:\n";
	outfile << "read / avs / POLY_FILE / mo_poly_work\n";
	//outfile << "read / avs / LINE_FILE / mo_line_work\n";

	outfile << "cmo / create / mo_poly_pts / / / triplane \n";
	outfile << "cmo / select / mo_poly_pts \n";
	outfile << "copypts / mo_poly_pts / mo_poly_work \n";
	outfile << "cmo / select / mo_poly_work \n";

	outfile << "#--# Infile: Triangulate a set of point and set imt, itetclr, itp\n";

	outfile << "triangulate/counterclockwise \n";
	outfile << "cmo / setatt / mo_poly_work / imt / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "cmo / setatt / mo_poly_work / itetclr / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "resetpts / itp\n";
	outfile << "#\n";
	outfile << "# Create parent/child points at material interfaces\n";
	outfile << "#\n";
	outfile << "settets\n";
	outfile << "#\n";
	outfile << "# Reset the itetclr values of nodes to be the same as the element imt\n";
	outfile << "#\n";
	outfile << "resetpts / imt \n";

	outfile << "quality\n";

	outfile << "#\n";
	outfile << "# Create isotropic mesh.\n";
	outfile << "#\n";

	outfile << "# This set of arguments will remove degenerate elements from a mesh by merging nodes that have\n"
			"# the same coordinate values ( within 1.e-9). Mesh edges longer than ell.E_target_edge_length \n"
			"# will be bisected. \n";

	// Version 20_09_2016
	outfile << "massage / "<< ell.E_target_edge_length << " / 1E-9 / 1.E-9 / nosmooth / strictmergelenth\n";


	outfile << "# associate the pset name allpts with all points\n";
	outfile << "pset/allpts/seq/1,0,0/\n";

	outfile << "# associate the name pboundary with the points whose type field(itp1)\n";
	outfile << "# has value greater than or equal to 10 (these would be boundary nodes)\n";

	outfile << "pset/pboundary/attribute/itp/1,0,0/10/ge\n";

	outfile << "# associate the pset name psmooth with all point except boundary nodes\n";
	outfile << "pset / psmooth / not / allpts pboundary\n";

	outfile << "# massage / H_SCALE / 1.e-5 / 1.e-5 / pset get pref / &\n";
	outfile << "#nosmooth / strictmergelenth\n";

	outfile << "massage / "<< ell.E_target_edge_length << " / 1E-9 / 1E-9 / nosmooth / strictmergelenth\n";

	outfile << "massage / 1E+20 / " << ell.E_target_edge_length * 0.5 << " / " << ell.E_target_edge_length*0.5 << " / &\n";
	outfile << "nosmooth / strictmergelenth\n";

	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";

	outfile << "\n";
	outfile << "# cmo poly final \n";
	outfile << "cmo / create / mo_poly_final / / / triplane \n";
	outfile << "cmo / select / mo_poly_final \n";
	outfile << "copypts / mo_poly_final / mo_poly_work / \n";
	outfile << "pset/allpts/seq/1,0,0/ \n";
	outfile << "pset/pboundary/attribute/itp/1,0,0/10/ge \n";
	outfile << "pset / psmooth / not / allpts pboundary \n";
	outfile << "rmpoint / pset get pboundary \n";

	outfile << "rmpoint / compress \n";

	outfile << "copypts / mo_poly_final / mo_poly_pts \n";
	outfile << "cmo / select / mo_poly_final \n";

	// NTD - Modification 02/11/2016
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	outfile << "smooth / position / esug / pset get psmooth; recon 0;\n";
	// NTD - End of modification 02/11/2016

	outfile << "resetpts / itp \n";

	outfile << "cmo / setatt / mo_poly_final / imt / 1 0 0 / ID \n";
	outfile << "cmo / setatt / mo_poly_final / itp / 1 0 0 / 0 \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";
	outfile << "trans/ 1 0 0 / zero / xyz \n";

	outfile << "cmo / setatt / mo_poly_final / zic / 1 0 0 / 0.0 \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";
	outfile << "connect \n";

	outfile << "trans / 1 0 0 / original / xyz \n";
	outfile << "cmo / printatt / mo_poly_final / -xyz- / minmax \n";

	outfile << "# Rotate\n";
	outfile << "rotateln / 1 0 0 / nocopy / X1, Y1, Z1 / X2, Y2, Z2 / THETA / 0.,0.,0.,/\n";

	// NTD 22/09/16
	//outfile << "compute / distance_field / mo_poly_final / mo_line_work / dfield\n";
	// End

	outfile << "cmo / delete / mo_poly_pts \n";
	outfile << "cmo / delete / mo_poly_work \n";
	//outfile << "cmo / delete / mo_line_work \n";
	outfile << "\n";

	outfile << "cmo / setatt / mo_poly_final / imt / 1 0 0 / " << ell.Parent_E_name << "\n";
	outfile << "cmo / setatt / mo_poly_final / itetclr / 1 0 0 / " << ell.Parent_E_name << "\n";

	outfile << "dump / gmv / OUTFILE_GMV / mo_poly_final\n";
	outfile << "dump / avs / OUTFILE_AVS / mo_poly_final\n";

	outfile << "quality\n";
	outfile << "cmo / status / brief\n";

	outfile << "finish\n";
	outfile.close();
}

template<typename Ell_>
void create_lagrit_scripts_list(std::list<Ell_> ell_list){

	typedef	typename std::list<Ell_>::iterator Ell_Iterator;
	Ell_			cEll_;			// Current ellipse
	int 			digits( CGAL::count_digits(ell_list.size()) );

	// Creates LaGriT script to be run for each polygon

	std::cout << "Writing LaGriT Control Files" << std::endl;

	for (int it = 0; it !=ell_list.size(); it++){
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		/*
		std::string paramsfile("output/parameters/parameters_" + CGAL::int2str_setw(cEll_.E_name,digits) + ".mlgi");
		std::string polyfile("output/lagrit/ellipse_cutoff_" + CGAL::int2str_setw(cEll_.E_name,digits) + ".inp");
		std::string mesh_polyfile_lagrit("output/lagrit/mesh_poly_" + CGAL::int2str_setw(cEll_.E_name,digits) + ".lgi");
		 */

		std::string paramsfile("parameters_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".mlgi");
		//std::string linefile("intersections/intersections_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp");
		std::string linefile("intersections_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp");
		std::string polyfile("ellipse_cutoff_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".inp");
		std::string mesh_polyfile_lagrit("output/lagrit/mesh_poly_" + CGAL::int2str_setw(cEll_.E_name + 1,digits) + ".lgi");

		if (ell_list.size() == 1){
			CGAL::create_lagrit_scripts_onefrac(paramsfile, linefile, polyfile, mesh_polyfile_lagrit, cEll_);
		}
		else{
			CGAL::create_lagrit_scripts(paramsfile, linefile, polyfile, mesh_polyfile_lagrit, cEll_);
		}
	}

	std::cout << "Writing LaGriT Control Files: Complete." << std::endl;
}

} // end of namespace

#endif /* INCLUDE_CREATE_LAGRIT_SCRIPTS_HH_ */
