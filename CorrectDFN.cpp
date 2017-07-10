//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Kernel option
//#include "include/Exact_predicates_inexact_constructions_kernel_ntd.h"

//#include "include/Exact_predicates_inexact_constructions_kernel_ntd_double.h"
#include "include/Exact_predicates_inexact_constructions_kernel_ntd_float.h"
#include "include/Exact_predicates_inexact_constructions_kernel_ntd_MP_Float.h"
#include "include/Exact_predicates_exact_constructions_kernel_ntd.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
//#include <CGAL/Polygon_mesh_processing/remesh.h>
//#include "include/Polygon_mesh_processing/remesh_ntd.h"

#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>

// Add some boost files
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/include/deref.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/include/begin.hpp>

#include <unistd.h>
#include <ctime>
#include <boost/utility.hpp>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>     /* std::next & std::prev */
#include <stdlib.h>		/* system, NULL, EXIT_FAILURE */
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>

//#include <tuple>
// parse parameters from input file
#include "include/parse_params.hh"
#include "include/parse_parameter_file.hh"
/// toto
//#include "write_c2t3_to_vtk_xml_file.h"
#include "include/read_off_write_vtk_xml_file.h"
#include "include/read_inp_write_vtk_xml_file.hh"
#include "include/read_inp_write_vtk_legacy_file.hh"
#include "include/read_inp_write_dgf_file.hh"
#include "include/write_polys_inp_file.hh"
#include "include/write_params_DFNWorks_txt_file.hh"
#include "include/write_full_mesh_for_DFNWorks.hh"
#include "include/removeNegativeVoronoiCells.hh"

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Cartesian_converter.h>

//#include <CGAL/Simple_cartesian.h>

// Transformation of the polygon
#include <CGAL/Aff_transformation_3.h>
//#include <CGAL/Aff_transformation_tags.h>

// Intersection
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>

// Order of the intersection points
#include "include/order_intersectionpoints_3d.hh"

// Function to check if all elements of boolean list are true
#include "include/All_bool_true.hh"

// Minimum & maximum of an integer list
#include "include/Min_Int_List.hh"
#include "include/Max_Int_List.hh"

// Minimum & maximum of a list of numbers
#include "include/Min_Max_List.hh"

// Transpose a list
#include "include/transpose_a_list.hh"

#include "include/remove_close_points.hh"

//Intersection between two planes on which two polygons lie
#include "include/plane_plane_intersection.hh"
#include "include/list_seg_index_unique.hh"
//#include "include/plane_plane_intersection_convert2MP.hh"
//#include "include/Triangulation_MP_Float.hh"

//#include <write_off_to_vtk_xml_file.h>
#include "include/Output_off_vtu.hh"
#include "include/Output_off_vtu_gmv.hh"
#include "include/export_polygon_points.hh"
#include "include/export_vtu_files_from_ell_list.hh"
#include "include/export_inp_files_from_ell_list.hh"
#include "include/write_params_file_for_lagrit.hh"
#include "include/create_parameter_mlgi_files.hh"
#include "include/create_lagrit_scripts.hh"
#include "include/create_merge_poly_files.hh"
#include "include/redefine_zones.hh"
#include "include/output_ell_to_vtu.hh"
#include "include/getNumLines.hh"

// create the DFN bounding box
#include "include/create_dfn_bb.hh"

// print the bounding box to *.vtu
#include "include/print_dfn_bb.hh"

#include "include/remove_out_of_bb_parts.hh"
#include "include/remove_out_of_bb_parts_rm_pts.hh"

// Cutting off polygons along the intersection
#include "include/cutting_off_polygons.hh"
#include "include/cutting_off_first_polygon.hh"
#include "include/cutting_off_first_polygon_rmpts.hh"
#include "include/cutting_off_two_polygons.hh"

// Triangulating surface meshes
//#include "include/isotropic_triangulation_sm.hh"

// Return the lexicographically smallest and largest coordinates of bb
#include "include/find_extrema_bb.hh"

// Library of messages to print out when running the program
#include "include/libmessages.hh"

// Read fracture parameters input file
#include "include/read_fracs_params_file.hh"
// Vector_3 operations
#include "include/Vector_3_Operations.hh"

// My global functions
#include "include/myglobal_functions.hh"

// My rotation matrix
#include "include/MyRotation.hh"

// My bounding box
#include <CGAL/Iso_cuboid_3.h>
#include <include/points_are_close.hh>
#include "include/My_Bounding_Box.hh"

// Remove disconnected polygons / polygons which do not lie in DFN_BB
#include "include/remove_polygons.hh"

// Get target edge length of the ellipse
#include "include/get_target_edge_length.hh"

// Run lagrit script to triangulate polygons
#include "include/triangulate_polygons.hh"

// Split the intersection file
#include "include/split_intersection_file.hh"

// Check_quality_voronoi_mesh
#include "include/check_quality_voronoi_mesh.hh"
// =====================================================================
// Define macro-variables
// =====================================================================

#define PI_ 3.141592653589793
#define CGAL_PMP_REMESHING_DEBUG

#define CGAL_PMP_REMESHING_VERBOSE
#define CGAL_PMP_REMESHING_VERY_VERBOSE

#define DEBUG_NTD 0
#define PRINT_BBOX 1
#define VERBOSE 0

#ifndef CGAL_REMOVE_POINTS
#define CGAL_REMOVE_POINTS 0
#endif

#ifndef USE_CONVERT_TO_EPECK
#define USE_CONVERT_TO_EPECK 1
#endif

#ifndef USE_EPECK
#define USE_EPECK 0
#endif

// =====================================================================

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

//typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT_MP;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


#if (USE_EPECK)
typedef CGAL::Exact_predicates_exact_constructions_kernel 				Kernel;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel_MP_Float 	Kernel;
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel 			I_Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel 				EPECK;


typedef CGAL::Simple_cartesian<NT_MP> Simple_cartesian;
//typedef CGAL::Simple_cartesian<double> Simple_cartesian;
typedef CGAL::Surface_mesh<Simple_cartesian::Point_3> Mesh;
//typedef Mesh::Vertex_index vertex_descriptor;
//typedef CGAL::Point_3<Kernel>                    Point_3;
typedef Kernel::Point_3                    Point_3;
typedef Kernel::Vector_3                   Vector_3;
typedef Kernel::Plane_3                    Plane_3;
typedef Kernel::Line_3                     Line_3;
typedef Kernel::Segment_3                  Segment_3;
typedef CGAL::Surface_mesh<Point_3>        Surface_mesh;
typedef Surface_mesh::Vertex_index 		   vertex_descriptor;

typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef Surface_mesh::Vertex_iterator vertex_iterator;

typedef std::list<Point_3>::iterator Pts_Iterator;
typedef std::list<Line_3>::iterator Line_Iterator;
typedef std::list<NT_MP>::iterator NT_Iterator;
typedef std::list<int>::iterator Int_Iterator;
typedef std::list<char>::iterator Char_Iterator;
typedef std::list<Vector_3>::iterator Vector_Iterator;

typedef std::list< std::list<Point_3> >::iterator ListPts_Iterator;

typedef CGAL::Inverse_index<Pts_Iterator> Inverse_Index;
typedef CGAL::Iso_cuboid_3<Kernel> Iso_cuboid_3;

typedef std::list< std::list<Point_3> >::iterator List_Pts_Iterator;

//==============================================================================================
// EPECK

typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef EPECK::Point_3                    			Point_3_epeck;
typedef EPECK::Line_3                     			Line_3_epeck;
typedef EPECK::Segment_3                  			Segment_3_epeck;
typedef EPECK::Plane_3                  			Plane_3_epeck;
typedef EPECK::Vector_3								Vector_3_epeck;
typedef EPECK::Intersect_3 							Intersect_3_epeck;
typedef CGAL::Surface_mesh<Point_3_epeck>       	Surface_mesh_epeck;
typedef Surface_mesh_epeck::Vertex_index 			vertex_descriptor_epeck;

typedef std::list<Point_3_epeck>::iterator 			Pts_Iterator_epeck;
typedef std::list<Line_3_epeck>::iterator 			Line_Iterator_epeck;
typedef std::list<Vector_3_epeck>::iterator 		Vector_Iterator_epeck;

typedef CGAL::Cartesian_converter<Kernel,EPECK > Converter;
//==============================================================================================

//#define NT_MP eps_ 1E-10

// Structure of the list of ellipse parsing from input file
struct Ell_list_in{
	// The ellipse will first be constructed on XY plane form E_center and E_xradius, E_yradius
	int 				num_Ell;				// Number of ellipses
	std::list<int> 		E_name;					// Name of the ellipse
	std::list<int> 		E_family;				// Family of the ellipse
	std::list<Point_3> 	E_center;				// Center of the ellipse
	std::list<NT_MP> 	E_xradius, E_yradius;	// x_radius, y_radius
	std::list<Vector_3> E_normalvect;			// Normal vector of the ellipse
};

// Structure of the list of polygon parsing from input file
struct Poly_list_in{
	// The ellipse will first be constructed on XY plane form E_center and E_xradius, E_yradius
	int 							num_Poly;		// Number of polygons
	std::list<int> 					P_name;			// Name of the polygon
	std::list<int> 					P_family;		// Family of the polygon
	std::list< std::list<Point_3> > P_points;		// List of points
	std::list<Vector_3> 			P_normalvect;	// Normal vector of the polygon
};

// Structure of ellipse with surface mesh
struct Ell_{
	int 			E_name;					// Name of the ellipse
	Vector_3 		E_normalvect;			// Normal vector of the ellipse
	Surface_mesh 	E_mesh;                 // Mesh of the ellipse
	double			E_target_edge_length;	// target edgge length of the ellipse
	int 			Parent_E_name;			// Name of the parent ellipse
};

typedef std::list<Ell_> Ell_list;
typedef Ell_list::const_iterator Ell_Iterator;

typedef std::list< std::list<Ell_> > Ell_list_list;
typedef Ell_list_list::const_iterator Ell_list_Iterator;

/*
// Structure of the list ellipse after creating surface mesh
struct Ell_list{
	std::list<int> 			E_name;					// Name of the ellipse
	std::list<Vector_3> 	E_normalvect;			// normal vector
	std::list<Surface_mesh> E_mesh;                  // Mesh of the ellipse
};
 */

// DFN bounding box
struct DFN_BB_{
	Surface_mesh X_L;				// Lower face along X_axis
	Surface_mesh X_U;				// Upper face along X_axis
	Surface_mesh Y_L;				// Lower face along Y_axis
	Surface_mesh Y_U;				// Upper face along Y_axis
	Surface_mesh Z_L;				// Lower face along Z_axis
	Surface_mesh Z_U;				// Upper face along Z_axis
};

// Intersection points and ellipse name
struct IP_Ell
{
	Point_3 pts_;
	char Ell_name_;
};

typedef std::list<IP_Ell> IPE_list;
typedef IPE_list::const_iterator IPE_Iterator;

struct Ell_interst{
	char *Name;						// Name of the ellipse
	std::list<Point_3> IntPoints;	// List of intersection points
};

// Output for the function to determine intersection points between two polygons
struct Output_IP
{
	int num_intersec_1_, num_intersec_2_;
	std::list<Point_3> intersection_points_1_, intersection_points_2_;
	std::list<int> list_seg_interst1_, list_seg_interst2_;

};

struct Output_IPL
{
	bool bool_intersecting_;
	int num_intersec_1_, num_intersec_2_;
	std::list<Point_3> intersection_points_1_, intersection_points_2_;
	std::list<int> list_seg_interst1_, list_seg_interst2_;
	Line_3 intersection_line_;
	Line_3_epeck intersection_line_epeck_;
};

struct Output_Child_Polys
{
	Surface_mesh poly_1_child_1_, poly_1_child_2_, poly_2_child_1_, poly_2_child_2_;
	std::list<Point_3> intersectionline_;
	bool no_polygons_interst_;
};


struct Output_Inside_BB
{
	Surface_mesh poly_inside_bb;
	bool inside_bb;
};

struct InputParams
{
	bool 	hasBlockDensity, hasBlockSpatialParams,
	hasBlockAddingPointsMethod, hasBlockDiscreMethod, hasBlockRemoveClosePoints;
	bool	hasn4,  hasaddingPointsMethod, hasdiscreMethod,
	hasremoveClosePoints, hasmergeCloseIntersectionPoints, hasminLengthRatio, hasdimScaling, hasintRound,
	hasmergeClosePointsRelCritLength;
	bool 	hasBlockCleanUpDirectory, hasBlockDataFile, hasdataFile, hasdesactFracsFile, hasBlockParametersLaGriT,
	hasEPS_FILTER, hasEPS_INT, hasdataForDFNWorks, hascorrectMesh, hasinverseSignTheta, hasparamFilePrecision,
	hasEPS_BOUNDARY;
	int 	cleanupDirectory;
	std::string dataFile, desactFracsFile;
	double EPS_INT, EPS_FILTER, EPS_BOUNDARY, minLengthRatio, mergeClosePointsRelCritLength;
	bool 	hasDFNBB_LowerLeftX,
	hasDFNBB_LowerLeftY,
	hasDFNBB_LowerLeftZ,
	hasDFNBB_UpperRightX,
	hasDFNBB_UpperRightY,
	hasDFNBB_UpperRightZ;
	int 	n4, addingPointsMethod, discreMethod, dimScaling, intRound;
	NT_MP 	DFNBB_LowerLeftX,
	DFNBB_LowerLeftY,
	DFNBB_LowerLeftZ,
	DFNBB_UpperRightX,
	DFNBB_UpperRightY,
	DFNBB_UpperRightZ;
	bool rm_close_pts, merge_close_intersection_pts, dataForDFNWorks, correctMesh, inverseSignTheta, paramFilePrecision;
};

struct stat info;

//using namespace std;

int main(int argc, char* argv[])
{
	time_t start, stop_input, stop_check_out_of_dfnbb, stop_check_intersection, stop_cutting_off, stop_lagrit_run, stop;
	time(&start);
	clock_t begin = clock();

	CGAL::ifp_logo_1();
	std::string fname, cmd, source_dir("/work/irlin104_1/ngotr/WorkSpaceEclipse460/Create_DFN/");

	bool runtime_status(true);

	const char* inputFile = (argc > 1) ? argv[1] : "dfnGen.input";
	const char* inputData = (argc > 2) ? argv[2] : "ellipse_input.dat";

	//std::cout  << "argc = " << argc << std::endl;

	cout << " ==> Correcting the mesh... " << endl;

	// [IMP] Be careful with the precision of number type!!!
	CGAL::Lazy_exact_nt<NT_MP>::set_relative_precision_of_to_double(1E-60);

	InputParams iparams;
	//iparams = CGAL::parse_parameter_file<InputParams>(argc, argv);
	iparams = CGAL::parse_parameter_file<InputParams>(inputFile);

	double EPS_FILTER(1E-4), EPS_INT(1E-5), EPS_BOUNDARY(1E-4);
	if (iparams.hasBlockParametersLaGriT && iparams.hasEPS_FILTER){
		EPS_FILTER = iparams.EPS_FILTER;
	}
	if (iparams.hasBlockParametersLaGriT && iparams.hasEPS_INT){
		EPS_INT = iparams.EPS_INT;
	}
	if (iparams.hasBlockParametersLaGriT && iparams.EPS_BOUNDARY){
		EPS_BOUNDARY = iparams.EPS_BOUNDARY;
	}

	// Path to lagrit executable
	std::string lagrit_path( std::string("/work/irlin104_1/ngotr/Outils/lagrit/lagrit_ulin3.2") );

	// =====================================================================================================
	//                                      Write input file for DFNWorks run
	// =====================================================================================================

	std::system("cp ./output/polymeshes/mesh*.inp .");

	int dimScaling = 1, intRound = 12;
	if (iparams.hasdimScaling){
		dimScaling = iparams.dimScaling;
	}

	// Run lagrit script to triangulate polygons
	bool dataForDFNWorks(false), correctMesh(true);
	if (iparams.hasBlockParametersLaGriT){
		if (iparams.hasdataForDFNWorks){
			dataForDFNWorks = iparams.dataForDFNWorks;
		}
		else{
			std::cout << "** Warning: BlockParametersLaGriT is set but the program failed to parse dataForDFNWorks value.\n"
					<< "** Default \"dataForDFNWorks = false\" will be set." << std::endl;
		}
	}

	if (dataForDFNWorks){
		//std::system("rm -rf part_1.lg");
		//CGAL::removeNegativeVoronoiCells(full_ell);
		std::cout << "Running lagrit < merge_poly_new.lgi ..." << std::endl;

		//cmd = lagrit_path + std::string(" < merge_poly_new.lgi > lagrit_logs/log_merge_poly_new");
		//std::system(cmd.c_str());

		cmd = lagrit_path + std::string(" < merge_rmpts.lgi > lagrit_logs/log_merge_all_new");
		std::system(cmd.c_str());

		std::cout << "Checking quality of Voronoi mesh: " << std::endl;
		bool check_quality_voronoi_mesh = CGAL::check_quality_voronoi_mesh();
		if ( check_quality_voronoi_mesh ){
			std::cout << "CGAL::check_quality_voronoi_mesh() = " << check_quality_voronoi_mesh << std::endl;
		}
		else{
			// Some close but not duplicated points are merged => reduce EPS_FILTER parameter
			std::cout << "CGAL::check_quality_voronoi_mesh() = " << check_quality_voronoi_mesh << std::endl;
			bool check_quality_voronoi_mesh_new = CGAL::check_quality_voronoi_mesh();
			while(!check_quality_voronoi_mesh_new){
				EPS_INT  = EPS_INT * 0.9; // reduce EPS_INT parameter by a factor of 0.9
				EPS_FILTER = EPS_FILTER * 0.9; // reduce EPS_FILTER parameter by a factor of 0.9
				CGAL::create_merge_rmpts_file(EPS_INT, EPS_FILTER);
				// Run merge_rmpts.
				CGAL::run_merge_mesh_scripts(lagrit_path);
				check_quality_voronoi_mesh_new = CGAL::check_quality_voronoi_mesh();
			}
		}
	}
	else{
		// Run merge_rmpts twice
		CGAL::run_merge_mesh_scripts(lagrit_path);

		std::cout << "Checking quality of Voronoi mesh: " << std::endl;
		bool check_quality_voronoi_mesh = CGAL::check_quality_voronoi_mesh();
		if ( check_quality_voronoi_mesh ){
			std::cout << "CGAL::check_quality_voronoi_mesh() = " << check_quality_voronoi_mesh << std::endl;
		}
		else{
			// Some close but not duplicated points are merged => reduce EPS_FILTER parameter
			std::cout << "CGAL::check_quality_voronoi_mesh() = " << check_quality_voronoi_mesh << std::endl;
			bool check_quality_voronoi_mesh_new = CGAL::check_quality_voronoi_mesh();
			while(!check_quality_voronoi_mesh_new){
				EPS_INT  = EPS_INT * 0.9; // reduce EPS_INT parameter by a factor of 0.9
				EPS_FILTER = EPS_FILTER * 0.9; // reduce EPS_FILTER parameter by a factor of 0.9
				CGAL::create_merge_rmpts_file(EPS_INT, EPS_FILTER);
				// Run merge_rmpts.
				CGAL::run_merge_mesh_scripts(lagrit_path);
				check_quality_voronoi_mesh_new = CGAL::check_quality_voronoi_mesh();
			}
		}

		CGAL::redefine_zones(lagrit_path,EPS_FILTER);
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Update "ell.Parent_E_name" for LaGriT DFNWorks - PTDFNTrans run
	std::list<int> parentName_Out_list;
	std::ifstream    parentName_infile("parentName_Out_list.dat");
	std::string line;

	while (!parentName_infile.eof()){
		std::getline(parentName_infile, line);
		std::istringstream ist(line);
		ist.clear();
		int dummy;
		ist >> dummy;
		parentName_Out_list.push_back(dummy);
	}
	parentName_Out_list.sort();
	parentName_Out_list.unique();

	/*
	for (std::list<int>::iterator i = parentName_Out_list.begin(); i!=parentName_Out_list.end(); i++){
		cout << "parentName_Out_list = " << *i << endl;
	}
	*/

	CGAL::write_full_mesh_for_DFNWorks<NT_MP>(std::string("full_mesh.inp"),std::string("full_mesh_DFNWorks.inp"),
			parentName_Out_list,dimScaling, lagrit_path);

	if (dataForDFNWorks){
		CGAL::redefine_zones_dimScaling(lagrit_path,EPS_FILTER,dimScaling);
	}

	// =====================================================================================================
	//                                      Clean up the temporary files
	// =====================================================================================================

	// Clean up # Move individual gmv, avs into folder 'meshes'
	//std::system("rm -rf output/meshes & mkdir output/meshes");
	//std::system("rm -rf output/polymeshes & mkdir output/polymeshes");

	std::system("rm -rf mesh*.inp");

	std::cout << "Writing the full mesh to vtk file... ";


	if (iparams.hasintRound){
		intRound = iparams.intRound;
	}

	//CGAL::read_inp_write_vtk_xml_file( std::string("full_mesh") );
	//CGAL::read_inp_write_vtk_legacy_file( std::string("full_mesh") );
	if (runtime_status){
		//runtime_status = CGAL::read_inp_write_vtk_legacy_file_with_matID( std::string("full_mesh") );
		//runtime_status = CGAL::read_inp_write_vtk_legacy_file_with_matID_dimScaling( std::string("full_mesh"), dimScaling );
		runtime_status = CGAL::read_inp_write_vtk_legacy_file_with_matID_dimScaling_round( std::string("full_mesh"), dimScaling, intRound);
	}

	std::cout << "Complete." << std::endl;

	std::cout << "Writing the full mesh to dgf file... ";
	// VTK -> MSH
	std::system("./Executables/vtk2msh full_mesh");
	// MSH -> DGF
	std::system("./Executables/gmsh2dgf full_mesh.msh");
	std::system("mv full_mesh.dgf full_mesh.dgf_");

	std::string apertureFile(std::string("aperture.dat")), permFile(std::string("permeability.dat"));

	//CGAL::read_inp_write_dgf_file_with_aper( std::string("full_mesh"), std::string("aperture.dat"));
	//CGAL::read_inp_write_dgf_file( std::string("full_mesh"));

	//std::cout << "dimScaling = " << dimScaling << std::endl;

	if (runtime_status){
		//runtime_status = CGAL::read_inp_write_dgf_file_with_aper( std::string("full_mesh"), apertureFile);
		//runtime_status = CGAL::read_inp_write_dgf_file_with_perm( std::string("full_mesh"), permFile);
		//runtime_status = CGAL::read_inp_write_dgf_file_with_perm_aper( std::string("full_mesh"), permFile, apertureFile);
		//runtime_status = CGAL::read_inp_write_dgf_file_with_perm_aper( std::string("full_mesh"), dimScaling, permFile, apertureFile);
		runtime_status = CGAL::read_inp_write_dgf_file_with_perm_aper( std::string("full_mesh"), dimScaling, permFile, apertureFile, intRound);
	}

	std::system("rm -rf full_mesh.dgf_");
	std::system("rm -rf boundary_*");
	std::system("rm -rf exe_create_parameter_mlgi_files*");

	if (iparams.hasBlockCleanUpDirectory){
		if (iparams.cleanupDirectory == 1){
			std::system ("rm -rf mesh_*.gmv mesh_*.lgi*");
			std::system ("rm -rf ellipse_cutoff_*.inp *.mlgi");
			std::system ("rm -rf intersections_*.inp tri_poly_fracture*.stor");
			//std::system("rm -rf ./output/meshes/ -rf ./output/polymeshes/");
			std::system ("rm -rf ./output/meshes/");
		}
		else if (iparams.cleanupDirectory == 2){
			std::system ("rm -rf mesh_*.gmv *.lgi");
			std::system ("rm -rf ellipse_cutoff_*.inp *.mlgi intersections_*.inp");
			std::system ("rm -rf intersections_*.inp tri_poly_fracture*.stor part_*.lg");
			std::system ("rm -rf ./output/meshes/ -rf ./output/polymeshes/");
		}
	}

	/*
	std::cout << " " << std::endl;
	cout << "total_intersection_number = " << total_intersection_number << endl;
	 */

	std::ifstream readinputfile("full_mesh.inp");
	std::string readinputline;
	std::getline(readinputfile,readinputline);
	std::istringstream ist(readinputline);
	int mesh_numvert, mesh_numelem;
	ist >> mesh_numvert >> mesh_numelem;

	std::ofstream outfile_info_mesh;
	outfile_info_mesh.open("output/info_mesh.txt", std::ios_base::app);
	outfile_info_mesh << mesh_numvert << std::endl;
	outfile_info_mesh << mesh_numelem;
	readinputfile.close();
	outfile_info_mesh.close();
	//std::cout << "\n" << std::endl;

	//printf("Finished in about %.0f seconds. \n", difftime(stop, start));

	if (runtime_status){
		std::cout << "Program run completed." << std::endl;
	}
	else{
		std::cout << "Program run terminated with errors." << std::endl;
	}

	return 0;
}
