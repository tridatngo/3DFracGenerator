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
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

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

// Define macro-variables
#include "include/macrovariables.hh"

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

// defaultUsageMessage
#include "include/defaultusagemessage.hh"

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
	const char* inputData = "ellipse_input.dat";

	int ncpu(1);
	bool getDataFileFromCmdL(false);

	for (int it = 0; it!= argc; it++ ){
		if (strcmp(argv[it], "-h") == 0 || strcmp(argv[it], "--help") == 0 ){
			std::cout << CGAL::defaultUsageMessage(std::string(argv[0])) << "\n";
			return 1;
		}
	}

	for (int it = 0; it!= argc; it++ ){
		//std::cout << "Arg " << it << " : " <<  argv[it] << std::endl;
		if (strcmp(argv[it], "-dataFile") == 0){
			getDataFileFromCmdL = true;
			inputData = argv[it+1];
		}
		if (strcmp(argv[it], "-np") == 0){
			ncpu = atoi(argv[it+1]);
		}
	}

	// Make output sub-directory
	//if ( boost::filesystem::exists( "output" )){
	if (stat("output", &info) == 0){
		//std::cout << "Folder output exists!" << std::endl;
		std::system("rm -rf output");
	}

	cmd = "rsync -r " + source_dir + "Executables/ Executables";
	std::system(cmd.c_str());
	std::system("rm -rf lagrit_logs* output*");
	std::system("rm -rf ellipse_cutoff* intersections* mesh_poly* parameters* tri_poly_fracture* *.zone");
	std::system("rm -rf mesh_poly* parameters* tri_poly_fracture* *.zone");

	// Remove pre-existing .lgi and .inp files
	std::system ("rm -rf full_mesh* *.lg *.lgi");
	std::system ("rm -rf *.inp *.mlgi");
	std::system ("rm -rf tri_fracture.stor lagrit_logs *boundary* ");
	std::system ("rm -rf logx3dgen outx3dge");

	// Create sub-directories
	std::system("mkdir output");
	std::system("mkdir output/ellipse_cutoff");
	std::system("mkdir output/ellipse_within_dfnbb");
	std::system("mkdir output/lagrit");
	std::system("mkdir output/parameters");
	std::system("mkdir lagrit_logs");

	//if ( boost::filesystem::exists( "intersections" )){
	if (stat("intersections", &info) == 0){	//std::cout << "Folder output exists!" << std::endl;
		std::system("rm -rf intersections/*.*");
		std::system("rm -rf intersections");
	}
	std::system("mkdir intersections");

	// [IMP] Be careful with the precision of number type!!!
	CGAL::Lazy_exact_nt<NT_MP>::set_relative_precision_of_to_double(1E-60);

	std::ifstream input(inputFile);
	Surface_mesh polygon1, polygon2, polygon3, current_poly, poly_empty;

	//Mesh mesh;
	int n4=160; // default value of n4

	InputParams iparams;
	//iparams = CGAL::parse_parameter_file<InputParams>(argc, argv);
	iparams = CGAL::parse_parameter_file<InputParams>(inputFile, argv);
	bool RM_CLOSE_PTS(false);
	bool MERGE_CLOSE_INTERSECTION_PTS(false);
	double mergeClosePointsRelCritLength(0.05);

	if (VERBOSE == 1){
		std::cout  << "iparams.hasn4 = " << iparams.hasn4 << std::endl;
		std::cout  << "iparams.n4 = " << iparams.n4 << std::endl;
		std::cout  << "iparams.DFNBB_LowerLeftX = " << iparams.DFNBB_LowerLeftX <<std::endl;
		std::cout  << "iparams.DFNBB_LowerLeftY = " << iparams.DFNBB_LowerLeftY <<std::endl;
		std::cout  << "iparams.DFNBB_LowerLeftZ = " << iparams.DFNBB_LowerLeftZ <<std::endl;
		std::cout  << "iparams.DFNBB_UpperRightX = " << iparams.DFNBB_UpperRightX <<std::endl;
		std::cout  << "iparams.DFNBB_UpperRightY = " << iparams.DFNBB_UpperRightY <<std::endl;
		std::cout  << "iparams.DFNBB_UpperRightZ = " << iparams.DFNBB_UpperRightZ <<std::endl;
		std::cout  << "iparams.addingPointsMethod = " << iparams.addingPointsMethod <<std::endl;
		std::cout  << "iparams.discreMethod = " << iparams.discreMethod <<std::endl;
		std::cout  << "iparams.cleanupDirectory = " << iparams.cleanupDirectory << std::endl;
		std::cout  << "iparams.hasdataFile = " << iparams.hasdataFile << std::endl;
		std::cout  << "iparams.dataFile = " << iparams.dataFile << std::endl;
		std::cout  << "iparams.rm_close_pts = " << iparams.rm_close_pts << std::endl;
	}

	if ( iparams.hasBlockRemoveClosePoints &&  iparams.hasremoveClosePoints ){
		RM_CLOSE_PTS = iparams.rm_close_pts;
	}

	if ( iparams.hasBlockRemoveClosePoints &&  iparams.hasmergeCloseIntersectionPoints ){
		MERGE_CLOSE_INTERSECTION_PTS = iparams.merge_close_intersection_pts;
	}
	else
	{
		std::cerr << "Warning: RemoveClosePoints block is set but the program failed to parse the mergeCloseIntersectionPoints value." << std::endl;
		std::cerr << "mergeCloseIntersectionPoints (default) = " << MERGE_CLOSE_INTERSECTION_PTS << std::endl;
	}

	if (MERGE_CLOSE_INTERSECTION_PTS){
		if (iparams.hasmergeClosePointsRelCritLength ){
			mergeClosePointsRelCritLength = iparams.mergeClosePointsRelCritLength;
		}
		else{
			std::cout << "Warning: mergeClosePointsRelCritLength value was not found in the input file. Default mergeClosePointsRelCritLength = " << mergeClosePointsRelCritLength
					  << " will be set." << std::endl;
		}
	}

	//std::cout << "** MERGE_CLOSE_INTERSECTION_PTS = " << MERGE_CLOSE_INTERSECTION_PTS << std::endl;
	//std::cout << "** mergeClosePointsRelCritLength = " << mergeClosePointsRelCritLength << std::endl;

	/*
	// Create a read-only string for the pointer iparams.dataFile
	std::string inputData_ = std::string(inputData);
	if (iparams.hasBlockDataFile && iparams.hasdataFile){
		inputData_ =  std::string(iparams.dataFile);
	}
	 */

	if ( iparams.hasn4 ){
		if(iparams.n4 == 0){
			std::cerr << "** Error : Density block is set but the program failed to parse n4 value."  << std::endl;
			std::cerr << "** Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else
		{
			n4 = iparams.n4;
		}
	}

	/*
	 * 	Choose the adding points method (addingPointsMethod)
	 *
	 *  1/ addingPointMethod = 1 : Along the intersection between two polygons P1 and P2, the intermediate points are added
	 *  with a density equal to 1/(min(dl1,dl2)), where dl1 dl2 denote the characteristic length of P1 and P2, respectively.
	 *
	 *	2/ addingPointMethod = 2: The intermediate points lie on only-P1-segment will be added with a density 1/dl1.
	 *	The intermediate points lie on only-P2-segment will be added with a density 1/dl1. The intermediate points lie on
	 *	common segments will be added with a density 1/(min(dl1,dl2)).
	 */

	int addingPointsMethod = 1; // default value

	// Read addingPointsMethod from input file
	if (iparams.hasaddingPointsMethod){
		if(iparams.addingPointsMethod == 0){
			std::cerr << "** Error: AddingPointsMethod block is set but the program failed to parse addingPointsMethod value."  << std::endl;
			std::cerr << "** The default value addingPointsMethod = 1 is used." << std::endl;
		}
		else
		{
			addingPointsMethod = iparams.addingPointsMethod;
		}
	}

	/*
	 * 	Choose the discretization method (discreMethod)
	 *
	 *  1/ discreMethod = 1 : Choose the number n4. Every ellipse will be approximated by a polygon consisting of n4 points.
	 *
	 *	2/ discreMethod = 2 : Choose the number n4. Compute valDens = distance(p_min,p_max) / n4. This is the characteristic
	 *	length for all polygon --> Homogeneous mesh.
	 */

	int discreMethod = 1; // default value

	// Read discreMethod from input file
	if (iparams.hasdiscreMethod){
		if(iparams.discreMethod == 0){
			std::cerr << "** Error: DiscreMethod block is set but the program failed to parse discreMethod value."  << std::endl;
			std::cerr << "** The default value addingPointsMethod = 1 is used." << std::endl;
		}
		else
		{
			discreMethod = iparams.discreMethod;
		}
	}

	double theta =2.*PI_/n4;

	//double target_edge_length = Lell/n4;
	// Create the DFN bounding box

	//Point_3 p_min(-4.0, -4.0, -4.0);
	//Point_3 p_max(2.0, 2.5, 4.0);
	//Point_3 p_min(-1.0, -0.5, -0.5);
	//Point_3 p_max(1.0, 0.5, 0.5);
	Point_3 p_min(iparams.DFNBB_LowerLeftX, iparams.DFNBB_LowerLeftY, iparams.DFNBB_LowerLeftZ);
	Point_3 p_max(iparams.DFNBB_UpperRightX, iparams.DFNBB_UpperRightY, iparams.DFNBB_UpperRightZ);

	//DFN_BB_ DFN_BB;

	DFN_BB_ DFN_BB = CGAL::create_dfn_bb<Surface_mesh, Point_3, DFN_BB_>(p_min,p_max);

	// print the bounding box to *.vtu
	if (PRINT_BBOX == 1){
		CGAL::print_dfn_bb<Surface_mesh, DFN_BB_> (DFN_BB);
	}


	/*========================================================================================================*/

	Ell_list_list	ell_child_list_list; // list of list of child ellipse
	Ell_list_in 	ell_list_in;
	Ell_list 		ell_list, ell_list_new, ell_child_list, ell_child_list_new, cEll_2_test_list, full_ell, full_ell_new;
	Ell_			ell_, ell_child, ell_empty;

	Poly_list_in	poly_list_in;
	//Output_IP 		output_ip;
	Output_IPL 		output_ipl;

	//std::string inputData (std::string("ellipse_input.dat"));

	std::cout << "Parameter file : " << inputFile << std::endl;

	//if (iparams.hasBlockDataFile && iparams.hasdataFile && boost::filesystem::exists( iparams.dataFile ) && argc <= 2){
	if (iparams.hasBlockDataFile && iparams.hasdataFile && stat(iparams.dataFile.c_str(), &info) == 0 && !getDataFileFromCmdL){
		std::cout << "Data file : " << iparams.dataFile << std::endl;

		//poly_list_in = CGAL::read_fracs_params_file_poly<Kernel, NT_MP, Poly_list_in> ( iparams.dataFile );

		if (iparams.hasdesactFracsFile){
			ell_list_in = CGAL::read_fracs_params_file_ell_desactfracs<Kernel, NT_MP, Ell_list_in> ( iparams.dataFile, iparams.desactFracsFile );
			poly_list_in = CGAL::read_fracs_params_file_poly_insidebbox_desactfracs<Kernel, NT_MP, Poly_list_in> ( iparams.dataFile, iparams.desactFracsFile, p_min, p_max );
		}
		else{
			ell_list_in = CGAL::read_fracs_params_file_ell<Kernel, NT_MP, Ell_list_in> ( iparams.dataFile );
			poly_list_in = CGAL::read_fracs_params_file_poly_insidebbox<Kernel, NT_MP, Poly_list_in> ( iparams.dataFile, p_min, p_max );
		}
	}
	else{
		if (iparams.hasdataFile){
			if (getDataFileFromCmdL){
				std::cout << "Get data from the command line." << std::endl;
			}
			else{
				std::cout << "Warning: The input data file \""<< iparams.dataFile <<
						"\" from parameter file does not exist. Choose the default file \"ellipse_input.dat\"." << std::endl;
			}
		}
		std::cout << "Data file : " << inputData << std::endl;

		//poly_list_in = CGAL::read_fracs_params_file_poly<Kernel, NT_MP, Poly_list_in> ( std::string(inputData) );
		if (iparams.hasdesactFracsFile){
			ell_list_in = CGAL::read_fracs_params_file_ell_desactfracs<Kernel, NT_MP, Ell_list_in> ( std::string(inputData), iparams.desactFracsFile );
			poly_list_in = CGAL::read_fracs_params_file_poly_insidebbox_desactfracs<Kernel, NT_MP, Poly_list_in> ( std::string(inputData), iparams.desactFracsFile, p_min, p_max );
		}
		else{
			ell_list_in = CGAL::read_fracs_params_file_ell<Kernel, NT_MP, Ell_list_in> ( std::string(inputData) );
			poly_list_in = CGAL::read_fracs_params_file_poly_insidebbox<Kernel, NT_MP, Poly_list_in> ( std::string(inputData), p_min, p_max );
		}
	}

	std::cout << "Number of ellipses from input file : " << ell_list_in.num_Ell << std::endl;
	std::cout << "Number of polygons from input file : " << poly_list_in.num_Poly << std::endl;

	//std::cout << "Number of lines from input file : " << CGAL::getnumline(iparams.dataFile) << std::endl;

	//exit(-1);

	// =====================================================================================================
	// Create surface mesh from data input (family 1: Ellipses)
	// =====================================================================================================
	if (discreMethod == 1){
		for (int it = 0; it !=ell_list_in.num_Ell; it++){
			Int_Iterator cE_name_in				=	ell_list_in.E_name.begin();		// Current ellipse
			Pts_Iterator cE_center_in			=	ell_list_in.E_center.begin();		// Center of the current ellipse
			NT_Iterator cE_xradius_in			=	ell_list_in.E_xradius.begin();		// current x_radius
			NT_Iterator cE_yradius_in			=	ell_list_in.E_yradius.begin();		// current y_radius
			Vector_Iterator cE_normalvect_in	=	ell_list_in.E_normalvect.begin(); 	// current normal vector

			std::advance(cE_name_in, it);
			std::advance(cE_center_in, it);
			std::advance(cE_xradius_in, it);
			std::advance(cE_yradius_in, it);
			std::advance(cE_normalvect_in, it);

			NT_MP Ex(*cE_xradius_in), Ey(*cE_yradius_in);	// x_radius, y_radius
			Vector_3 normal_vect = *cE_normalvect_in;				// normal vector
			Vector_3 trans_vect = Vector_3((*cE_center_in).x(), (*cE_center_in).y(), (*cE_center_in).z()); // translation vector

			for (int i = 0; i<n4; i++){
				Point_3 point_to_add(Point_3(Ex*cos(i*theta),Ey*sin(i*theta),0));
				//Point_3 point_to_add_rot = CGAL::PointRotation<Kernel, NT_MP>(point_to_add, normal_vect, trans_vect);
				Point_3 point_to_add_rot = CGAL::PointInverseRotation<Kernel, NT_MP>(point_to_add, normal_vect, trans_vect);
				current_poly.add_vertex(point_to_add_rot);
			}

			double val1 = 2.*(pow(CGAL::to_double(Ex),2.0)+ pow(CGAL::to_double(Ey),2.0));
			double val2 = 0.5*pow(CGAL::to_double(Ex)-CGAL::to_double(Ey),2.0);
			double Lell = PI_* sqrt(val1 - val2);
			double target_edge_length = Lell/n4;

			Surface_mesh::Vertex_range current_m_range = current_poly.vertices();
			current_poly.add_face(current_m_range);

			std::ostringstream cE_name_in_ss;
			cE_name_in_ss << *cE_name_in;

			//fname = "output/ellipse_" + cE_name_in_ss.str();
			//CGAL::output_off_vtu <Surface_mesh>(fname, current_poly, false, true);

			/*
			ell_list.E_name.push_back(*cE_name_in);
			ell_list.E_mesh.push_back(current_poly);
			ell_list.E_normalvect.push_back(*cE_normalvect_in);
			 */
			ell_.E_name					=	*cE_name_in;			// Name of the ellipse
			ell_.E_normalvect			=	*cE_normalvect_in;		// Normal vector of the ellipse
			ell_.E_mesh					=	current_poly;			// Mesh of the ellipse
			ell_.E_target_edge_length 	= 	target_edge_length;
			ell_.Parent_E_name			=	*cE_name_in;			// Name of the ellipse

			ell_list.push_back(ell_);
			current_poly 	= poly_empty; 	// clean up the polygon mesh
			ell_			= ell_empty;	// clean up the polygon
		}
	}
	else if (discreMethod ==2){

		for (int it = 0; it !=ell_list_in.num_Ell; it++){

			Int_Iterator cE_name_in				=	ell_list_in.E_name.begin();		// Current ellipse
			Pts_Iterator cE_center_in			=	ell_list_in.E_center.begin();		// Center of the current ellipse
			NT_Iterator cE_xradius_in			=	ell_list_in.E_xradius.begin();		// current x_radius
			NT_Iterator cE_yradius_in			=	ell_list_in.E_yradius.begin();		// current y_radius
			Vector_Iterator cE_normalvect_in	=	ell_list_in.E_normalvect.begin(); 	// current normal vector

			std::advance(cE_name_in, it);
			std::advance(cE_center_in, it);
			std::advance(cE_xradius_in, it);
			std::advance(cE_yradius_in, it);
			std::advance(cE_normalvect_in, it);

			NT_MP Ex(*cE_xradius_in), Ey(*cE_yradius_in);	// x_radius, y_radius
			Vector_3 normal_vect = *cE_normalvect_in;				// normal vector
			Vector_3 trans_vect = Vector_3((*cE_center_in).x(), (*cE_center_in).y(), (*cE_center_in).z()); // translation vector

			double target_edge_length = CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(p_min, p_max))) / n4;

			double val1 = 2.*(pow(CGAL::to_double(Ex),2.0)+ pow(CGAL::to_double(Ey),2.0));
			double val2 = 0.5*pow(CGAL::to_double(Ex)-CGAL::to_double(Ey),2.0);
			double Lell = PI_* sqrt(val1 - val2);
			double numPoints_ = round(Lell/target_edge_length);
			double theta_ = 2.*PI_/numPoints_;

			for (int i = 0; i<numPoints_; i++){
				Point_3 point_to_add(Point_3(Ex*cos(i*theta_),Ey*sin(i*theta_),0));
				//Point_3 point_to_add_rot = CGAL::PointRotation<Kernel, NT_MP>(point_to_add, normal_vect, trans_vect);
				Point_3 point_to_add_rot = CGAL::PointInverseRotation<Kernel, NT_MP>(point_to_add, normal_vect, trans_vect);
				current_poly.add_vertex(point_to_add_rot);
			}


			Surface_mesh::Vertex_range current_m_range = current_poly.vertices();
			current_poly.add_face(current_m_range);

			std::ostringstream cE_name_in_ss;
			cE_name_in_ss << *cE_name_in;

			//fname = "output/ellipse_" + cE_name_in_ss.str();
			//CGAL::output_off_vtu <Surface_mesh>(fname, current_poly, false, true);

			/*
			ell_list.E_name.push_back(*cE_name_in);
			ell_list.E_mesh.push_back(current_poly);
			ell_list.E_normalvect.push_back(*cE_normalvect_in);
			 */
			ell_.E_name					=	*cE_name_in;			// Name of the ellipse
			ell_.E_normalvect			=	*cE_normalvect_in;		// Normal vector of the ellipse
			ell_.E_mesh					=	current_poly;			// Mesh of the ellipse
			ell_.E_target_edge_length 	= 	target_edge_length;
			ell_.Parent_E_name			=	*cE_name_in;			// Name of the ellipse

			ell_list.push_back(ell_);
			current_poly 	= poly_empty; 	// clean up the polygon mesh
			ell_			= ell_empty;	// clean up the polygon
		}
	}

	// =====================================================================================================
	// Create surface mesh from data input (family 2: Polygons)
	// =====================================================================================================

	for (int it = 0; it != poly_list_in.num_Poly; it++){
		Int_Iterator cP_name_in				=	poly_list_in.P_name.begin();		// Current polygon
		Vector_Iterator cP_normalvect_in	=	poly_list_in.P_normalvect.begin(); 	// current normal vector
		List_Pts_Iterator cP_list_points_in	=	poly_list_in.P_points.begin(); 	// current normal vector
		//Pts_Iterator cP_points_in				=	poly_list_in.P_normalvect.begin(); 	// current normal vector

		std::advance(cP_name_in, it);
		std::advance(cP_normalvect_in, it);
		std::advance(cP_list_points_in, it);

		Pts_Iterator first_point_it =  (*cP_list_points_in).begin();

		Vector_3 normal_vect = *cP_normalvect_in;				// normal vector
		Vector_3 trans_vect = Vector_3((*first_point_it).x(), (*first_point_it).y(), (*first_point_it).z()); // translation vector

		double target_edge_length;
		if (discreMethod == 1){
			target_edge_length = CGAL::get_poly_perimeter<Kernel>(*cP_list_points_in)/n4;
		}
		else if (discreMethod == 2){
			//double target_edge_length = Lell/n4;
			target_edge_length = CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(p_min, p_max))) / n4;
		}

		//for (int i = 0; i<n4; i++){
		// Add new points
		std::list<Point_3> cP_list_points_update;
		cP_list_points_update.push_back(*(*cP_list_points_in).begin());
		unsigned int ip_ = 0;
		for( ; ip_ < (*cP_list_points_in).size(); ++ip_) {
			Pts_Iterator it_p1, it_p2;
			if (ip_ == (*cP_list_points_in).size()-1){
				it_p1 = (*cP_list_points_in).begin();
				std::advance(it_p1, ip_);

				it_p2 = (*cP_list_points_in).begin();
			}
			else{
				it_p1 = (*cP_list_points_in).begin();
				std::advance(it_p1, ip_);
				it_p2 = (*cP_list_points_in).begin();
				std::advance(it_p2, ip_+1);
			}

			//int numAddedPoints = round(CGAL::abs(*it_p2 - *it_p1)/ target_edge_length);
			double distance_ =  CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(*it_p2, *it_p1)));
			int numAddedPoints(0);
			//numAddedPoints = ceil(distance_ / target_edge_length)-1;
			numAddedPoints = round(distance_ / target_edge_length);

			if (VERBOSE == 1){
				std::cout << "*(*cP_list_points_in).begin() = " << *(*cP_list_points_in).begin() << std::endl;
				std::cout << "Point 1 = " << *it_p1 << std::endl;
				std::cout << "Point 2 = " << *it_p2 << std::endl;
				std::cout << "Number of added points (poly_1) = " << numAddedPoints << std::endl;
			}

			for (int i = 0; i!= (numAddedPoints+1); i++)
				cP_list_points_update.push_back(*it_p1 + (i+1) * (*it_p2 - *it_p1)/(numAddedPoints+1));

		}
		/*
			for (int i = 0; i<(*cP_list_points_in).size(); i++){
				Pts_Iterator point_it =  (*cP_list_points_in).begin();
				std::advance(point_it, i);

				Point_3 point_to_add(*point_it);
				current_poly.add_vertex(point_to_add);
			}
		 */

		for (int i = 0; i<cP_list_points_update.size(); i++){
			Pts_Iterator point_it = cP_list_points_update.begin();
			std::advance(point_it, i);

			Point_3 point_to_add(*point_it);
			current_poly.add_vertex(point_to_add);
		}

		Surface_mesh::Vertex_range current_m_range = current_poly.vertices();
		current_poly.add_face(current_m_range);

		std::ostringstream cP_name_in_ss;
		cP_name_in_ss << *cP_name_in;

		//std::cout << "*cP_normalvect_in = " <<  *cP_normalvect_in << std::endl;

		ell_.E_name					=	*cP_name_in;			// Name of the ellipse
		ell_.E_normalvect			=	*cP_normalvect_in;		// Normal vector of the ellipse
		ell_.E_mesh					=	current_poly;			// Mesh of the ellipse
		ell_.E_target_edge_length 	= 	target_edge_length;
		ell_.Parent_E_name			=	*cP_name_in;			// Name of the ellipse

		ell_list.push_back(ell_);
		current_poly 	= poly_empty; 	// clean up the polygon mesh
		ell_			= ell_empty;	// clean up the polygon
	}

	std::cout << "ell_list.size() = " <<  ell_list.size() << std::endl;

	//exit(1);

    time(&stop_input);
	clock_t end_input = clock();
	double CPU_input_secs = double(end_input - begin) / CLOCKS_PER_SEC;

	// =====================================================================================================
	// Test ellipse
	// =====================================================================================================
	if(DEBUG_NTD){
		fname = "test_ellipse_1";
		for (int it = 1; it !=2; it++){
			int 			cE_name;					// Name of the current ellipse
			Vector_3 		cE_normalvect;				// Normal vector of the current ellipse
			Surface_mesh 	cE_mesh;                 	// Mesh of the current ellipse
			Ell_			cEll_, Ell_empty;			// Current ellipse
			Ell_Iterator	cEll_iter = ell_list.begin();
			std::advance(cEll_iter, it);
			cEll_ = *cEll_iter;
			std::cout << "cEll_.E_normalvect = " << cEll_.E_normalvect << std::endl;
			CGAL::output_off_vtu <Surface_mesh>(fname, cEll_.E_mesh, false, true);


			unsigned int i(0), end( cEll_.E_mesh.number_of_vertices());

			std::fstream fs, fsx, fsy, fsz;
			fsx.open ("rennesX.dfn",  std::fstream::in | std::fstream::out | std::fstream::app);
			for( ; i < end; ++i) {
				vertex_descriptor vb_test(i);
				//std::cout << "cEll_.E_mesh.point(vb_test).x() = " << cEll_.E_mesh.point(vb_test).x() << std::endl;
				fsx << cEll_.E_mesh.point(vb_test).x() << std::endl;
			}
			fsx.close();

			fsy.open ("rennesY.dfn",  std::fstream::in | std::fstream::out | std::fstream::app);
			i = 0;
			for( ; i < end; ++i) {
				vertex_descriptor vb_test(i);
				//std::cout << "cEll_.E_mesh.point(vb_test).y() = " << cEll_.E_mesh.point(vb_test).y() << std::endl;
				fsy << cEll_.E_mesh.point(vb_test).y() << std::endl;
			}
			fsy.close();


			fsz.open ("rennesZ.dfn",  std::fstream::in | std::fstream::out | std::fstream::app);
			i = 0;
			for( ; i < end; ++i) {
				vertex_descriptor vb_test(i);
				//std::cout << "cEll_.E_mesh.point(vb_test).z() = " << cEll_.E_mesh.point(vb_test).z() << std::endl;
				fsz << cEll_.E_mesh.point(vb_test).z() << std::endl;
			}
			fsz.close();

			fs.open ("rennesnbSommets.dfn", std::fstream::out);
			fs << cEll_.E_mesh.number_of_vertices();
			fs.close();
		}

		// If you want to stop here
		exit(-1);
	}

	// =====================================================================================================
	//                              Find out outside-dfn-bounding-box polygons
	// =====================================================================================================
	std::list<int> should_be_removed, should_be_removed_new;

	for (int it = 0; it !=ell_list.size(); it++){
		int 			cE_name;					// Name of the current ellipse
		Vector_3 		cE_normalvect;				// Normal vector of the current ellipse
		Surface_mesh 	cE_mesh;                 	// Mesh of the current ellipse
		Ell_			cEll_, Ell_empty;			// Current ellipse
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		//std::cout << "cEll_.normal_vector =" << cEll_.E_normalvect << std::endl;

		Iso_cuboid_3 isoC = CGAL::My_Bounding_Box<Kernel, NT_MP>(cEll_.E_mesh);

		if (VERBOSE == 1){
			std::cout << "xmin of the BB of ellipse " << cEll_.E_name << " is " << isoC.xmin() << std::endl;
			std::cout << "xmax of the BB of ellipse " << cEll_.E_name << " is " << isoC.xmax() << std::endl;
			std::cout << "ymin of the BB of ellipse " << cEll_.E_name << " is " << isoC.ymin() << std::endl;
			std::cout << "ymax of the BB of ellipse " << cEll_.E_name << " is " << isoC.ymax() << std::endl;
			std::cout << "zmin of the BB of ellipse " << cEll_.E_name << " is " << isoC.zmin() << std::endl;
			std::cout << "zmax of the BB of ellipse " << cEll_.E_name << " is " << isoC.zmax() << std::endl;

			std::cout << "xmin of the DFN_BB is " << p_min[0] << std::endl;
			std::cout << "xmax of the DFN_BB is " << p_max[0] << std::endl;
			std::cout << "ymin of the DFN_BB is " << p_min[1] << std::endl;
			std::cout << "ymax of the DFN_BB is " << p_max[1] << std::endl;
			std::cout << "zmin of the DFN_BB is " << p_min[2] << std::endl;
			std::cout << "zmax of the DFN_BB is " << p_max[2] << std::endl;
		}

		//if ( isoC.xmin() >= p_max[0] || isoC.xmax() <= p_min[0] || isoC.ymin() >= p_max[1] || isoC.ymax() <= p_min[1] || isoC.zmin() >= p_max[2] || isoC.zmax() <= p_min[2]){
		if ( 	isoC.xmin() > p_max[0] || (isoC.xmin() == p_max[0] && isoC.xmin() < isoC.xmax()) ||
				isoC.xmax() < p_min[0] || (isoC.xmax() == p_min[0] && isoC.xmax() > isoC.xmin()) ||
				isoC.ymin() > p_max[1] || (isoC.ymin() == p_max[1] && isoC.ymin() < isoC.ymax()) ||
				isoC.ymax() < p_min[1] || (isoC.ymax() == p_min[1] && isoC.ymax() > isoC.ymin()) ||
				isoC.zmin() > p_max[2] || (isoC.zmin() == p_max[2] && isoC.zmin() < isoC.zmax()) ||
				isoC.zmax() < p_min[2] || (isoC.zmax() == p_min[2] && isoC.zmax() > isoC.zmin())){
			std::cout << "The polygon " << cEll_.E_name << " is not inside the DFN_BB and should be removed." << std::endl;
			should_be_removed.push_back(it);
		}
	}

	for (Int_Iterator i = should_be_removed.begin(); i != should_be_removed.end(); i++)
	{
		std::cout << "Polygon " << *i << " will be removed." << std::endl;
	}

	// update ell_list : Remove disconnected polygons / polygons which do not lie in DFN_BB
	ell_list = CGAL::update_polygon_numbering <Kernel, Ell_>(ell_list);
	ell_list = CGAL::remove_polygons < Kernel, Ell_> (ell_list, should_be_removed);
	//std::cout << "Number of ellipses after removing outside-dfn-bounding-box polygons = " <<  ell_list.size() << std::endl;

	if (ell_list.size() == 0){
		std::cout << "No fracture inside the DFN bounding box. Program terminated." << std::endl;
		exit(-1);
	}
	//std::cout << std::endl << "** should_be_removed.size() = " << should_be_removed.size() << std::endl;

	// Remove the outside-bounding-box part of the polygons
	for (int it = 0; it !=ell_list.size(); it++){
		Ell_			cEll_, Ell_empty;			// Current ellipse
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;

		if (!RM_CLOSE_PTS){
			cEll_.E_mesh = CGAL::remove_out_of_bb_parts <Kernel, Output_Child_Polys,IP_Ell, Output_IP, DFN_BB_>
			(DFN_BB, cEll_.E_mesh, cEll_.E_target_edge_length);
		}

		else{
			NT_MP tol_length = cEll_.E_target_edge_length * 0.1; // default tol_length

			if (iparams.hasminLengthRatio){
				//std::cout << "iparams.minLengthRatio = " << iparams.minLengthRatio << std::endl;
				tol_length = cEll_.E_target_edge_length * iparams.minLengthRatio;
			}
			cEll_.E_mesh = (CGAL::remove_out_of_bb_parts_rm_pts <Kernel, Output_Inside_BB, Output_Child_Polys, IP_Ell, Output_IP, DFN_BB_, NT_MP>
			(DFN_BB, cEll_.E_mesh, cEll_.E_target_edge_length, tol_length)).poly_inside_bb;

			if (!(CGAL::remove_out_of_bb_parts_rm_pts <Kernel, Output_Inside_BB, Output_Child_Polys, IP_Ell, Output_IP, DFN_BB_, NT_MP>
			(DFN_BB, cEll_.E_mesh, cEll_.E_target_edge_length, tol_length)).inside_bb){
				std::cout << "Mark the fracture as should_be_removed_new " << std::endl;
				should_be_removed_new.push_back(it);
			}
		}

		if (0){
			std::ostringstream cE_name_ss;
			cE_name_ss << cEll_.E_name;
			fname = "output/ellipse_within_dfnbb/ellipse_within_dfnbb_" + cE_name_ss.str();
			CGAL::output_off_vtu <Surface_mesh>(fname, cEll_.E_mesh, false, true);
		}

		ell_list_new.push_back(cEll_);
	}

	//std::cout << "** should_be_removed_new.size() = " << should_be_removed_new.size() << std::endl;

	ell_list = ell_list_new;
	ell_list_new.clear();

	// update ell_list : Remove disconnected polygons / polygons which do not lie in DFN_BB
	ell_list = CGAL::update_polygon_numbering <Kernel, Ell_>(ell_list);
	ell_list = CGAL::remove_polygons < Kernel, Ell_> (ell_list, should_be_removed_new); // NTD 09Feb2017 desactFracsFile.dat discrepancies
	std::cout << "Number of ellipses after removing outside-dfn-bounding-box polygons = " <<  ell_list.size() << std::endl;

	if (ell_list.size() == 0){
		std::cout << "No fracture inside the DFN bounding box. Program terminated." << std::endl;
		exit(-1);
	}

	for (int it = 0; it !=ell_list.size(); it++){
		int 			cE_name;					// Name of the current ellipse
		Vector_3 		cE_normalvect;				// Normal vector of the current ellipse
		Surface_mesh 	cE_mesh;                 	// Mesh of the current ellipse
		Ell_			cEll_, Ell_empty;			// Current ellipse
		Ell_Iterator	cEll_iter = ell_list.begin();
		std::advance(cEll_iter, it);
		cEll_ = *cEll_iter;


		// Export the list of fractures inside the DFN bounding box
		std::ofstream outfile_insidebbox;
		outfile_insidebbox.open("output/insidebbox_frac.txt", std::ios_base::app);
		//std::cout << "Polygon parent name " << cEll_.Parent_E_name << std::endl;
		outfile_insidebbox << cEll_.Parent_E_name << std::endl;
		outfile_insidebbox.close();

		// Export the fractures inside the DFN bounding box to *.vtu file to visualise
		std::ostringstream cE_name_ss;
		cE_name_ss << cEll_.E_name;
		fname = "output/ellipse_within_dfnbb/ellipse_within_dfnbb_" + cE_name_ss.str();
		//CGAL::output_off_vtu <Surface_mesh>(fname, cEll_.E_mesh, false, true);
		CGAL::output_ell_to_vtu <Kernel,NT_MP>(fname, cEll_.E_mesh);
	}

	// Update parent polygon numbering
	ell_list = CGAL::update_parent_polygon_numbering <Kernel, Ell_>(ell_list);

	if (VERBOSE == 1){
		for (int i = 0; i !=ell_list.size(); i++){
			Ell_			cEll_1;			// Current ellipse
			Ell_Iterator	cEll_1_iter = ell_list.begin();
			std::advance(cEll_1_iter, i);
			cEll_1 = *cEll_1_iter;
			std::cout << "Polygon name " << cEll_1.E_name << std::endl;
			std::cout << "Polygon parent name " << cEll_1.Parent_E_name << std::endl;
		}
	}

	//exit(1);

	// Write polys.inp file
	CGAL::write_polys_inp_file<Kernel, Ell_, NT_MP>(ell_list, p_min, p_max, std::string("output/polys.inp"));

	// Write params_DFNWorks.txt file
	CGAL::write_params_DFNWorks_txt_file<Kernel, Ell_, NT_MP>(ell_list, p_min, p_max, std::string("params_DFNWorks.txt"));


    time(&stop_check_out_of_dfnbb);
	clock_t end_check_out_of_dfnbb = clock();
	double CPU_out_of_dfnbb_secs = double(end_check_out_of_dfnbb - end_input) / CLOCKS_PER_SEC;

	// =====================================================================================================
	//                                     Quick-checking intersection between polygons
	// =====================================================================================================

	int total_intersection_number(0);
	for (int i = 0; i !=ell_list.size(); i++){
		//for (int i = 0; i !=1; i++){
		Ell_			cEll_1, cEll_2, cEll_child, Ell_empty;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;

		std::cout << "cEll_1.Parent_E_name: " << cEll_1.Parent_E_name  << std::endl;

		CGAL::silent_line('=',80);
		//std::cout << "Find the intersection of the polygon " << cEll_1.E_name << std::endl;

		std::list<int> 		checkintersection;

		std::list<Line_3> 			listIntLine;
		std::list<Line_3_epeck> 	listIntLine_epeck;
		std::list<int> 		listIntEll;

		for (int j = 0; j !=ell_list.size(); j++){

			/*
			 * @todo NTD Optimize the algorithm by testing only the intersection
			 * between polygon i and polygon j with j > i. In the inverse case,
			 * result should be taken from the reciprocal test.
			 */

			if (j == i){
				checkintersection.push_back(0);
			}

			if (j != i){
				Ell_Iterator 	cEll_2_iter = ell_list.begin();
				std::advance(cEll_2_iter, j);
				cEll_2 = *cEll_2_iter;
				if (VERBOSE == 1){
					std::cout << "Test intersection of the polygon " << cEll_1.E_name
							<< " with the polygon "<< cEll_2.E_name <<std::endl;
				}
				// 05/01/2017 Octree approach for intersection determination between two ellipses

				Iso_cuboid_3 isoE1 = CGAL::My_Bounding_Box<Kernel, NT_MP>(cEll_1.E_mesh);
				Iso_cuboid_3 isoE2 = CGAL::My_Bounding_Box<Kernel, NT_MP>(cEll_2.E_mesh);

				if (VERBOSE == 1){
					std::cout << "xmin of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.xmin() << std::endl;
					std::cout << "xmax of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.xmax() << std::endl;
					std::cout << "ymin of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.ymin() << std::endl;
					std::cout << "ymax of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.ymax() << std::endl;
					std::cout << "zmin of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.zmin() << std::endl;
					std::cout << "zmax of the BB of ellipse " << cEll_1.E_name << " is " << isoE1.zmax() << std::endl;

					std::cout << "xmin of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.xmin() << std::endl;
					std::cout << "xmax of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.xmax() << std::endl;
					std::cout << "ymin of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.ymin() << std::endl;
					std::cout << "ymax of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.ymax() << std::endl;
					std::cout << "zmin of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.zmin() << std::endl;
					std::cout << "zmax of the BB of ellipse " << cEll_2.E_name << " is " << isoE2.zmax() << std::endl;
				}

				if ( 	isoE1.xmin() > isoE2.xmax() ||	isoE1.xmax() < isoE2.xmin() ||
						isoE1.ymin() > isoE2.ymax() ||	isoE1.ymax() < isoE2.ymin() ||
						isoE1.zmin() > isoE2.zmax() ||	isoE1.zmax() < isoE2.zmin()){
					checkintersection.push_back(0);
				}
				else{
					// NTD 21/09/2016: Imprecision of intersection between line consisting of close-to-zero point
					output_ipl = CGAL::intersections_of_two_polys_exp_line<Kernel,Output_IPL> (cEll_1.E_mesh, cEll_2.E_mesh);
					//output_ipl = CGAL::intersections_of_two_polys_exp_line_MP<NT_MP, Surface_mesh,Output_IPL> (cEll_1.E_mesh, cEll_2.E_mesh);

					if (output_ipl.num_intersec_1_ + output_ipl.num_intersec_2_ == 3){
						/*
						 * @todo NTD - two tangent polygons
						 */
						std::cout << "Two polygons " << cEll_1.E_name << " and "
								<< cEll_2.E_name << " are tangent." <<std::endl;
						checkintersection.push_back(1);
					}
					else{
						if(output_ipl.num_intersec_1_ + output_ipl.num_intersec_2_ == 4)
						{
							std::cout << "The polygon " << cEll_1.E_name << " intersects the polygon "
									<< cEll_2.E_name << std::endl;
							checkintersection.push_back(2);

						}
						else{
							checkintersection.push_back(0);
						}
					}

					if (output_ipl.bool_intersecting_){
						if (VERBOSE == 1){
							std::cout << "[Main program] Found intersection between the plane " << cEll_1.E_name << " and the plane " << cEll_2.E_name << std::endl;
							std::cout << std::setprecision(30) << output_ipl.intersection_line_ << std::endl;
						}

						//listIntLine.push_back(output_ipl.intersection_line_);

						listIntLine.push_back(output_ipl.intersection_line_);
						listIntLine_epeck.push_back(output_ipl.intersection_line_epeck_);

						listIntEll.push_back(cEll_2.E_name);
					}
					/*
				else{
					listIntLine.push_back( Line_3( Point_3(0, 0, 0) , Point_3( 0, 0, 0) ));
					listIntEll.push_back(77777777);
				}
					 */
				}

			}
		}

		fname = "intersections/intersection_orig_" + CGAL::int2str(cEll_1.E_name) + ".dat";
		std::ofstream    fname_ofs;
		//fname_ofs.open(fname.c_str());
		fname_ofs.open(fname.c_str());

		for (Int_Iterator i = checkintersection.begin(); i != checkintersection.end(); i++){
			fname_ofs <<  *i << "\n";
			if(*i > 0){
				total_intersection_number++;
			}
		}

		fname = "intersections/intersection_line_" + CGAL::int2str(cEll_1.E_name) + ".dat";
		std::ofstream    fnamel_ofs;
		fnamel_ofs.open(fname.c_str());

		int inum = 0;
		if(!USE_CONVERT_TO_EPECK){
			for (Int_Iterator it = listIntEll.begin(); it != listIntEll.end(); it++){
				Line_Iterator listIntEll_iter(listIntLine.begin());
				std::advance(listIntEll_iter,inum);
				inum++;

				fnamel_ofs << *it << " " << std::setprecision(30)  << std::uppercase << std::scientific << *listIntEll_iter << "\n";
			}
		}
		else{
			for (Int_Iterator it = listIntEll.begin(); it != listIntEll.end(); it++){
				Line_Iterator_epeck listIntEll_epeck_iter(listIntLine_epeck.begin());
				std::advance(listIntEll_epeck_iter,inum);
				inum++;

				fnamel_ofs << *it << " " << std::setprecision(30)  << std::uppercase << std::scientific << *listIntEll_epeck_iter << "\n";
			}
		}

		checkintersection.clear();
		fname_ofs.close();
	}

	total_intersection_number = total_intersection_number /2;

	std::ofstream    fname_total_inum_ofs("output/info_mesh.txt");
	fname_total_inum_ofs << "# ell_list.size() --> total_intersection_number --> number_of_vertices --> number_of_elements" << std::endl;
	fname_total_inum_ofs << ell_list.size() << std::endl;
	fname_total_inum_ofs << total_intersection_number << std::endl;
	fname_total_inum_ofs.close();

	fname = "intersections/intersections.inp";
	std::ofstream    fname_il_ofs;
	fname_il_ofs.open(fname.c_str());
	fname_il_ofs << "x, y, z, ID of fracture1, ID of fracture 2" << "\n";

	/*
	for (int i = 0; i !=1; i++){
		Ell_			cEll_1, cEll_2, Ell_empty;			// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;

		CGAL::silent_line('=',80);
		//std::cout << "Find the intersection of the polygon " << cEll_1.E_name << std::endl;

	}
	 */

	time(&stop_check_intersection);
	clock_t end_check_intersection = clock();
	double CPU_intersection_secs = double(end_check_intersection - end_check_out_of_dfnbb) / CLOCKS_PER_SEC;

	// =====================================================================================================
	//                                     Loop for cutting off polygons
	// =====================================================================================================

	for (int i = 0; i !=ell_list.size(); i++){
		//for (int i = 0; i !=1; i++){
		CGAL::silent_line('=',80);

		Ell_			cEll_1, cEll_2, cEll_2_test, cEll_child, cEll_child_pre, cEll_child_re, Ell_empty;	// Current ellipse
		Ell_Iterator	cEll_1_iter = ell_list.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;
		Output_Child_Polys child_polys;
		ell_child_list.push_back(cEll_1);

		fname = "intersections/intersection_orig_" + CGAL::int2str(cEll_1.E_name) + ".dat";
		std::ifstream    infile(fname.c_str());

		if (! infile.is_open()) {
			std::cerr << "Failed to open the input file " << fname << "." << std::endl;
			return -1;
		}

		std::list<int> checkintersection_;
		int intersection_;

		for (int it = 0; it != ell_list.size(); it++){
			std::string line;
			std::getline(infile, line);
			int dummy;
			std::istringstream ist(line);
			ist.clear();
			ist >> dummy;
			checkintersection_.push_back(dummy);
		}

		for (int j = 0; j != ell_list.size(); j++){
			std::cout << "Finding intersection between the ellipses " << i << " and " << j << std::endl;

			CGAL::silent_line('-',60);

			Int_Iterator  intersection_ = checkintersection_.begin();
			std::advance(intersection_, j);

			Ell_Iterator cEll_2_iter = ell_list.begin();
			std::advance(cEll_2_iter, j);
			cEll_2 = *cEll_2_iter;

			std::cout << "i = "<< i <<", j = "<< j <<", *intersection_ = " << *intersection_
					<< ", ell_child_list.size() = " << ell_child_list.size()<< std::endl;

			//if (*intersection_ > 0 && j > i){
			if (*intersection_ > 0){

				// We have to check intersection between ellipse i and the child-ellipse j_n
				int k(0);
				while (k < 1E+20){
					std::cout << "k = " << k << std::endl;
					//std::cout << "ell_child_list.size() = " << ell_child_list.size() << std::endl;

					if ( k  > ell_child_list.size() -1){
						break;
					}

					Ell_Iterator cEll_child_iter = ell_child_list.begin();
					std::advance(cEll_child_iter, k);
					cEll_child = *cEll_child_iter;

					/*
						if (k >0){
							for (int l = 0; l < (k-1); l++ ){
								std::cout << "k = " << k << ", l = " << l << std::endl;
								ell_child_list_new.push_back(cEll_child);
							}
						}
					 */

					for (int l = 0; l < k; l++ ){
						std::cout << "Push back the previous ellipses : k = " << k << ", l = " << l << ", ell_child_list.size() = " <<ell_child_list.size() << std::endl;
						Ell_Iterator cEll_child_iter_pre = ell_child_list.begin();
						std::advance(cEll_child_iter_pre, l);
						cEll_child_pre = *cEll_child_iter_pre;
						ell_child_list_new.push_back(cEll_child_pre);
					}
					//double target_edge_length_1 = target_edge_length * 0.5;

					//child_polys = CGAL::cutting_off_polygons <Kernel, Output_Child_Polys, IP_Ell, Output_IP>
					//( cEll_child.E_mesh, cEll_2.E_mesh, target_edge_length * 0.5);
					if (VERBOSE == 1){
						std::cout << "cEll_child.E_mesh.num_vertices() = " << cEll_child.E_mesh.num_vertices() << std::endl;
						std::cout << "cEll_2.E_mesh.num_vertices() = " << cEll_2.E_mesh.num_vertices() << std::endl;
					}
					// Determine the target edge length to add points
					double target_edge_length = CGAL::min(cEll_child.E_target_edge_length, cEll_2.E_target_edge_length);
					if (VERBOSE == 1){
						std::cout << "cEll_child.E_target_edge_length = " << cEll_child.E_target_edge_length << std::endl;
						std::cout << "cEll_2.E_target_edge_length = " << cEll_2.E_target_edge_length << std::endl;
					}

					if (addingPointsMethod ==1){
						if (!RM_CLOSE_PTS){
							child_polys = CGAL::cutting_off_first_polygon <Kernel, Output_Child_Polys, Ell_, IP_Ell, Output_IPL, NT_MP>
							(ell_list, cEll_child, cEll_2, target_edge_length);
						}
						else{
							NT_MP tol_length = target_edge_length * 0.1;

							if (iparams.hasminLengthRatio){
								//std::cout << "iparams.minLengthRatio = " << iparams.minLengthRatio << std::endl;
								tol_length = target_edge_length * iparams.minLengthRatio;
							}

							child_polys = CGAL::cutting_off_first_polygon_rmpts <Kernel, Output_Child_Polys, Ell_, IP_Ell, Output_IPL, NT_MP>
							(ell_list, cEll_child, cEll_2, target_edge_length, tol_length, MERGE_CLOSE_INTERSECTION_PTS, mergeClosePointsRelCritLength);
						}
					}
					else if  (addingPointsMethod ==2){
						child_polys = CGAL::cutting_off_two_polygons <Kernel, Output_Child_Polys, IP_Ell, Output_IP>
						( cEll_child.E_mesh, cEll_2.E_mesh, cEll_child.E_target_edge_length, cEll_2.E_target_edge_length);
					}

					if (VERBOSE == 1){
						std::cout << "ell_child_list_new.size() = " << ell_child_list_new.size() << std::endl;
						std::cout << "ell_child_list.size() = " << ell_child_list.size() << std::endl;
					}
					// std::cout << "child_polys.no_polygons_interst_ = " << child_polys.no_polygons_interst_ << std::endl;

					if (child_polys.no_polygons_interst_ == false){
						// Add the child polygons to the new list
						ell_child.E_name				=	cEll_1.E_name;					// Name of the ellipse
						ell_child.E_target_edge_length	=	cEll_1.E_target_edge_length;	// target edge length of the ellipse
						ell_child.E_normalvect			=	cEll_1.E_normalvect;			// Normal vector of the ellipse
						ell_child.E_mesh				=	child_polys.poly_1_child_1_;	// Mesh of the ellipse
						ell_child.Parent_E_name			= 	cEll_1.Parent_E_name;			// Name of the parent ellipse

						if (VERBOSE == 1){
							std::cout << "child_polys.poly_1_child_1_.num_vertices() = "
									<< child_polys.poly_1_child_1_.num_vertices() << std::endl;
						}

						ell_child_list_new.push_back(ell_child);

						ell_child.E_mesh			=	child_polys.poly_1_child_2_;			// Mesh of the ellipse
						ell_child_list_new.push_back(ell_child);

						CGAL::update_polygon_numbering <Kernel, Ell_>(ell_child_list_new);

						if (VERBOSE == 1){
							std::cout << "ell_child_list_new.size() = " << ell_child_list_new.size() << std::endl;
							//ell_child_list.clear();

							//ell_child_list = ell_child_list_new;

							//std::cout << "ell_child_list.size() = " << ell_child_list.size() << std::endl;
							//ell_child_list_new.clear();
							std::cout << "child_polys.intersectionline_.size() = " << child_polys.intersectionline_.size() << std::endl;
						}
						for (Pts_Iterator it = child_polys.intersectionline_.begin(); it !=child_polys.intersectionline_.end(); it++){
							// 23/09/2016 fname_il_ofs << (*it) << " " << i << " " << j <<"\n";
							fname_il_ofs  << std::setprecision(30)  << (*it) << " " << i << " " << j <<"\n";
						}
					}

					else{
						// Add the polygon to the new list
						ell_child_list_new.push_back(cEll_child);
					}

					for (int l = k+1; l < ell_child_list.size(); l++ ){
						std::cout << "Push back the remaining ellipses : k = " << k << ", l = " << l << ", ell_child_list.size() = " <<ell_child_list.size() << std::endl;
						Ell_Iterator cEll_child_iter_re = ell_child_list.begin();
						std::advance(cEll_child_iter_re, l);
						cEll_child_re = *cEll_child_iter_re;
						ell_child_list_new.push_back(cEll_child_re);
					}

					//CGAL::export_vtu_files_from_ell_list<Ell_, Surface_mesh>(ell_child_list_new, std::string("output/ellipse_test_") );

					ell_child_list.clear();
					ell_child_list = ell_child_list_new;

					std::cout << "ell_child_list.size() = " << ell_child_list.size() << std::endl;

					ell_child_list_new.clear();


					if (child_polys.no_polygons_interst_ == false){
						k=k+2;
					}
					else{
						k++;
					}
				}
				CGAL::silent_line('-',60);

			}
		}
		ell_child_list_list.push_back(ell_child_list);
		//ell_child_list_new.clear();
		ell_child_list.clear();
		CGAL::silent_line('=',80);
	}

	fname_il_ofs.close();

	//std::cout << "ell_child_list.size() = " << ell_child_list.size() << std::endl;
	//std::cout << "ell_child_list_list.size() = " << ell_child_list_list.size() << std::endl;

	if (ell_list.size() == 1){
		full_ell.push_back(*ell_list.begin());
	}
	else
	{
		// Add fractures to the final mesh and remove disconnected fractures
		for (Ell_list_Iterator it1 = ell_child_list_list.begin(); it1 != ell_child_list_list.end(); it1++){
			std::cout << "*it1.size() = " << (*it1).size() << std::endl;
			if ((*it1).size() > 1){
				for (Ell_Iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++){
					full_ell.push_back(*it2);
				}
			}
		}

		std::cout << "full_ell.size() = " << full_ell.size() << std::endl;

		if (full_ell.size() == 0){
			std::cerr << "No connected fractures in the domain. Change DFN domain size in the parameter file \"" << inputFile <<
					"\" or input data of fractures \"" << inputData << "\"."<< std::endl;
			exit(-1);
		}

	}

	// Update polygon numbering of full_ell
	full_ell = CGAL::update_polygon_numbering <Kernel, Ell_>(full_ell);

	// Update target edge length of full_ell

	if (!RM_CLOSE_PTS){
		full_ell = CGAL::update_target_edge_length <Kernel, Ell_>(full_ell);
	}
	/*
	else{
		if (iparams.hasminLengthRatio){
			full_ell = CGAL::update_target_edge_length_rm_pts <Kernel, Ell_>(full_ell, 1 + iparams.minLengthRatio);
		}
		else{
			full_ell = CGAL::update_target_edge_length_rm_pts <Kernel, Ell_>(full_ell, 1.1);
		}
	}
	*/

	time(&stop_cutting_off);
	clock_t end_cutting_off = clock();
	double CPU_cutting_off_secs = double(end_cutting_off - end_check_intersection) / CLOCKS_PER_SEC;

	// =====================================================================================================
	//                                 Prepare input files for lagrit run
	// =====================================================================================================

	// Split the intersection file
	fname = std::string("intersections/intersections.inp");
	int numParentPoly = ell_list.size();
	CGAL::split_intersection_file <NT_MP> (fname, numParentPoly);

	double min_target_edge_length = (*full_ell.begin()).E_target_edge_length;

	for (int i = 0; i !=full_ell.size(); i++){
		Ell_			cEll_1;			// Current ellipse
		Ell_Iterator	cEll_1_iter = full_ell.begin();
		std::advance(cEll_1_iter, i);
		cEll_1 = *cEll_1_iter;

		if (VERBOSE == 1){
			// Print out the name of polygon and it parent
			CGAL::silent_line('=',80);
			std::cout << "Polygon name " << cEll_1.E_name << std::endl;
			std::cout << "Polygon parent name " << cEll_1.Parent_E_name << std::endl;
			std::cout << "Polygon normal vector " << cEll_1.E_normalvect << std::endl;
		}
		std::string	parentfilename = "intersections/intersections_parent_" + CGAL::int2str_setw(cEll_1.Parent_E_name, CGAL::count_digits(ell_list.size())) + ".inp";
		std::string	childfilename =  "intersections/intersections_" + CGAL::int2str_setw(i + 1, CGAL::count_digits(full_ell.size())) + ".inp";

		cmd = "cp " + parentfilename + " " + childfilename;
		std::system (cmd.c_str());

		// The case where there is only a fracture in the DFN

		if (full_ell.size() == 1){
			std::ofstream of_stream_file (childfilename.c_str());
			of_stream_file << "0 0 0 0 0\n";
		}
		if (min_target_edge_length > cEll_1.E_target_edge_length){
			min_target_edge_length = cEll_1.E_target_edge_length;
		}

		std::ofstream    parentName_file;
		parentName_file.open("output/poly_reorderedparentName.dat", std::ios_base::app);
		parentName_file << cEll_1.Parent_E_name << std::endl;
		parentName_file.close();
	}

	std::cout << "min_target_edge_length = " << min_target_edge_length << std::endl;

	//std::system ("rm -rf intersections/intersections_parent_*");

	std::string poly_file(std::string("output/lagrit/ellipse_cutoff_"));

	// Export vtu files from the ellipse list
	CGAL::export_vtu_files_from_ell_list<Ell_, Surface_mesh>(full_ell, std::string("output/ellipse_cutoff/ellipse_cutoff_") );
	//CGAL::export_inp_files_from_ell_list<Ell_, Surface_mesh>(full_ell, std::string("output/lagrit/ellipse_cutoff_orig_"), CGAL::count_digits(full_ell.size()));

	CGAL::write_params_file_for_lagrit<Ell_, Kernel, NT_MP>(full_ell, p_min, p_max, std::string("output/params.txt"));
	CGAL::write_normal_vertors<Ell_, Kernel, NT_MP>(full_ell, std::string("output/lagrit/normal_vectors.txt"));
	CGAL::export_inp_files_from_ell_list_xy_plane<Ell_, Kernel, NT_MP>(full_ell, poly_file, p_min, p_max, CGAL::count_digits(full_ell.size()));

	std::cout <<"iparams.hasinverseSignTheta = " << iparams.hasinverseSignTheta << std::endl;
	std::cout <<"iparams.inverseSignTheta = " << iparams.inverseSignTheta << std::endl;
	std::cout <<"iparams.hasparamFilePrecision = " << iparams.hasparamFilePrecision << std::endl;
	std::cout <<"iparams.paramFilePrecision = " << iparams.paramFilePrecision << std::endl;

	//exit(-1);

	if (iparams.hasparamFilePrecision && iparams.paramFilePrecision){
		if (iparams.hasinverseSignTheta && !iparams.inverseSignTheta){
			//CGAL::create_parameter_mlgi_files<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
			//std::cout << "==> Running dfn generator for vercors reservoir." << std::endl;
			CGAL::create_parameter_mlgi_files_vercors_setprecision<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
		}
		else{
			//cout << "vercors: 0" << endl;
			CGAL::create_parameter_mlgi_files_setprecision<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
		}
	}
	else{

		if (iparams.hasinverseSignTheta && !iparams.inverseSignTheta){
			//CGAL::create_parameter_mlgi_files<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
			//std::cout << "vercors: 1" << std::endl;
			CGAL::create_parameter_mlgi_files_vercors<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
		}
		else{
			//cout << "vercors: 0" << endl;
			CGAL::create_parameter_mlgi_files<Ell_, NT_MP>(full_ell, std::string("output/params.txt"));
		}
	}

	std::string checkfolder("output/parameters");
	if ( stat( checkfolder.c_str(), &info) != 0){
		std::cout << "stat( checkfolder.c_str(), &info) != 0" << std::endl;

		if (iparams.hasparamFilePrecision && iparams.paramFilePrecision){
			if (iparams.hasinverseSignTheta && !iparams.inverseSignTheta){
				cout << "==> Running correcting mesh part of dfn generator for vercors reservoir." << endl;
				//std::system("pwd");
				//std::system("ln -s Executables/exe_create_mlgi_vercors_setpr exe_create_mlgi_vercors_setpr");
				/*
				std::system("cp Executables/exe_create_mlgi_vercors_setpr .");
				std::cout << "OK 1" << std::endl;
				std::system("./exe_create_mlgi_vercors_setpr");
				std::cout << "OK 2" << std::endl;
				*/

				std::system("ln -s Executables/exe_create_mlgi_vercors_setpr exe_create_mlgi_vercors_setpr");
				std::system("./exe_create_mlgi_vercors_setpr");

				//std::system("rm -rf exe_create_parameter_mlgi_files_vercors");
			}
			else{
				std::system("ln -s Executables/exe_create_mlgi_setpr exe_create_mlgi_setpr");
				std::system("./exe_create_mlgi_setpr");
				//std::system("rm -rf create_parameter_mlgi_files");
			}
		}
		else{
			if (iparams.hasinverseSignTheta && !iparams.inverseSignTheta){
				//cout << "==> Running dfn generator for vercors reservoir." << endl;
				//std::system("pwd");
				std::system("ln -s Executables/exe_create_mlgi_vercors exe_create_mlgi_vercors");
				std::system("./exe_create_mlgi_vercors");
				//std::system("rm -rf exe_create_parameter_mlgi_files_vercors");
			}
			else{
				std::system("ln -s Executables/exe_create_mlgi exe_create_mlgi");
				std::system("./exe_create_mlgi");
				//std::system("rm -rf create_parameter_mlgi_files");
			}
		}
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Update "ell.Parent_E_name" for LaGriT DFNWorks - PTDFNTrans run
	std::list<int> parentName_In_list, parentName_Out_list;
	for (int i = 0; i !=full_ell.size(); i++){
		Ell_Iterator	cEll_1_iter = full_ell.begin();
		std::advance(cEll_1_iter, i);
		parentName_In_list.push_back((*cEll_1_iter).Parent_E_name);
	}

	parentName_In_list.sort();
	parentName_In_list.unique();

	std::cout << "ell_list.size() = " <<  ell_list.size() << std::endl;

	for (int i = 0; i !=parentName_In_list.size(); i++){
		Int_Iterator parentName_iter = parentName_In_list.begin();
		std::advance(parentName_iter, i);
		//std::cout << "parentName_In_list[" << i << "] = " << *parentName_iter << std::endl;
	}
	// Return disconnected fractures
	for (int j=1; j < *parentName_In_list.begin(); j++){
		parentName_Out_list.push_back(j);
	}

	for (int i = 0; i !=parentName_In_list.size()-1; i++){
		Int_Iterator parentName_curr = parentName_In_list.begin();
		std::advance(parentName_curr, i);
		Int_Iterator parentName_next = parentName_In_list.begin();
		std::advance(parentName_next, i+1);

		for (int j=1; j < (*parentName_next - *parentName_curr); j++){
			parentName_Out_list.push_back(*parentName_curr + j);
		}
	}

	for (int j=*parentName_In_list.end(); j < ell_list.size(); j++){
		parentName_Out_list.push_back(j+1);
	}

	std::ofstream    parentName_Out_list_file("parentName_Out_list.dat");

	for (int i = 0; i !=parentName_Out_list.size(); i++){
		Int_Iterator parentName_iter = parentName_Out_list.begin();
		std::advance(parentName_iter, i);
		//std::cout << "parentName_Out_list[" << i << "] = " << *parentName_iter << std::endl;
		parentName_Out_list_file << *parentName_iter << std::endl;

	}
	parentName_Out_list_file.close();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// 									    [IMPORTANT]
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// Create lagrit control files     <--------------------- [IMP]
	CGAL::create_lagrit_scripts_list(full_ell);

	// Path to lagrit executable
	std::string lagrit_path( std::string("/work/irlin104_1/ngotr/Outils/lagrit/lagrit_ulin3.2") );


	// Using multiprocessing
	int npoly(full_ell.size());
	int part_size(npoly);
	if (npoly%ncpu == 0){
		part_size = npoly/ncpu;
	}
	else{
		part_size = round(npoly/ncpu) + 1;
	}
	std::list<int> endis;
	endis.push_back(0);

	for (int i = 0; i!= ncpu-1; i++){
		endis.push_back(part_size * (i+1));
	}
	endis.push_back(npoly);

	if (VERBOSE == 1){
		std::cout << "ncpu = " << ncpu << std::endl;
		std::cout << "npoly = " << npoly << std::endl;
		std::cout << "part_size = " << part_size << std::endl;

		for (Int_Iterator it= endis.begin(); it != endis.end(); it++){
			std::cout << "endis = " << *it << std::endl;
		}
	}

	// Run lagrit script to triangulate polygons
	bool dataForDFNWorks(false), correctMesh(true);
	if (iparams.hasBlockParametersLaGriT){

		if (iparams.hascorrectMesh){
			correctMesh = iparams.correctMesh;
		}
		else{
			std::cout << "** Warning: BlockParametersLaGriT is set but the program failed to parse correctMesh value. "
					  << "** Default \"correctMesh = true\" will be set." << std::endl;
		}

		if (iparams.hasdataForDFNWorks){
			dataForDFNWorks = iparams.dataForDFNWorks;
		}
		else{
			std::cout << "** Warning: BlockParametersLaGriT is set but the program failed to parse dataForDFNWorks value.\n"
					  << "** Default \"dataForDFNWorks = false\" will be set." << std::endl;
		}
	}

	//std::cout << "** hasdataForDFNWorks = "<< iparams.hasdataForDFNWorks << std::endl;
	//std::cout << "** dataForDFNWorks = "<< dataForDFNWorks << std::endl;
	//std::cout << "** iparams.dataForDFNWorks = "<< iparams.dataForDFNWorks << std::endl;

	CGAL::triangulate_polygons(lagrit_path, full_ell.size(), correctMesh);

	// Create merge_poly file
	CGAL::create_merge_poly_files_list_multi(full_ell, endis);

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

	// Create merge_rmpts
	CGAL::create_merge_rmpts_file_multi(EPS_INT, EPS_FILTER, endis);

	// Run merge_rmpts twice
	CGAL::run_merge_mesh_scripts_multi(lagrit_path,endis);

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
			CGAL::create_merge_rmpts_file_multi(EPS_INT, EPS_FILTER, endis);
			// Run merge_rmpts.
			CGAL::run_merge_mesh_scripts_multi(lagrit_path, endis);
			check_quality_voronoi_mesh_new = CGAL::check_quality_voronoi_mesh();
		}
	}

	CGAL::redefine_zones(lagrit_path,EPS_FILTER);

	time(&stop_lagrit_run);
	clock_t end_lagrit_run = clock();
	double CPU_lagrit_run_secs = double(end_lagrit_run - end_cutting_off) / CLOCKS_PER_SEC;

	// =====================================================================================================
	//                                      Write input file for DFNWorks run
	// =====================================================================================================

	int dimScaling = 1, intRound = 12;
	if (iparams.hasdimScaling){
		dimScaling = iparams.dimScaling;
	}

	if (dataForDFNWorks){
		std::system("rm -rf part_*.lg");
		/*
		CGAL::removeNegativeVoronoiCells_multi(full_ell);
		cmd = lagrit_path + std::string(" < merge_poly_new.lgi > lagrit_logs/log_merge_poly_new");
		std::system(cmd.c_str());

		cmd = lagrit_path + std::string(" < merge_rmpts.lgi > lagrit_logs/log_merge_all_new");
		std::system(cmd.c_str());
		 */

		CGAL::removeNegativeVoronoiCells(full_ell);
		for (int j = 1; j!=endis.size(); j++){
			cmd = lagrit_path + std::string(" < merge_poly_part_" + CGAL::int2str(j) + "_new.lgi > lagrit_logs/log_merge_poly_part_" + CGAL::int2str(j) +"_new");
			std::system(cmd.c_str());
		}

		cmd = lagrit_path + std::string(" < merge_rmpts_multi.lgi > lagrit_logs/log_merge_all_new");
		std::system(cmd.c_str());

	}

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
	std::system("mkdir output/meshes");
	std::system("mkdir output/polymeshes");
	std::system("mv mesh*.inp ./output/polymeshes");
	std::system("mv mesh*.gmv ./output/meshes");
	std::system("mv logx3dgen ./output");
	std::system("mv outx3dgen ./output");

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

	// Write aperture and permeability parameters to files
	if (1){
		std::ofstream    outfile_aper(apertureFile.c_str());
		for (int i  = 0; i != ell_list.size(); i++){
			outfile_aper << "1E-03\n";
		}
		outfile_aper.close();

		std::ofstream    outfile_perm(permFile.c_str());
		for (int i  = 0; i != ell_list.size(); i++){
			outfile_perm << "1E-11 1E-11 1E-11\n";
		}
		outfile_perm.close();

		std::string apertureFile_DFNWorks(std::string("aperture_DFNWorks.dat"));
		std::string permFile_DFNWorks(std::string("permeability_DFNWorks.dat"));

		std::ofstream    outfile_aper_DFNWorks(apertureFile_DFNWorks.c_str());
		outfile_aper_DFNWorks << "aperture\n";
		for (int i  = 0; i != ell_list.size(); i++){
			outfile_aper_DFNWorks << "-" << i+7 << " 0 0  1E-03\n";
		}
		outfile_aper_DFNWorks.close();

		std::ofstream    outfile_perm_DFNWorks(permFile_DFNWorks.c_str());
		outfile_perm_DFNWorks << "permeability\n";
		for (int i  = 0; i != ell_list.size(); i++){
			outfile_perm_DFNWorks << "-" << i+7 << " 0 0  1E-11 1E-11 1E-11\n";
		}
		outfile_perm_DFNWorks.close();
	}

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
			std::system ("rm -rf intersections_*.inp *.stor");
			//std::system("rm -rf ./output/meshes/ -rf ./output/polymeshes/");
			std::system ("rm -rf ./output/meshes/");
			std::system ("rm -rf full_mesh_orig.vtk full_mesh.gmvF");
		}
		else if (iparams.cleanupDirectory == 2){
			std::system ("rm -rf mesh_*.gmv *.lgi");
			std::system ("rm -rf ellipse_cutoff_*.inp *.mlgi intersections_*.inp");
			std::system ("rm -rf intersections_*.inp *.stor *.zone part_*.lg");
			std::system ("rm -rf ./output/meshes/ -rf ./output/polymeshes/");
			std::system ("rm -rf full_mesh_orig.vtk full_mesh.gmvF");

		}
	}

	/*
	std::cout << " " << std::endl;
	cout << "total_intersection_number = " << total_intersection_number << endl;
	*/

	// =====================================================================================================
	//                                      Export CPU times
	// =====================================================================================================
    time(&stop);
	clock_t end = clock();
	double CPU_post_treatment_secs = double(end - end_lagrit_run) / CLOCKS_PER_SEC;
	double CPU_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "\n";
	std::cout << "CPU_secs = " << CPU_secs << std::endl;
	/*
	std::cout << "CPU_input_secs = " << CPU_input_secs << std::endl;
	std::cout << "CPU_out_of_dfnbb_secs = " << CPU_out_of_dfnbb_secs << std::endl;
	std::cout << "CPU_intersection_secs = " << CPU_intersection_secs << std::endl;
	std::cout << "CPU_cutting_off_secs = " << CPU_cutting_off_secs << std::endl;
	std::cout << "CPU_lagrit_run_secs = " << CPU_lagrit_run_secs << std::endl;
	std::cout << "CPU_post_treatment_secs = " << CPU_post_treatment_secs << std::endl;
	*/

	std::ofstream    fname_CPU_time("output/CPU_time.txt");
	fname_CPU_time << "# CPU_input_secs" << std::endl;
	fname_CPU_time << "# CPU_out_of_dfnbb_secs" << std::endl;
	fname_CPU_time << "# CPU_intersection_secs" << std::endl;
	fname_CPU_time << "# CPU_cutting_off_secs" << std::endl;
	fname_CPU_time << "# CPU_lagrit_run_secs" << std::endl;
	fname_CPU_time << "# CPU_post_treatment_secs" << std::endl;
	fname_CPU_time << "# CPU_secs" << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_input_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_out_of_dfnbb_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_intersection_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_cutting_off_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_lagrit_run_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_post_treatment_secs << std::endl;
	fname_CPU_time << std::setprecision(12) << std::scientific << CPU_secs << std::endl;
	fname_CPU_time.close();

	std::ofstream    fname_CPU_time_percent("output/CPU_time_percent.txt");
	fname_CPU_time_percent << "# CPU_input_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_out_of_dfnbb_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_intersection_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_cutting_off_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_lagrit_run_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_post_treatment_secs" << std::endl;
	fname_CPU_time_percent << "# CPU_secs" << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_input_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_out_of_dfnbb_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_intersection_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_cutting_off_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_lagrit_run_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_post_treatment_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent << std::setprecision(6) << CPU_secs/CPU_secs*100 << std::endl;
	fname_CPU_time_percent.close();

    std::ofstream    fname_elapsed_time("output/elapsed_time.txt");
    fname_elapsed_time << "# elapsed_input_secs" << std::endl;
    fname_elapsed_time << "# elapsed_out_of_dfnbb_secs" << std::endl;
    fname_elapsed_time << "# elapsed_intersection_secs" << std::endl;
    fname_elapsed_time << "# elapsed_cutting_off_secs" << std::endl;
    fname_elapsed_time << "# elapsed_lagrit_run_secs" << std::endl;
    fname_elapsed_time << "# elapsed_post_treatment_secs" << std::endl;
    fname_elapsed_time << "# elapsed_secs" << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop_input, start) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop_check_out_of_dfnbb, stop_input) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop_check_intersection, stop_check_out_of_dfnbb) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop_cutting_off, stop_check_intersection) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop_lagrit_run,stop_cutting_off) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop, stop_lagrit_run) << std::endl;
    fname_elapsed_time << std::setprecision(12) << std::scientific << difftime(stop, start) << std::endl;
    fname_elapsed_time.close();

    std::ofstream    fname_elapsed_time_percent("output/elapsed_time_percent.txt");
    fname_elapsed_time_percent << "# elapsed_input_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_out_of_dfnbb_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_intersection_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_cutting_off_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_lagrit_run_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_post_treatment_secs" << std::endl;
    fname_elapsed_time_percent << "# elapsed_secs" << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop_input, start)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop_check_out_of_dfnbb, stop_input)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop_check_intersection, stop_check_out_of_dfnbb)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop_cutting_off, stop_check_intersection)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop_lagrit_run,stop_cutting_off)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop, stop_lagrit_run)/difftime(stop, start)*100 << std::endl;
    fname_elapsed_time_percent << std::setprecision(6) << difftime(stop, start)/difftime(stop, start)*100 << std::endl;

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
