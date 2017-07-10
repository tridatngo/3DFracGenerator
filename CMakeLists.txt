# This is the CMake script for compiling a CGAL application.

project( Create_DFN_example )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  #create_single_source_cgal_program( "mesh_a_sphere_by_surface_mesher.cpp" )
  #create_single_source_cgal_program( "create_an_ellipse.cpp" )
  #create_single_source_cgal_program( "create_an_ellipse_implicit.cpp" )
  #create_single_source_cgal_program( "example_polygon_2D.cpp" )
  #create_single_source_cgal_program( "example_polygon_2D3D.cpp" )
  #create_single_source_cgal_program( "sm_iterators.cpp" )
  #create_single_source_cgal_program( "triangulate_faces_example.cpp" )
  #create_single_source_cgal_program( "mesh_a_3d_gray_image.cpp" )
  #create_single_source_cgal_program( "mesh_an_implicit_function.cpp" )
  #create_single_source_cgal_program( "refine_fair_example.cpp")
  #create_single_source_cgal_program( "isotropic_remeshing_example.cpp" )
  #create_single_source_cgal_program( "Surface_Mesh_triangulate_example.cpp" )
  #create_single_source_cgal_program( "Intersection_SM_triangulate_example.cpp" )
  #create_single_source_cgal_program( "AABB_Intersection_SM.cpp" )
  create_single_source_cgal_program( "exeCreateDFN.cpp" )
  create_single_source_cgal_program( "CorrectDFN.cpp" )

  if (EIGEN3_FOUND)
     include( ${EIGEN3_USE_FILE} )
     create_single_source_cgal_program( "test_intersection.cpp" )
  else()
     message(STATUS "NOTICE: These examples require the Eigen library, version 3.2 or later and will not be compiled.")
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

