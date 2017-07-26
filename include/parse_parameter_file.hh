/*
 * parse_parameter_file.hh
 *
 *  Created on: Sep 19, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_PARSE_PARAMETER_FILE_HH_
#define INCLUDE_PARSE_PARAMETER_FILE_HH_

#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<cstring>
#include <iomanip> // setprecision
#include "defaultusagemessage.hh"

namespace CGAL{
template <typename InputParams>
InputParams parse_parameter_file(const char* inputFile, char ** argv){
	InputParams iparams;
	bool    hasBlockDensity (false), hasBlockSpatialParams(false),
			hasBlockAddingPointsMethod(false), hasBlockDiscreMethod(false),
			hasBlockRemoveClosePoints;
	bool 	hasDFNBB_LowerLeftX(false),
			hasDFNBB_LowerLeftY(false),
			hasDFNBB_LowerLeftZ(false),
			hasDFNBB_UpperRightX(false),
			hasDFNBB_UpperRightY(false),
			hasDFNBB_UpperRightZ(false);
	bool 	hasBlockCleanUpDirectory(false), hascleanupDirectory(false),
			hasBlockDataFile(false), hasdataFile(false), hasdesactFracsFile(false),
			hasBlockParametersLaGriT(false), hasEPS_INT(false), hasEPS_FILTER(false), hasEPS_BOUNDARY(false),
			hasminLengthRatio(false);

	bool	hasn4(false),
			hasaddingPointsMethod(false),
			hasdiscreMethod(false),
			hasdimScaling(false),
			hasintRound(false),
			hasremoveClosePoints(false),
			hasmergeCloseIntersectionPoints(false),
			hasmergeClosePointsRelCritLength(false),
			hascorrectMesh(false),
			hasdataForDFNWorks(false),
			hasinverseSignTheta(false),
			hasparamFilePrecision(false);

	int n4, addingPointsMethod, discreMethod, dimScaling, intRound, cleanupDirectory,
		int_rm_close_pts, int_merge_close_intersection_pts, int_correctMesh,
		int_dataForDFNWorks, int_inverseSignTheta, int_paramFilePrecision;
	bool rm_close_pts, merge_close_intersection_pts(false), correctMesh, dataForDFNWorks, inverseSignTheta, paramFilePrecision;
	std::string dataFile, desactFracsFile("defaultFile.dat");
	double 	DFNBB_LowerLeftX,
	DFNBB_LowerLeftY,
	DFNBB_LowerLeftZ,
	DFNBB_UpperRightX,
	DFNBB_UpperRightY,
	DFNBB_UpperRightZ;

	double EPS_INT, EPS_FILTER, EPS_BOUNDARY, minLengthRatio, mergeClosePointsRelCritLength;

	std::ifstream    infile(inputFile);
	if (! infile.is_open()) {
		std::cerr << "Failed to open the input file." << std::endl;
		std::cout << CGAL::defaultUsageMessage(std::string(argv[0])) << "\n";
		exit(-1);
	}

	std::string line, linenospace;
	std::istringstream ist(line);

	//infile.open(inputFile);
	while (!infile.eof())
	{
		std::getline(infile, line);
		linenospace = line;
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
		linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());
		//std::cout << " Line : " << line << std::endl;
		//std::cout << " Line no space : " << linenospace << std::endl;

		if (linenospace.substr (0,1) != "#" ){
			// Search block
			if (linenospace.find("[Density]") != std::string::npos){
				hasBlockDensity = true;
			}

			if (linenospace.find("[SpatialParams]") != std::string::npos){
				hasBlockSpatialParams = true;
			}

			if (linenospace.find("[AddingPointsMethod]") != std::string::npos){
				hasBlockAddingPointsMethod = true;
			}

			if (linenospace.find("[DiscreMethod]") != std::string::npos){
				hasBlockDiscreMethod = true;
			}

			if (linenospace.find("[CleanUpDirectory]") != std::string::npos){
				hasBlockCleanUpDirectory = true;
			}

			if (linenospace.find("[DataFile]") != std::string::npos){
				hasBlockDataFile = true;
			}

			if (linenospace.find("[ParametersLaGriT]") != std::string::npos){
				hasBlockParametersLaGriT = true;
			}

			if (linenospace.find("[RemoveClosePoints]") != std::string::npos){
				hasBlockRemoveClosePoints = true;
			}

			// Search parameter
			if (linenospace.find("n4=") != std::string::npos){
				hasn4 = true;
			}

			if (linenospace.find("DFNBB_LowerLeftX=") != std::string::npos){
				hasDFNBB_LowerLeftX = true;
			}


			if (linenospace.find("DFNBB_LowerLeftY=") != std::string::npos){
				hasDFNBB_LowerLeftY = true;
			}

			if (linenospace.find("DFNBB_LowerLeftZ=") != std::string::npos){
				hasDFNBB_LowerLeftZ = true;
			}

			if (linenospace.find("DFNBB_UpperRightX=") != std::string::npos){
				hasDFNBB_UpperRightX = true;
			}

			if (linenospace.find("DFNBB_UpperRightY=") != std::string::npos){
				hasDFNBB_UpperRightY = true;
			}


			if (linenospace.find("DFNBB_UpperRightZ=") != std::string::npos){
				hasDFNBB_UpperRightZ = true;
			}


			if (linenospace.find("addingPointsMethod") != std::string::npos){
				hasaddingPointsMethod = true;
			}

			if (linenospace.find("discreMethod") != std::string::npos){
				hasdiscreMethod = true;
			}

			if (linenospace.find("dimScaling") != std::string::npos){
				hasdimScaling = true;
			}

			if (linenospace.find("intRound") != std::string::npos){
				hasintRound = true;
			}

			if (linenospace.find("cleanupDirectory") != std::string::npos){
				hascleanupDirectory = true;
			}

			if (linenospace.find("dataFile") != std::string::npos){
				hasdataFile = true;
			}

			if (linenospace.find("desactFracsFile") != std::string::npos){
				hasdesactFracsFile = true;
			}

			if (linenospace.find("EPS_INT") != std::string::npos){
				hasEPS_INT = true;
			}

			if (linenospace.find("EPS_FILTER") != std::string::npos){
				hasEPS_FILTER = true;
			}

			if (linenospace.find("EPS_BOUNDARY") != std::string::npos){
				hasEPS_BOUNDARY = true;
			}

			if (linenospace.find("removeClosePoints") != std::string::npos){
				hasremoveClosePoints = true;
			}

			if (linenospace.find("mergeCloseIntersectionPoints") != std::string::npos){
				hasmergeCloseIntersectionPoints = true;
			}

			if (linenospace.find("mergeClosePointsRelCritLength") != std::string::npos){
				hasmergeClosePointsRelCritLength = true;
			}

			if (linenospace.find("minLengthRatio") != std::string::npos){
				hasminLengthRatio = true;
			}

			if (linenospace.find("correctMesh") != std::string::npos){
				hascorrectMesh = true;
			}

			if (linenospace.find("dataForDFNWorks") != std::string::npos){
				hasdataForDFNWorks = true;
			}

			if (linenospace.find("inverseSignTheta") != std::string::npos){
				hasinverseSignTheta = true;
			}

			if (linenospace.find("paramFilePrecision") != std::string::npos){
				hasparamFilePrecision = true;
			}
		}
	}

	infile.close();

	// Parse n4
	if (hasn4){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			//std::cout << " linenospace.substr (0,1) : " << linenospace.substr (0,1) << std::endl;
			//std::cout << " Line : " << line << std::endl;
			//std::cout << " Line no space : " << linenospace << std::endl;

			if (linenospace.substr (0,1) != "#" && linenospace.find("n4=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}
				/*
				std::cout << " first : " << first << std::endl;
				std::cout << " last : " << last << std::endl;
				 */

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> n4;

				break;
			}
		}
		infile.close();
	}

	// Parse spatial parameters
	if (hasDFNBB_LowerLeftX){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_LowerLeftX=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_LowerLeftX;

				break;
			}
		}
		infile.close();
	}

	if (hasDFNBB_LowerLeftY){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_LowerLeftY=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_LowerLeftY;

				break;
			}
		}
		infile.close();
	}

	if (hasDFNBB_LowerLeftZ){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_LowerLeftZ=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_LowerLeftZ;

				break;
			}
		}
		infile.close();
	}

	if (hasDFNBB_UpperRightX){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_UpperRightX=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_UpperRightX;

				break;
			}
		}
		infile.close();
	}

	if (hasDFNBB_UpperRightY){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_UpperRightY=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_UpperRightY;

				break;
			}
		}
		infile.close();
	}

	if (hasDFNBB_UpperRightZ){
		infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			linenospace = line;
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
			linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

			if (linenospace.substr (0,1) != "#" && linenospace.find("DFNBB_UpperRightZ=") != std::string::npos){
				unsigned first = line.find("=");
				unsigned last = line.size();
				if (linenospace.find("#") != std::string::npos){
					last = line.find("#")-1;
				}

				std::string strNew = line.substr (first+1,last-first);
				std::istringstream sstream(strNew);
				sstream >> DFNBB_UpperRightZ;

				break;
			}
		}
		infile.close();
	}

	// Parse addingPointsMethod
	if (hasBlockAddingPointsMethod){
		if(hasaddingPointsMethod){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("addingPointsMethod=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> addingPointsMethod;

					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cerr << "Error : AddingPointsMethod block is set but the program failed to parse addingPointsMethod value. "<<std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
	}

	// Parse discreMethod
	if (hasBlockDiscreMethod){
		if(hasdiscreMethod){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("discreMethod=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> discreMethod;

					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cerr << "Error : DiscreMethod block is set but the program failed to parse discreMethod value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
	}

	// Parse cleanupDirectory
	if (hasBlockCleanUpDirectory){
		if(hascleanupDirectory){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("cleanupDirectory=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> cleanupDirectory;
					/*
						if (intcleanupDirectory == 0){
							cleanupDirectory = false;
						}
						else if (intcleanupDirectory != 0){
							cleanupDirectory = true;
						}
					 */
					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cerr << "Error : CleanUpDirectory block is set but the program failed to parse cleanupDirectory value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
	}

	// Parse rm_close_pts
	if (hasBlockRemoveClosePoints){
		if(hasremoveClosePoints){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("removeClosePoints=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_rm_close_pts;

					if (int_rm_close_pts == 0){
						rm_close_pts = false;
					}
					else if (int_rm_close_pts != 0){
						rm_close_pts = true;
					}

					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cerr << "Error : RemoveClosePoints block is set but the program failed to parse removeClosePoints value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}

		if(hasminLengthRatio){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("minLengthRatio=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> minLengthRatio;
					break;
				}
			}
			infile.close();
		}

		if(hasmergeCloseIntersectionPoints){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("mergeCloseIntersectionPoints=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_merge_close_intersection_pts;

					if (int_merge_close_intersection_pts == 0){
						merge_close_intersection_pts = false;
					}
					else if (int_merge_close_intersection_pts != 0){
						merge_close_intersection_pts = true;
					}

					break;
				}
			}
			infile.close();
		}

		if(hasmergeClosePointsRelCritLength){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("mergeClosePointsRelCritLength=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> mergeClosePointsRelCritLength;
					break;
				}
			}
			infile.close();
		}
	}

	// Parse dataFile
	if (hasBlockDataFile){
		if(hasdataFile){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("dataFile=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					strNew.erase(remove(strNew.begin(), strNew.end(), '\t'), strNew.end());
					strNew.erase(remove(strNew.begin(), strNew.end(), ' '), strNew.end());
					//dataFile = strNew.c_str();
					dataFile = strNew;
					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cerr << "Error : DataFile block is set but the program failed to parse dataFile value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}


		if(hasdesactFracsFile){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("desactFracsFile=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					strNew.erase(remove(strNew.begin(), strNew.end(), '\t'), strNew.end());
					strNew.erase(remove(strNew.begin(), strNew.end(), ' '), strNew.end());
					desactFracsFile = strNew;
					break;
				}
			}
			infile.close();
		}
	}

	// Parse LaGriT parameters
	if (hasBlockParametersLaGriT){
		if(hasEPS_FILTER){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("EPS_FILTER=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> EPS_FILTER;
					//std::cout << "EPS_FILTER = " << std::setprecision(3) <<std::scientific << EPS_FILTER << std::endl;
					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cout << "Warning : ParametersLaGriT block is set but the program failed to parse EPS_FILTER value. " << std::endl;
			std::cout << "The default value will be taken or please correct the input file (" << inputFile << ")." << std::endl;
			//exit(-1);
		}


		if(hasEPS_BOUNDARY){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("EPS_BOUNDARY=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> EPS_BOUNDARY;
					break;
				}
			}
			infile.close();
		}

		if(hasEPS_INT){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("EPS_INT=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> EPS_INT;
					//std::cout << "EPS_FILTER = " << std::setprecision(3) <<std::scientific << EPS_FILTER << std::endl;
					break;
				}
			}
			infile.close();
		}
		else
		{
			std::cout << "Warning : ParametersLaGriT block is set but the program failed to parse EPS_INT value. " << std::endl;
			std::cout << "The default value will be taken or please correct the input file (" << inputFile << ")." << std::endl;
			//exit(-1);
		}

		if(hascorrectMesh){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				//std::cout << "linenospace: " << linenospace <<std::endl;

				if (linenospace.substr (0,1) != "#" && linenospace.find("correctMesh=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_correctMesh;

					//std::cout << "int_correctMesh = " << int_correctMesh <<std::endl;

					if (int_correctMesh == 0){
						correctMesh = false;
					}
					else if (int_correctMesh != 0){
						correctMesh = true;
					}
					break;
				}
			}

			infile.close();
		}
		else
		{
			std::cout << "Warning : ParametersLaGriT block is set but the program failed to parse correctMesh value. " << std::endl;
			std::cout << "The default value will be taken or please correct the input file (" << inputFile << ")." << std::endl;
			//exit(-1);
		}

		if(hasdataForDFNWorks){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				//std::cout << "linenospace: " << linenospace <<std::endl;

				if (linenospace.substr (0,1) != "#" && linenospace.find("dataForDFNWorks=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_dataForDFNWorks;

					//std::cout << "int_dataForDFNWorks = " << int_correctMesh <<std::endl;

					if (int_dataForDFNWorks == 0){
						dataForDFNWorks = false;
					}
					else if (int_correctMesh != 0){
						dataForDFNWorks = true;
					}
					break;
				}
			}

			infile.close();
		}

		if(hasinverseSignTheta){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				//std::cout << "linenospace: " << linenospace <<std::endl;

				if (linenospace.substr (0,1) != "#" && linenospace.find("inverseSignTheta=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_inverseSignTheta;

					//std::cout << "int_inverseSignTheta = " << int_inverseSignTheta <<std::endl;

					if (int_inverseSignTheta == 0){
						inverseSignTheta = false;
					}
					else if (int_inverseSignTheta != 0){
						inverseSignTheta = true;
					}
					break;
				}
			}

			infile.close();
		}
		else
		{
			std::cout << "Warning : ParametersLaGriT block is set but the program failed to parse inverseSignTheta value. " << std::endl;
			std::cout << "The default value will be taken or please correct the input file (" << inputFile << ")." << std::endl;
			//exit(-1);
		}



		if(hasparamFilePrecision){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				//std::cout << "linenospace: " << linenospace <<std::endl;

				if (linenospace.substr (0,1) != "#" && linenospace.find("paramFilePrecision=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> int_paramFilePrecision;


					if (int_paramFilePrecision == 0){
						paramFilePrecision = false;
					}
					else if (int_paramFilePrecision != 0){
						paramFilePrecision = true;
					}
					break;
				}
			}

			infile.close();
		}
		else
		{
			std::cout << "Warning : ParametersLaGriT block is set but the program failed to parse paramFilePrecision value. " << std::endl;
			std::cout << "The default value will be taken or please correct the input file (" << inputFile << ")." << std::endl;
			//exit(-1);
		}

		if(hasdimScaling){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("dimScaling=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> dimScaling;

					break;
				}
			}
			infile.close();
		}

		if(hasintRound){
			infile.open(inputFile);
			while (!infile.eof())
			{
				std::getline(infile, line);
				linenospace = line;
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), '\t'), linenospace.end());
				linenospace.erase(remove(linenospace.begin(), linenospace.end(), ' '), linenospace.end());

				if (linenospace.substr (0,1) != "#" && linenospace.find("intRound=") != std::string::npos){
					unsigned first = line.find("=");
					unsigned last = line.size();
					if (linenospace.find("#") != std::string::npos){
						last = line.find("#")-1;
					}

					std::string strNew = line.substr (first+1,last-first);
					std::istringstream sstream(strNew);
					sstream >> intRound;

					break;
				}
			}
			infile.close();
		}
	}

	// Errors
	if(hasBlockDensity && ! hasn4){
		std::cerr << "Error : Density block is set but the program failed to parse n4 value. " << std::endl;
		std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
		exit(-1);
	}

	if(hasBlockSpatialParams){
		if(!hasDFNBB_LowerLeftX){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_LowerLeftX value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else if (!hasDFNBB_LowerLeftY){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_LowerLeftY value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else if (!hasDFNBB_LowerLeftZ){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_LowerLeftZ value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else if (!hasDFNBB_UpperRightX){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_UpperRightX value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else if (!hasDFNBB_UpperRightY){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_UpperRightY value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
		else if (!hasDFNBB_UpperRightZ){
			std::cerr << "Error : SpatialParams block is set but the program failed to parse DFNBB_UpperRightZ value. " << std::endl;
			std::cerr << "Please correct the input file (" << inputFile << ")." << std::endl;
			exit(-1);
		}
	}

	// Return block
	iparams.hasBlockDensity = hasBlockDensity;
	iparams.hasBlockSpatialParams = hasBlockSpatialParams;
	iparams.hasBlockAddingPointsMethod = hasBlockAddingPointsMethod;
	iparams.hasBlockDiscreMethod = hasBlockDiscreMethod;
	iparams.hasBlockCleanUpDirectory = hasBlockCleanUpDirectory;
	iparams.hasBlockDataFile = hasBlockDataFile;
	iparams.hasBlockParametersLaGriT = hasBlockParametersLaGriT;
	iparams.hasBlockRemoveClosePoints = hasBlockRemoveClosePoints;

	// Return parameter
	iparams.hasn4 = hasn4;
	iparams.hasDFNBB_LowerLeftX = hasDFNBB_LowerLeftX;
	iparams.hasDFNBB_LowerLeftY = hasDFNBB_LowerLeftY;
	iparams.hasDFNBB_LowerLeftZ = hasDFNBB_LowerLeftZ;
	iparams.hasDFNBB_UpperRightX = hasDFNBB_UpperRightX;
	iparams.hasDFNBB_UpperRightY = hasDFNBB_UpperRightY;
	iparams.hasDFNBB_UpperRightZ = hasDFNBB_UpperRightZ;
	iparams.hasaddingPointsMethod = hasaddingPointsMethod;
	iparams.hasdiscreMethod = hasdiscreMethod;
	iparams.hasdataFile = hasdataFile;
	iparams.hasdesactFracsFile = hasdesactFracsFile;
	iparams.hasEPS_INT	  = hasEPS_INT;
	iparams.hasEPS_FILTER = hasEPS_FILTER;
	iparams.hasEPS_BOUNDARY = hasEPS_BOUNDARY;
	iparams.hasremoveClosePoints = hasremoveClosePoints;
	iparams.hasmergeCloseIntersectionPoints = hasmergeCloseIntersectionPoints;
	iparams.hasmergeClosePointsRelCritLength = hasmergeClosePointsRelCritLength;
	iparams.hascorrectMesh = hascorrectMesh;
	iparams.hasdataForDFNWorks = hasdataForDFNWorks;
	iparams.hasparamFilePrecision = hasparamFilePrecision;
	iparams.hasinverseSignTheta = hasinverseSignTheta;
	iparams.hasminLengthRatio = hasminLengthRatio;
	iparams.hasdimScaling = hasdimScaling;
	iparams.hasintRound = hasintRound;

	iparams.n4 = n4;
	iparams.DFNBB_LowerLeftX = DFNBB_LowerLeftX;
	iparams.DFNBB_LowerLeftY = DFNBB_LowerLeftY;
	iparams.DFNBB_LowerLeftZ = DFNBB_LowerLeftZ;
	iparams.DFNBB_UpperRightX = DFNBB_UpperRightX;
	iparams.DFNBB_UpperRightY = DFNBB_UpperRightY;
	iparams.DFNBB_UpperRightZ = DFNBB_UpperRightZ;
	iparams.addingPointsMethod = addingPointsMethod;
	iparams.discreMethod = discreMethod;
	iparams.cleanupDirectory = cleanupDirectory;;
	iparams.dataFile = dataFile;
	iparams.desactFracsFile = desactFracsFile;
	iparams.EPS_INT = EPS_INT;
	iparams.EPS_FILTER = EPS_FILTER;
	iparams.EPS_BOUNDARY = EPS_BOUNDARY;
	iparams.rm_close_pts = rm_close_pts;
	iparams.merge_close_intersection_pts = merge_close_intersection_pts;
	iparams.mergeClosePointsRelCritLength = mergeClosePointsRelCritLength;
	iparams.minLengthRatio = minLengthRatio;
	iparams.correctMesh = correctMesh;
	iparams.dataForDFNWorks = dataForDFNWorks;
	iparams.inverseSignTheta = inverseSignTheta;
	iparams.paramFilePrecision = paramFilePrecision;
	iparams.dimScaling = dimScaling;
	iparams.intRound = intRound;

	return iparams;
}

} // end of namespace

#endif /* INCLUDE_PARSE_PARAMETER_FILE_HH_ */
