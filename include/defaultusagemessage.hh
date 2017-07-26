/*
 * defaultusagemessage.hh
 *
 *  Created on: Jul 26, 2017
 *      Author: ngotr
 */

#ifndef INCLUDE_DEFAULTUSAGEMESSAGE_HH_
#define INCLUDE_DEFAULTUSAGEMESSAGE_HH_

#include <string>

namespace CGAL
{

inline std::string defaultUsageMessage(const std::string& programName)
{
	return  "Usage: " + programName + " [options] \n"
			"Options usually are parameters given to the simulation, \n"
			"and have to be specified with this syntax: \n"
			"\t-GroupName.ParameterName VALUE, for example -dataFile ellipse_input.dat\n"
			"\n"
			"If -ParameterFile is specified, parameters can also be defined there. In this case,\n"
			"lines of the form \n"
			"\n"
			"GroupName.dataFile = VALUE # comment \n"
			"\n"
			"have to be used. More conveniently, group names can be specified in square brackets, \n"
			"such that each following parameter name belongs to that group, \n"
			"\n"
			"[GroupName] \n"
			"dataFile = VALUE \n"
			"\n"
			"See dfnGen.input for an example of a parameter file. \n"
			"\n"
			"Parameters specified on the command line have priority over those in the parameter file.\n"
			"If no parameter file name is given, 'dfnGen.input' is chosen as default.\n"
			"\n"
			"Important options include:\n"
			"\t-h, --help                        Print this usage message and exit\n"
			"\t-PrintParameters [true|false]     Print the run-time modifiable parameters _after_ \n"
			"\t                                  the simulation [default: true]\n"
			"\t-PrintProperties [true|false]     Print the compile-time parameters _before_ \n"
			"\t                                  the simulation [default: false]\n"
			"\t-dataFile FILENAME                File with parameter definitions\n"
			"\t-np NumProc                       Number of processes\n"
			"\n\n";
}
} // end namespace CGAL

#endif /* INCLUDE_DEFAULTUSAGEMESSAGE_HH_ */
