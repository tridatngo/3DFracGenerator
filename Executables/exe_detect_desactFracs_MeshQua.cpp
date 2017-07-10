/*
 * exe_create_parameter_mlgi_files.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: ngotr
 */

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iomanip> // setprecision
#include <list>
#include <unistd.h>
#include <ctime>
#include <boost/utility.hpp>
#include <math.h>
#include <algorithm>
#include <iterator>     /* std::next & std::prev */
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <cstdlib>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/include/prior.hpp>
#include <boost/lexical_cast.hpp>

#define VERY_VERBOSE 0

int count_digits (const int &the_integer){
	int nDigits = std::floor(std::log10(std::abs(the_integer))) + 1;

	return nDigits;
}


std::string int2str_setw (const int &a,const int &w){
	std::string out;
	std::ostringstream out_ss;
	out_ss << std::setfill('0') << std::setw(w) << a;
	out = out_ss.str();

	return out;
}
//int main(int argc, char** argv) {

class CSVRow
{
public:
	std::string const& operator[](std::size_t index) const
	{
		return m_data[index];
	}
	std::size_t size() const
	{
		return m_data.size();
	}
	void readNextRow(std::istream& str)
	{
		std::string         line;
		std::getline(str, line);

		std::stringstream   lineStream(line);
		std::string         cell;

		m_data.clear();
		while(std::getline(lineStream, cell, ','))
		{
			m_data.push_back(cell);
		}
		// This checks for a trailing comma with no data after it.
		if (!lineStream && cell.empty())
		{
			// If there was a trailing comma then add an empty element.
			m_data.push_back("");
		}
	}
private:
	std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}

class CSVIterator
{
public:
	typedef std::input_iterator_tag     iterator_category;
	typedef CSVRow                      value_type;
	typedef std::size_t                 difference_type;
	typedef CSVRow*                     pointer;
	typedef CSVRow&                     reference;

	CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
	CSVIterator()                   :m_str(NULL) {}

	// Pre Increment
	CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
	// Post increment
	CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
	CSVRow const& operator*()   const       {return m_row;}
	CSVRow const* operator->()  const       {return &m_row;}

	bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
	bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
private:
	std::istream*       m_str;
	CSVRow              m_row;
};

void find_and_replace(std::string& source, std::string const& find, std::string const& replace)
{
	for(std::string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
	{
		source.replace(i, find.length(), replace);
		i += replace.length();
	}
}

int main(int argc, char** argv) {

	std::string filename_ = "";
	std::string cwd = "";
	double mq_thres(10);
	if( argc == 3 ) {
		filename_ = argv[1];
		mq_thres = boost::lexical_cast<double>(argv[2]);
	}
	else if ( argc == 4 ){
		filename_ = argv[1];
		mq_thres = boost::lexical_cast<double>(argv[2]);
		cwd = argv[3];
	}
	else {
		std::cout << "Usage: \"./exefile CSVInputFile Threshold\" or \"./exefile CSVInputFile Threshold CurrentWorkingDirectory\"\n";
		std::cerr << "Failed to open the input files." << std::endl;
		return -1;
	}

	std::ifstream       file(filename_.c_str());

	CSVRow              row;
	int num_lines(0), mat_ID_col(0), MeshQuality_col(0);
	std::list<int>	Material_ID;
	std::list<double> MeshQuality;

	while(file >> row)
	{
		// Read the header
		if(num_lines == 0){
			for (int i = 0; i!= row.size(); i++){
				//std::cout << row[i].c_str() << "\n";
				if (std::strcmp(row[i].c_str(), "\"Material_ID\"") == 0){
					mat_ID_col = i;
				}
				else if (std::strcmp(row[i].c_str(), "\"Quality\"") == 0){
					MeshQuality_col = i;
				}
			}
		}
		else{
			Material_ID.push_back(boost::lexical_cast<int>(row[mat_ID_col]));
			MeshQuality.push_back(boost::lexical_cast<double>(row[MeshQuality_col]));
		}
		num_lines++;
	}

	//std::cout << "Number of lines = " << num_lines << "\n";
	//std::cout << "Material_ID_col = " << mat_ID_col << "\n";
	//std::cout << "Quality_col = " << MeshQuality_col << "\n";
	std::cout << "Quality.size() = " << MeshQuality.size() << "\n";
	std::cout << "Material_ID.size() = " << Material_ID.size() << "\n";

	std::string line, linenospace;
	std::list<int> fracToRemoveList, oldfracParentNameList, finalfracParentNameList;

	// List of fractures which are to be removed
	int mq_iter(0);

	std::cout << "Creating the list of fractures that will be removed...\n";
	for (std::list<int>::iterator matID_iter = Material_ID.begin(); matID_iter != Material_ID.end(); matID_iter++){
		std::list<double>::iterator mq_iter_list(MeshQuality.begin());
		std::advance(mq_iter_list,mq_iter);
		if(*mq_iter_list > mq_thres){
			fracToRemoveList.push_back(*matID_iter);
		}
		//std::cout << "mq_iter = " << mq_iter << std::endl;
		mq_iter++;
	}

	std::cout << "Creating the list of fractures that will be removed... Complete.\n";

	std::string filename3(cwd + "/output/insidebbox_frac.txt");
	std::string filename4(cwd + "/desactFracsFile.dat");
	std::system("rm -rf filename4");

	//std::replace( filename3.begin(), filename3.end(), "..o", "../o");
	//std::replace( filename4.begin(), filename4.end(), "..d", "../d");
	//std::replace( filename3.begin(), filename3.end(), '\/\/', '\/');
	//std::replace( filename4.begin(), filename4.end(), '\/\/', '\/');

	find_and_replace( filename3, "..o", "../o");
	find_and_replace( filename4, "..d", "../d");
	find_and_replace( filename3, "//", "/");
	find_and_replace( filename4, "//", "/");

	std::ifstream    infile3(filename3.c_str());
	if (! infile3.is_open()) {
		std::cerr << "Error: Failed to open the input file \""<< filename3 << "\"."<< std::endl;
		exit(-1);
	}
	while(!infile3.eof()){
		std::getline(infile3, line);
		int dummy_int;
		std::istringstream ist(line);
		ist.str(line);
		ist.clear();
		ist >> dummy_int;
		oldfracParentNameList.push_back(dummy_int);
	}

	for (std::list<int>::iterator it = fracToRemoveList.begin(); it != fracToRemoveList.end(); it++){
		std::cout << "FracParentName = " << *it << std::endl;
		std::list<int>::iterator oldfracParent_iter(oldfracParentNameList.begin());
		std::advance(oldfracParent_iter, *it-1);
		finalfracParentNameList.push_back(*oldfracParent_iter);
	}
	finalfracParentNameList.sort();
	finalfracParentNameList.unique();

	infile3.close();

	if(finalfracParentNameList.size() == 0){
		std::cout << "All grid cells has an aspect ratio lower than " << mq_thres << "." << std::endl;
	}
	else{
		std::ofstream    outfile(filename4.c_str());
		for (std::list<int>::iterator it = finalfracParentNameList.begin(); it != finalfracParentNameList.end(); it++){
			// NTD 08Feb2017
			std::cout << "FinalFracParentName = " << *it << std::endl;
			outfile << *it << std::endl;
		}
	}

	std::cout << std::endl <<"Program terminated with success." << std::endl << std::endl;

	return 1;
}



