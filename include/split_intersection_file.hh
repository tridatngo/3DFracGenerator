/*
 * split_intersection_file.hh
 *
 *  Created on: Sep 21, 2016
 *      Author: ngotr
 */

#ifndef INCLUDE_SPLIT_INTERSECTION_FILE_HH_
#define INCLUDE_SPLIT_INTERSECTION_FILE_HH_

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
#include<list>
#include "myglobal_functions.hh"

namespace CGAL{
template <typename NT>
void split_intersection_file(std::string &fname, int & nEll){

	typedef typename std::list<NT>::iterator NT_Iterator;

	std::cout << "Writing intersection files (*.inp) for the ellipses. " << std::endl;
	int digits(CGAL::count_digits(nEll));

	for (int it = 0; it != nEll; it++) {
		NT x, y, z;
		std::list<NT> vecx, vecy, vecz;
		std::list<int> vecEll2;

		int nEll1, nEll2, dummy, num_lines(0);

		std::ifstream    infile(fname.c_str());
		if (! infile.is_open()) {
			std::cerr << "Failed to open the input file " << fname << "." << std::endl;
			exit(-1);
		}

		std::string line;
		std::getline(infile, line);

		//infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			std::istringstream ist(line);
			ist >> x >> y >> z >> nEll1 >> nEll2;
			if (nEll1 == it) {
				vecx.push_back(x);
				vecy.push_back(y);
				vecz.push_back(z);
				vecEll2.push_back(nEll2);
			}
		}
		infile.close();

		//std::cout << "vecx.size() = " << vecx.size() << std::endl;

		std::string		filename = "intersections/intersections_parent_" + CGAL::int2str_setw(it + 1, digits) + ".inp";
		std::ofstream   outfile(filename.c_str());

		if (vecx.size() > 0){
			// std::cout << "vecx.size() = " << vecx.size() << std::endl;
			for (int j = 0; j != vecx.size()-1; j++) {
				std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
				std::list<int>::iterator Ell_next_iter = vecEll2.begin();

				std::advance(Ell_curr_iter, j);
				std::advance(Ell_next_iter, j+1);

				if (*Ell_curr_iter == *Ell_next_iter) {
					num_lines++;
				}
			}

			outfile << vecx.size() << " " << num_lines << " 2 0 0\n";

			for (int j = 0; j != vecx.size(); j++) {
				NT_Iterator xj_iter = vecx.begin();
				std::advance(xj_iter, j);

				NT_Iterator yj_iter = vecy.begin();
				std::advance(yj_iter, j);

				NT_Iterator zj_iter = vecz.begin();
				std::advance(zj_iter, j);

				outfile << j + 1 << " " << std::setprecision(20) << *xj_iter << " " << *yj_iter << " " << *zj_iter << "\n";
			}

			int elem_cnt = 1;
			for (int j = 0; j != vecx.size() - 1; j++) {
				std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
				std::list<int>::iterator Ell_next_iter = vecEll2.begin();

				std::advance(Ell_curr_iter, j);
				std::advance(Ell_next_iter, j + 1);

				if (*Ell_curr_iter == *Ell_next_iter) {
					outfile << elem_cnt << " " << it + 1 << " line " << j + 1 << " " << j + 2 << "\n";
					elem_cnt += 1;
				}
			}


			outfile << "2 1 1\n";
			outfile << "a_b, integer\n";
			outfile << "b_a, integer\n";

			for (int j = 0; j != vecx.size(); j++) {
				std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
				std::advance(Ell_curr_iter, j) ;
				outfile << j + 1 << " " << it << " " << *Ell_curr_iter << "\n";
			}
			vecx.clear();
			vecy.clear();
			vecz.clear();
			vecEll2.clear();

			//std::cout << "Writing the intersection *.inp file for the ellipse " << it <<" : Complete." << std::endl;
			outfile.close();
		}
		else{
			outfile.close();
		}

	}

	std::cout << "Writing intersection files *.inp for ellipses: Complete." << std::endl;
}

// TODO : implement splitting intersection file with sorting point
template <typename NT, typename PT>
void split_intersection_file_sort(std::string &fname, int & nEll){

	typedef typename std::list<NT>::iterator NT_Iterator;
	typedef typename std::list<PT>			 PointList;
	typedef typename std::list<PT>::iterator Pts_Iterator;

	int digits(CGAL::count_digits(nEll));

	NT x, y, z;
	std::list<NT> vecx, vecy, vecz;
	std::list<int> vecEll2;

	int nEll1, nEll2, dummy, num_lines(0);

	for (int it = 0; it != nEll; it++) {
		std::ifstream    infile(fname.c_str());
		if (! infile.is_open()) {
			std::cerr << "Failed to open the input file " << fname << "." << std::endl;
			exit(-1);
		}

		std::string line;
		std::getline(infile, line);

		//infile.open(inputFile);
		while (!infile.eof())
		{
			std::getline(infile, line);
			std::istringstream ist(line);
			ist >> x >> y >> z >> nEll1 >> nEll2;
			if (nEll1 == it) {
				vecx.push_back(x);
				vecy.push_back(y);
				vecz.push_back(z);
				vecEll2.push_back(nEll2);
			}
		}
		infile.close();

		std::string		filename = "intersections/intersections_parent_" + CGAL::int2str_setw(it + 1, digits) + ".inp";
		std::ofstream   outfile(filename.c_str());

		for (int j = 0; j != vecx.size()-1; j++) {
			std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
			std::list<int>::iterator Ell_next_iter = vecEll2.begin();

			std::advance(Ell_curr_iter, j);
			std::advance(Ell_next_iter, j+1);

			if (*Ell_curr_iter == *Ell_next_iter) {
				num_lines++;
			}
		}

		outfile << vecx.size() << " " << num_lines << " 2 0 0\n";

		for (int j = 0; j != vecx.size(); j++) {
			NT_Iterator xj_iter = vecx.begin();
			std::advance(xj_iter, j);

			NT_Iterator yj_iter = vecy.begin();
			std::advance(yj_iter, j);

			NT_Iterator zj_iter = vecz.begin();
			std::advance(zj_iter, j);

			outfile << j + 1 << " " << std::setprecision(20) << *xj_iter << " " << *yj_iter << " " << *zj_iter << "\n";
		}

		int elem_cnt = 1;
		for (int j = 0; j != vecx.size() - 1; j++) {
			std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
			std::list<int>::iterator Ell_next_iter = vecEll2.begin();

			std::advance(Ell_curr_iter, j);
			std::advance(Ell_next_iter, j + 1);

			if (*Ell_curr_iter == *Ell_next_iter) {
				outfile << elem_cnt << " " << it + 1 << " line " << j + 1 << " " << j + 2 << "\n";
				elem_cnt += 1;
			}
		}


		outfile << "2 1 1\n";
		outfile << "a_b, integer\n";
		outfile << "b_a, integer\n";

		for (int j = 0; j != vecx.size(); j++) {
			std::list<int>::iterator Ell_curr_iter = vecEll2.begin();
			std::advance(Ell_curr_iter, j) ;
			outfile << j + 1 << " " << it << " " << *Ell_curr_iter << "\n";
		}
		vecx.clear();
		vecy.clear();
		vecz.clear();
		vecEll2.clear();

		//std::cout << "Writing the intersection *.inp file for the ellipse " << it <<" : Complete." << std::endl;
		//outfile.close();
	}

	std::cout << "Writing intersection files (*.inp) for the ellipses: Complete." << std::endl;
}

} // end of namespace

#endif /* INCLUDE_SPLIT_INTERSECTION_FILE_HH_ */

/*
def dump_intersection_avs(i, x, y, z, ii):

        filename = 'intersections_'+str(i).zfill(digits)+'.inp'
        f = open(filename,'w+')
        num_lines = 0
        for j in range(len(x) - 1):
                if ii[j] == ii[j+1]:
                        num_lines += 1

        f.write('%d %d 2 0 0\n'%(len(x),num_lines))
        for j in range(len(x)):
                f.write('%d %f %f %f\n'%(j+1, x[j], y[j], z[j]))
        elem_cnt = 1
        for j in range(len(x) - 1):
                if ii[j] == ii[j+1]:
                        f.write('%d %d line %d %d\n'%(elem_cnt,i,j+1,j+2))
                        elem_cnt += 1
        f.write('2 1 1\n')
        f.write('a_b, integer\n')
        f.write('b_a, integer\n')
        for j in range(len(x)):
                f.write('%d %d %d\n'%(j+1, i, ii[j]))
        f.close()

        cmd = 'mv ' + filename + ' intersections/'
        os.system(cmd)

def split_intersection_file(i):

        data = genfromtxt(intersection_file, skip_header = 1)
        # Find lines in the data for fracture i
        locations = where(data[:,3] == i)[0][:]
        # Collect the rest of that data
        x = data[locations,0]
        y = data[locations,1]
        z = data[locations,2]
        ii = data[locations,4]
        dump_intersection_avs(i, x, y, z, ii)

def split_intersection_file_header(intersection_file, nCPU, nPoly):
        print "Splitting the Intersection File:", intersection_file

        os.system('rm -rf intersections')
        os.mkdir('intersections')
        jobnames = range(1, nPoly + 1)
        pool = mp.Pool(nCPU)
        pool.map(split_intersection_file, jobnames)
        pool.close()
        pool.join()
        pool.terminate()

        print "splitting the intersection file: complete\n"
 */
