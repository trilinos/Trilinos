//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_EpetraMatcher_hpp_
#define _Isorropia_EpetraMatcher_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#ifdef ISORROPIA_HAVE_OMP
#include <omp.h>
#endif
#include <ctime>
#include <stack>
#include <set>

namespace std{
class Node
{
public:
	int layer_num;
	int scanned;
	vector<int> edgelist;
};

class matcher {
public:
	//input is MTX format
	//graph is compressed row storage of the input
	//vlayered is a vector of v's already in any layer.
	//mateU and mateV are vectors which stores the matched vertices in U and V

	vector<Node> LU;
	vector<Node> LV;
	vector<int> vlayered;
	vector<vector<int> > del_m;
	vector<vector<int> > graph;
	vector<int> mateU;
	vector<int> mateV;
	vector<int> vlist;
	set<int> nvlist;
	int *CRS_Pointers;
	int *CRS_Indices;
	double *CRS_Vals;
	int U,V,E;
	int k_star;
	int icm;
	bool finish;

	/// Other Functions
	matcher(char *);
	virtual ~matcher();
	void process_mtx_compressed(char *);
	int vlayer_clear();
	int is_intersect(int);
	void delete_matched_v();
	bool remove_edge(int, int);
	void add_edge(int, int);

	// HK functions
	vector<int> get_matching();
	int construct_layered_graph();
	void find_set_del_M();
	int recursive_path_finder(int, int, vector<int>*);
	int iterative_path_finder(int, int, vector<int>*);
	int augment_matching();
	void update_vlayered(int);
	bool U0_empty();
};
} //namespace std
#endif

