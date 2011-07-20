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

#include <Isorropia_Epetra.hpp>

#ifdef HAVE_EPETRA
//#ifdef HAVE_MPI
//#include <Epetra_MpiComm.h>
//#else
#include <Epetra_SerialComm.h>
//#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#endif
#ifdef ISORROPIA_HAVE_OMP
#include <omp.h>
#endif
#endif

#include <Isorropia_ConfigDefs.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <ctime>
#include <stack>
#include <set>

#ifdef MIN
#undef MIN
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))

namespace std{
class Isorropia_EpetraMatcher {	
private:
	int* mateU_;
	int* mateV_;
	int* Queue_;
	int* LV_;
	int* LU_;
	int* unmatchedU_;
	int *parent_;
	int *lookahead_;
	int *CRS_pointers_;
	int *CRS_indices_;
	double *CRS_vals_;
	bool finish_;
	int U_,V_,E_,avgDegU_,k_star_,icm_,BFSInd_,numThread_,Qst_,Qend_,matched_,choice_;
	
	#ifdef ISORROPIA_HAVE_OMP
	omp_lock_t *scannedV_;
	#endif
	

public:
	/// Interface Functions
	Isorropia_EpetraMatcher(const Epetra_CrsMatrix*, const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
	Isorropia_EpetraMatcher(Teuchos::RCP<const Epetra_CrsMatrix>,const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
	Isorropia_EpetraMatcher(const Epetra_CrsGraph *,const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
	Isorropia_EpetraMatcher(Teuchos::RCP<const Epetra_CrsGraph>,const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
	void extractRowPermutationCopy(int, int&, int* ) const;
	void extractColumnPermutationCopy(int, int&, int*) const;
	void getMatchedEdges(int,int&,int*) const;
	int getNumberOfMatchedVertices();
	Epetra_Map* getPermutedRowMap();
	Epetra_Map* getPermutedColumnMap();
	
	virtual ~Isorropia_EpetraMatcher();
	
	//Matching Functions
	void delete_matched_v();
	void filler();
	int match();
	int match_dfs();
	int match_hk();
	int construct_layered_graph();
	int find_set_del_M();
	int recursive_path_finder(int, int);
	int dfs_path_finder(int);
	int dfs_augment();
	int augment_matching(int);
	int DW_phase();
};
} //namespace std
#endif

