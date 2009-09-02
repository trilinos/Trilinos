/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Partitioner.hpp"
#include "Tifpack_OverlappingPartitioner.hpp"
#include "Tifpack_METISPartitioner.hpp"
#include "Tifpack_Graph.hpp"
#include "Tifpack_Graph_Tpetra_CrsGraph.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

// may need to change this for wierd installations
typedef int idxtype;
#ifdef HAVE_TIFPACK_METIS
extern "C" {
  void METIS_EstimateMemory(int *, idxtype *, idxtype *, int *, int *, int *);
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, 
			   idxtype *, int *, int *, int *, int *, int *,
			   idxtype *);
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, 
				idxtype *, idxtype *, int *, int *, int *, 
				int *, int *, idxtype *);

}
#endif

//==============================================================================
// NOTE:
// - matrix is supposed to be localized, and passes through the
// singleton filter. This means that I do not have to look
// for Dirichlet nodes (singletons). Also, all rows and columns are 
// local.
int Tifpack_METISPartitioner::ComputePartitions()
{

  int ierr;
#ifdef HAVE_TIFPACK_METIS
  int nbytes = 0;
  int edgecut;
#endif

  Teuchos::RefCountPtr<Tpetra_CrsGraph> SymGraph ;
  Teuchos::RefCountPtr<Tpetra_Map> SymMap;
  Teuchos::RefCountPtr<Tifpack_Graph_Tpetra_CrsGraph> SymTIFPACKGraph;
  Teuchos::RefCountPtr<Tifpack_Graph> TIFPACKGraph = Teuchos::rcp( (Tifpack_Graph*)Graph_, false );

  int Length = 2 * MaxNumEntries();
  int NumIndices;
  vector<int> Indices;
  Indices.resize(Length);

  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  vector<idxtype> wgtflag;
  wgtflag.resize(4);

  vector<int> options;
  options.resize(4);
  
  int numflag;

  if (UseSymmetricGraph_) {

    // need to build a symmetric graph. 
    // I do this in two stages:
    // 1.- construct an Tpetra_CrsMatrix, symmetric
    // 2.- convert the Tpetra_CrsMatrix into METIS format
    SymMap = Teuchos::rcp( new Tpetra_Map(NumMyRows(),0,Graph_->Comm()) );
    SymGraph = Teuchos::rcp( new Tpetra_CrsGraph(Copy,*SymMap,0) );

    for (int i = 0; i < NumMyRows() ; ++i) {

      ierr = Graph_->ExtractMyRowCopy(i, Length, NumIndices, 
				      &Indices[0]);
      TIFPACK_CHK_ERR(ierr);

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj != i) {
	  SymGraph->InsertGlobalIndices(i,1,&jj);
	  SymGraph->InsertGlobalIndices(jj,1,&i);
	}
      }      
    }
    TIFPACK_CHK_ERR(SymGraph->FillComplete());
    SymTIFPACKGraph = Teuchos::rcp( new Tifpack_Graph_Tpetra_CrsGraph(SymGraph) );
    TIFPACKGraph = SymTIFPACKGraph;
  }

  // now work on TIFPACKGraph, that can be the symmetric or
  // the non-symmetric one

  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  vector<idxtype> xadj;
  xadj.resize(NumMyRows() + 1);

  vector<idxtype> adjncy;
  adjncy.resize(NumMyNonzeros());
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows() ; ++i) {

    xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */

    ierr = TIFPACKGraph->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
    TIFPACK_CHK_ERR(ierr);

    for (int j = 0 ; j < NumIndices ; ++j) {
      int jj = Indices[j];
      if (jj != i) {
	adjncy[count++] = jj;
	xadj[count2+1]++;
      }
    }
    count2++;
  }

  vector<idxtype> NodesInSubgraph;
  NodesInSubgraph.resize(NumLocalParts_);

  // some cases can be handled separately
  
  int ok;

  if (NumLocalParts() == 1) {

    for (int i = 0 ; i < NumMyRows() ; ++i) 
      Partition_[i] = 0;
    
  } else if (NumLocalParts() == NumMyRows()) {

    for (int i = 0 ; i < NumMyRows() ; ++i) 
      Partition_[i] = i;
  
  } else {

    ok = 0;

    // sometimes METIS creates less partitions than specified.
    // ok will check this problem, and recall metis, asking
    // for NumLocalParts_/2 partitions
    while (ok == 0) {
      
      for (int i = 0 ; i < NumMyRows() ; ++i) 
	Partition_[i] = -1;
    
#ifdef HAVE_TIFPACK_METIS
      int j = NumMyRows();
      if (NumLocalParts_ < 8) {

	int i = 1; /* optype in the METIS manual */
	numflag = 0;
	METIS_EstimateMemory(&j, &xadj[0], &adjncy[0], 
			     &numflag, &i, &nbytes );
	
	METIS_PartGraphRecursive(&j, &xadj[0], &adjncy[0],
				 NULL, NULL,
				 &wgtflag[0], &numflag, &NumLocalParts_, 
				 &options[0], &edgecut, &Partition_[0]);
      } else {

	numflag = 0;
	
	METIS_PartGraphKway (&j, &xadj[0], &adjncy[0], 
			     NULL, 
			     NULL, &wgtflag[0], &numflag, 
			     &NumLocalParts_, &options[0],
			     &edgecut, &Partition_[0]);
      }
#else
      numflag = numflag * 2; // avoid warning for unused variable
      if (Graph_->Comm().MyPID() == 0) {
	cerr << "METIS was not linked; now I put all" << endl;
	cerr << "the local nodes in the same partition." << endl;
      }
      for (int i = 0 ; i < NumMyRows() ; ++i) 
	Partition_[i] = 0;
      NumLocalParts_ = 1;
#endif
      
      ok = 1;
      
      for (int i = 0 ; i < NumLocalParts() ; ++i) 
	NodesInSubgraph[i] = 0;

      for (int i = 0 ; i < NumMyRows() ; ++i) {
	int j = Partition_[i];
	if ((j < 0) || (j>= NumLocalParts())) {
	  ok = 0;
	  break;
	} 
	else NodesInSubgraph[j]++;
      }
      
      for (int i = 0 ; i < NumLocalParts() ; ++i) {
	if( NodesInSubgraph[i] == 0 ) {
	  ok = 0;
	  break;
	}
      }
      
      if (ok == 0) {
	cerr << "Specified number of subgraphs ("
	     << NumLocalParts_ << ") generates empty subgraphs." << endl;
	cerr << "Now I recall METIS with NumLocalParts_ = "
	     << NumLocalParts_ / 2 << "..." << endl;
	NumLocalParts_ = NumLocalParts_/2;
      }
      
      if (NumLocalParts() == 0) {
	TIFPACK_CHK_ERR(-10); // something went wrong
      }
      
      if (NumLocalParts() == 1) {
	for (int i = 0 ; i < NumMyRows() ; ++i) 
	  Partition_[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( NumLocalParts_ == 1 ) */

  return(0);
} 
