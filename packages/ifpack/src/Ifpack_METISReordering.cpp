/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Reordering.h"
#include "Ifpack_METISReordering.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

typedef int idxtype;
#ifdef HAVE_IFPACK_METIS
extern "C" {
  void METIS_NodeND(int *n, idxtype *xadj, idxtype *adjncy, 
		    int *numflag, int *options, int *perm, int *iperm);
}
#endif

//==============================================================================
Ifpack_METISReordering::Ifpack_METISReordering() :
  UseSymmetricGraph_(false),
  NumMyRows_(0),
  IsComputed_(false)
{}

//==============================================================================
// Mainly copied from Ifpack_METISPartitioner.cpp
//
// NOTE:
// - matrix is supposed to be localized, and passes through the
// singleton filter. This means that I do not have to look
// for Dirichlet nodes (singletons). Also, all rows and columns are 
// local.
int Ifpack_METISReordering::Compute(const Ifpack_Graph& Graph)
{

  NumMyRows_ = Graph.NumMyRows();
  Reorder_.resize(NumMyRows_);
  InvReorder_.resize(NumMyRows_);

  int ierr;

  Teuchos::RefCountPtr<Epetra_CrsGraph> SymGraph;
  Teuchos::RefCountPtr<Epetra_Map> SymMap;
  Teuchos::RefCountPtr<Ifpack_Graph_Epetra_CrsGraph> SymIFPACKGraph;
  Teuchos::RefCountPtr<Ifpack_Graph> IFPACKGraph = Teuchos::rcp( (Ifpack_Graph*)&Graph, false );

  int Length = 2 * Graph.MaxMyNumEntries();
  int NumIndices;
  std::vector<int> Indices;
  Indices.resize(Length);

  std::vector<int> options;
  options.resize(8);
  options[0] = 0; // default values

#ifdef HAVE_IFPACK_METIS  
  int numflag = 0; // C style
#endif

  if (UseSymmetricGraph_) {

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
    // need to build a symmetric graph. 
    // I do this in two stages:
    // 1.- construct an Epetra_CrsMatrix, symmetric
    // 2.- convert the Epetra_CrsMatrix into METIS format
    SymMap = Teuchos::rcp( new Epetra_Map(NumMyRows_,0,Graph.Comm()) );
    SymGraph = Teuchos::rcp( new Epetra_CrsGraph(Copy,*SymMap,0) );
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(SymGraph->RowMap().GlobalIndicesInt()) {
    for (int i = 0; i < NumMyRows_ ; ++i) {

      ierr = Graph.ExtractMyRowCopy(i, Length, NumIndices, 
				      &Indices[0]);
      IFPACK_CHK_ERR(ierr);

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj != i) {
          // insert A(i,j), then A(j,i)
	  SymGraph->InsertGlobalIndices(i,1,&jj);
	  SymGraph->InsertGlobalIndices(jj,1,&i);
	}
      }      
    }
	}
	else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(SymGraph->RowMap().GlobalIndicesLongLong()) {
    for (int i = 0; i < NumMyRows_ ; ++i) {
      long long i_LL = i;

      ierr = Graph.ExtractMyRowCopy(i, Length, NumIndices, 
				      &Indices[0]);
      IFPACK_CHK_ERR(ierr);

      for (int j = 0 ; j < NumIndices ; ++j) {
	long long jj = Indices[j];
	if (jj != i) {
          // insert A(i,j), then A(j,i)
	  SymGraph->InsertGlobalIndices(i_LL,1,&jj);
	  SymGraph->InsertGlobalIndices(jj,1,&i_LL);
	}
      }      
    }
    }
	else
#endif
    throw "Ifpack_METISReordering::Compute: GlobalIndices type unknown";

    IFPACK_CHK_ERR(SymGraph->OptimizeStorage());
    IFPACK_CHK_ERR(SymGraph->FillComplete());
    SymIFPACKGraph = Teuchos::rcp( new Ifpack_Graph_Epetra_CrsGraph(SymGraph) );
    IFPACKGraph = SymIFPACKGraph;
  }

  // convert to METIS format
  std::vector<idxtype> xadj;
  xadj.resize(NumMyRows_ + 1);

  std::vector<idxtype> adjncy;
  adjncy.resize(Graph.NumMyNonzeros());
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows_ ; ++i) {

    xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */

    ierr = IFPACKGraph->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int j = 0 ; j < NumIndices ; ++j) {
      int jj = Indices[j];
      if (jj != i) {
	adjncy[count++] = jj;
	xadj[count2+1]++;
      }
    }
    count2++;
  }

#ifdef HAVE_IFPACK_METIS
  // vectors from METIS. The second last is `perm', the last is `iperm'.
  // They store the fill-reducing permutation and inverse-permutation.
  // Let A be the original matrix and A' the permuted matrix. The
  // arrays perm and iperm are defined as follows. Row (column) i of A'
  // if the perm[i] row (col) of A, and row (column) i of A is the
  // iperm[i] row (column) of A'. The numbering starts from 0 in our case.
  METIS_NodeND(&NumMyRows_, &xadj[0], &adjncy[0],
	       &numflag, &options[0],
	       &InvReorder_[0], &Reorder_[0]);
#else
  cerr << "Please configure with --enable-ifpack-metis" << endl;
  cerr << "to use METIS Reordering." << endl;
  exit(EXIT_FAILURE);
#endif
      
  return(0);
} 

//==============================================================================
int Ifpack_METISReordering::Compute(const Epetra_RowMatrix& Matrix)
{
  Ifpack_Graph_Epetra_RowMatrix Graph(Teuchos::rcp(&Matrix, false));

  IFPACK_CHK_ERR(Compute(Graph));

  return(0);
}

//==============================================================================
int Ifpack_METISReordering::Reorder(const int i) const
{
#ifdef IFPACK_ABC
  if (!IsComputed())
    IFPACK_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK_CHK_ERR(-1);
#endif

  return(Reorder_[i]);
}

//==============================================================================
int Ifpack_METISReordering::InvReorder(const int i) const
{
#ifdef IFPACK_ABC
  if (!IsComputed())
    IFPACK_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK_CHK_ERR(-1);
#endif

  return(InvReorder_[i]);
}
//==============================================================================
int Ifpack_METISReordering::P(const Epetra_MultiVector& Xorig,
			    Epetra_MultiVector& X) const
{  
  int NumVectors = X.NumVectors();

  for (int j = 0 ; j < NumVectors ; ++j) {
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      int np = Reorder_[i];
      X[j][np] = Xorig[j][i];
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_METISReordering::Pinv(const Epetra_MultiVector& Xorig,
				 Epetra_MultiVector& X) const
{
  int NumVectors = X.NumVectors();

  for (int j = 0 ; j < NumVectors ; ++j) {
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      int np = Reorder_[i];
      X[j][i] = Xorig[j][np];
    }
  }

  return(0);
}

//==============================================================================
ostream& Ifpack_METISReordering::Print(std::ostream& os) const
{
  os << "*** Ifpack_METISReordering" << endl << endl;
  if (!IsComputed())
    os << "*** Reordering not yet computed." << endl;
  
  os << "*** Number of local rows = " << NumMyRows_ << endl;
  os << "Local Row\tReorder[i]\tInvReorder[i]" << endl;
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    os << '\t' << i << "\t\t" << Reorder_[i] << "\t\t" << InvReorder_[i] << endl;
  }
   
  return(os);
}

