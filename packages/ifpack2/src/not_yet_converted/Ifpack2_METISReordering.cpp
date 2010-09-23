/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Reordering.hpp"
#include "Ifpack2_METISReordering.hpp"
#include "Ifpack2_Graph.hpp"
#include "Ifpack2_Graph_Tpetra_CrsGraph.hpp"
#include "Ifpack2_Graph_Tpetra_RowMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_ParameterList.hpp"

typedef int idxtype;
#ifdef HAVE_IFPACK2_METIS
extern "C" {
  void METIS_NodeND(int *n, idxtype *xadj, idxtype *adjncy, 
		    int *numflag, int *options, int *perm, int *iperm);
}
#endif

//==============================================================================
Ifpack2_METISReordering::Ifpack2_METISReordering() :
  UseSymmetricGraph_(false),
  NumMyRows_(0),
  IsComputed_(false)
{}

//==============================================================================
// Mainly copied from Ifpack2_METISPartitioner.cpp
//
// NOTE:
// - matrix is supposed to be localized, and passes through the
// singleton filter. This means that I do not have to look
// for Dirichlet nodes (singletons). Also, all rows and columns are 
// local.
int Ifpack2_METISReordering::Compute(const Ifpack2_Graph& Graph)
{

  NumMyRows_ = Graph.NumMyRows();
  Reorder_.resize(NumMyRows_);
  InvReorder_.resize(NumMyRows_);

  int ierr;

  Teuchos::RCP<Tpetra_CrsGraph> SymGraph;
  Teuchos::RCP<Tpetra_Map> SymMap;
  Teuchos::RCP<Ifpack2_Graph_Tpetra_CrsGraph> SymTIFPACKGraph;
  Teuchos::RCP<Ifpack2_Graph> TIFPACKGraph = Teuchos::rcp( (Ifpack2_Graph*)&Graph, false );

  int Length = 2 * Graph.MaxMyNumEntries();
  int NumIndices;
  vector<int> Indices;
  Indices.resize(Length);

  vector<int> options;
  options.resize(8);
  options[0] = 0; // default values

#ifdef HAVE_IFPACK2_METIS  
  int numflag = 0; // C style
#endif

  if (UseSymmetricGraph_) {

    // need to build a symmetric graph. 
    // I do this in two stages:
    // 1.- construct an Tpetra_CrsMatrix, symmetric
    // 2.- convert the Tpetra_CrsMatrix into METIS format
    SymMap = Teuchos::rcp( new Tpetra_Map(NumMyRows_,0,Graph.Comm()) );
    SymGraph = Teuchos::rcp( new Tpetra_CrsGraph(Copy,*SymMap,0) );

    for (int i = 0; i < NumMyRows_ ; ++i) {

      ierr = Graph.ExtractMyRowCopy(i, Length, NumIndices, 
				      &Indices[0]);
      IFPACK2_CHK_ERR(ierr);

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj != i) {
          // insert A(i,j), then A(j,i)
	  SymGraph->InsertGlobalIndices(i,1,&jj);
	  SymGraph->InsertGlobalIndices(jj,1,&i);
	}
      }      
    }
    IFPACK2_CHK_ERR(SymGraph->OptimizeStorage());
    IFPACK2_CHK_ERR(SymGraph->FillComplete());
    SymTIFPACKGraph = Teuchos::rcp( new Ifpack2_Graph_Tpetra_CrsGraph(SymGraph) );
    TIFPACKGraph = SymTIFPACKGraph;
  }

  // convert to METIS format
  vector<idxtype> xadj;
  xadj.resize(NumMyRows_ + 1);

  vector<idxtype> adjncy;
  adjncy.resize(Graph.NumMyNonzeros());
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows_ ; ++i) {

    xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */

    ierr = TIFPACKGraph->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
    IFPACK2_CHK_ERR(ierr);

    for (int j = 0 ; j < NumIndices ; ++j) {
      int jj = Indices[j];
      if (jj != i) {
	adjncy[count++] = jj;
	xadj[count2+1]++;
      }
    }
    count2++;
  }

#ifdef HAVE_IFPACK2_METIS
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
int Ifpack2_METISReordering::Compute(const Tpetra_RowMatrix& Matrix)
{
  Ifpack2_Graph_Tpetra_RowMatrix Graph(Teuchos::rcp(&Matrix, false));

  IFPACK2_CHK_ERR(Compute(Graph));

  return(0);
}

//==============================================================================
int Ifpack2_METISReordering::Reorder(const int i) const
{
#ifdef IFPACK2_ABC
  if (!IsComputed())
    IFPACK2_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK2_CHK_ERR(-1);
#endif

  return(Reorder_[i]);
}

//==============================================================================
int Ifpack2_METISReordering::InvReorder(const int i) const
{
#ifdef IFPACK2_ABC
  if (!IsComputed())
    IFPACK2_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK2_CHK_ERR(-1);
#endif

  return(InvReorder_[i]);
}
//==============================================================================
int Ifpack2_METISReordering::P(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xorig,
			    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const
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
int Ifpack2_METISReordering::Pinv(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xorig,
				 Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const
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
ostream& Ifpack2_METISReordering::Print(std::ostream& os) const
{
  os << "*** Ifpack2_METISReordering" << endl << endl;
  if (!IsComputed())
    os << "*** Reordering not yet computed." << endl;
  
  os << "*** Number of local rows = " << NumMyRows_ << endl;
  os << "Local Row\tReorder[i]\tInvReorder[i]" << endl;
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    os << '\t' << i << "\t\t" << Reorder_[i] << "\t\t" << InvReorder_[i] << endl;
  }
   
  return(os);
}

