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
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

// may need to change this for wierd installations
typedef int idxtype;
#ifdef HAVE_IFPACK_METIS
extern "C" {
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
int Ifpack_METISPartitioner::ComputePartitions()
{
  using std::cerr;
  using std::endl;

  int ierr;
#ifdef HAVE_IFPACK_METIS
  int edgecut;
#endif

  Teuchos::RefCountPtr<Epetra_CrsGraph> SymGraph ;
  Teuchos::RefCountPtr<Epetra_Map> SymMap;
  Teuchos::RefCountPtr<Ifpack_Graph_Epetra_CrsGraph> SymIFPACKGraph;
  Teuchos::RefCountPtr<Ifpack_Graph> IFPACKGraph = Teuchos::rcp( (Ifpack_Graph*)Graph_, false );

  int Length = 2 * MaxNumEntries();
  int NumIndices;
  std::vector<int> Indices;
  Indices.resize(Length);

  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  std::vector<idxtype> wgtflag;
  wgtflag.resize(4);

  std::vector<int> options;
  options.resize(4);

  int numflag;

  if (UseSymmetricGraph_) {

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
    // need to build a symmetric graph.
    // I do this in two stages:
    // 1.- construct an Epetra_CrsMatrix, symmetric
    // 2.- convert the Epetra_CrsMatrix into METIS format
    SymMap = Teuchos::rcp( new Epetra_Map(NumMyRows(),0,Graph_->Comm()) );
    SymGraph = Teuchos::rcp( new Epetra_CrsGraph(Copy,*SymMap,0) );
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      if(SymGraph->RowMap().GlobalIndicesInt()) {
        for (int i = 0; i < NumMyRows() ; ++i) {

          ierr = Graph_->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
          IFPACK_CHK_ERR(ierr);

          for (int j = 0 ; j < NumIndices ; ++j) {
            int jj = Indices[j];
            if (jj != i) {
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
        for (int i = 0; i < NumMyRows() ; ++i) {
          long long i_LL = i;

          ierr = Graph_->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
          IFPACK_CHK_ERR(ierr);

          for (int j = 0 ; j < NumIndices ; ++j) {
            long long jj = Indices[j];
            if (jj != i_LL) {
              SymGraph->InsertGlobalIndices(i_LL,1,&jj);
              SymGraph->InsertGlobalIndices(jj,1,&i_LL);
            }
          }
        }
      }
      else
#endif
        throw "Ifpack_METISPartitioner::ComputePartitions: GlobalIndices type unknown";

    IFPACK_CHK_ERR(SymGraph->FillComplete());
    SymIFPACKGraph = Teuchos::rcp( new Ifpack_Graph_Epetra_CrsGraph(SymGraph) );
    IFPACKGraph = SymIFPACKGraph;
  }

  // now work on IFPACKGraph, that can be the symmetric or
  // the non-symmetric one

  /* set parameters */

  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */

  std::vector<idxtype> xadj;
  xadj.resize(NumMyRows() + 1);

  std::vector<idxtype> adjncy;
  adjncy.resize(NumMyNonzeros());

  int count = 0;
  int count2 = 0;
  xadj[0] = 0;

  for (int i = 0; i < NumMyRows() ; ++i) {

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

  std::vector<idxtype> NodesInSubgraph;
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

#ifdef HAVE_IFPACK_METIS
      int j = NumMyRows();
      if (NumLocalParts_ < 8) {

        numflag = 0;

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
        IFPACK_CHK_ERR(-10); // something went wrong
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
