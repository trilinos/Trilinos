//@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#include "Ifpack_IlukGraph.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"

#ifdef HAVE_IFPACK_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>
#endif

//==============================================================================
Ifpack_IlukGraph::Ifpack_IlukGraph(const Epetra_CrsGraph & Graph, int LevelFill, int LevelOverlap)
  : Graph_(Graph),
    DomainMap_(Graph.DomainMap()),
    RangeMap_(Graph.RangeMap()),
    Comm_(Graph.Comm()),
    OverlapGraph_(0),
    OverlapRowMap_(0),
    OverlapImporter_(0),
    LevelFill_(LevelFill),
    LevelOverlap_(LevelOverlap),
    L_Graph_(0),
    U_Graph_(0),
    IndexBase_(Graph.IndexBase()),
    NumGlobalRows_(Graph.NumGlobalRows()),
    NumGlobalCols_(Graph.NumGlobalCols()),
    NumGlobalBlockRows_(Graph.NumGlobalBlockRows()),
    NumGlobalBlockCols_(Graph.NumGlobalBlockCols()),
    NumGlobalBlockDiagonals_(0),
    NumGlobalNonzeros_(0),
    NumGlobalEntries_(0),
    NumMyBlockRows_(Graph.NumMyBlockRows()),
    NumMyBlockCols_(Graph.NumMyBlockCols()),
    NumMyRows_(Graph.NumMyRows()),
    NumMyCols_(Graph.NumMyCols()),
    NumMyBlockDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyEntries_(0)
{
}

//==============================================================================
Ifpack_IlukGraph::Ifpack_IlukGraph(const Ifpack_IlukGraph & Graph) 
  : Graph_(Graph.Graph_),
    DomainMap_(Graph.DomainMap()),
    RangeMap_(Graph.RangeMap()),
    Comm_(Graph.Comm()),
    OverlapGraph_(Graph.OverlapGraph_),
    OverlapRowMap_(Graph.OverlapRowMap_),
    OverlapImporter_(Graph.OverlapImporter_),
    LevelFill_(Graph.LevelFill_),
    LevelOverlap_(Graph.LevelOverlap_),
    L_Graph_(0),
    U_Graph_(0),
    IndexBase_(Graph.IndexBase_),
    NumGlobalRows_(Graph.NumGlobalRows_),
    NumGlobalCols_(Graph.NumGlobalCols_),
    NumGlobalBlockRows_(Graph.NumGlobalBlockRows_),
    NumGlobalBlockCols_(Graph.NumGlobalBlockCols_),
    NumGlobalBlockDiagonals_(Graph.NumGlobalBlockDiagonals_),
    NumGlobalNonzeros_(Graph.NumGlobalNonzeros_),
    NumGlobalEntries_(Graph.NumGlobalEntries_),
    NumMyBlockRows_(Graph.NumMyBlockRows_),
    NumMyBlockCols_(Graph.NumMyBlockCols_),
    NumMyRows_(Graph.NumMyRows_),
    NumMyCols_(Graph.NumMyCols_),
    NumMyBlockDiagonals_(Graph.NumMyBlockDiagonals_),
    NumMyNonzeros_(Graph.NumMyNonzeros_),
    NumMyEntries_(Graph.NumMyEntries_)
{
  Epetra_CrsGraph & L_Graph_In = Graph.L_Graph();
  Epetra_CrsGraph & U_Graph_In = Graph.U_Graph();
  L_Graph_ = new Epetra_CrsGraph(L_Graph_In);
  U_Graph_ = new Epetra_CrsGraph(U_Graph_In);
}

//==============================================================================
Ifpack_IlukGraph::~Ifpack_IlukGraph()
{
  delete L_Graph_;
  delete U_Graph_;
  if (OverlapGraph_!=&Graph_) delete OverlapGraph_;
  if (OverlapRowMap_!=&Graph_.RowMap()) delete OverlapRowMap_;
  if (OverlapImporter_!=0) delete OverlapImporter_;
}

#ifdef HAVE_IFPACK_TEUCHOS
//==============================================================================
int Ifpack_IlukGraph::SetParameters(const Teuchos::ParameterList& parameterlist,
                                    bool cerr_warning_if_unused)
{
  Ifpack::param_struct params;
  params.int_params[Ifpack::level_fill-FIRST_INT_PARAM] = LevelFill_;
  params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM] = LevelOverlap_;

  Ifpack::set_parameters(parameterlist, params, cerr_warning_if_unused);

  LevelFill_ = params.int_params[Ifpack::level_fill-FIRST_INT_PARAM];
  LevelOverlap_ = params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM];
  return(0);
}
#endif

//==============================================================================
int Ifpack_IlukGraph::ConstructOverlapGraph() {

  OverlapGraph_ = (Epetra_CrsGraph *) &Graph_;
  OverlapRowMap_ = (Epetra_BlockMap *) &Graph_.RowMap();
  
  if (LevelOverlap_==0 || !Graph_.DomainMap().DistributedGlobal()) return(0); // Nothing to do

  Epetra_CrsGraph * OldGraph;
  Epetra_BlockMap * OldRowMap;
  Epetra_BlockMap * DomainMap = (Epetra_BlockMap *) &Graph_.DomainMap();
  Epetra_BlockMap * RangeMap = (Epetra_BlockMap *) &Graph_.RangeMap();
  for (int level=1; level <= LevelOverlap_; level++) {
    OldGraph = OverlapGraph_; 
    OldRowMap = OverlapRowMap_;

    OverlapImporter_ = (Epetra_Import *) OldGraph->Importer();
    OverlapRowMap_ = new Epetra_BlockMap(OverlapImporter_->TargetMap());

    
    if (level<LevelOverlap_)
      OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0);
    else
      // On last iteration, we want to filter out all columns except those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, *OverlapRowMap_, 0);

    EPETRA_CHK_ERR(OverlapGraph_->Import( Graph_, *OverlapImporter_, Insert));
    if (level<LevelOverlap_) {
      EPETRA_CHK_ERR(OverlapGraph_->TransformToLocal(DomainMap, RangeMap));
    }
    else {
      // Copy last OverlapImporter because we will use it later
      OverlapImporter_ = new Epetra_Import(*OverlapRowMap_, *DomainMap);
      EPETRA_CHK_ERR(OverlapGraph_->TransformToLocal(DomainMap, RangeMap));
    }

    if (OldGraph!=&Graph_) delete OldGraph;
    if (OldRowMap!=&Graph_.RowMap()) delete OldRowMap;
  }

    NumMyBlockRows_ = OverlapGraph_->NumMyBlockRows();
    NumMyBlockCols_ = OverlapGraph_->NumMyBlockCols();
    NumMyRows_ = OverlapGraph_->NumMyRows();
    NumMyCols_ = OverlapGraph_->NumMyCols();

  return(0);
}

//==============================================================================
int Ifpack_IlukGraph::ConstructFilledGraph() {
  int ierr = 0;
  int i, j;
  int * In=0, * L=0, * U = 0;
  int NumIn, NumL, NumU;
  bool DiagFound;

  
  EPETRA_CHK_ERR(ConstructOverlapGraph());

  L_Graph_ = new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0);
  U_Graph_ = new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0);


  // Get Maximun Row length
  int MaxNumIndices = OverlapGraph_->MaxNumIndices();

  L = new int[MaxNumIndices];
  U = new int[MaxNumIndices];
    

  // First we copy the user's graph into L and U, regardless of fill level

  for (i=0; i< NumMyBlockRows_; i++) {


    OverlapGraph_->ExtractMyRowView(i, NumIn, In); // Get Indices

    
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = In[j];

      if (k<NumMyBlockRows_) { // Ignore column elements that are not in the square matrix

	if (k==i) DiagFound = true;

	else if (k < i) {
	  L[NumL] = k;
	  NumL++;
	}
	else {
	  U[NumU] = k;
	  NumU++;
	}
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumMyBlockDiagonals_++;
    if (NumL) L_Graph_->InsertMyIndices(i, NumL, L);
    if (NumU) U_Graph_->InsertMyIndices(i, NumU, U);
    
  }

  delete [] L;
  delete [] U;

  if (LevelFill_ > 0) {

    // Complete Fill steps
    Epetra_BlockMap * L_DomainMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
    Epetra_BlockMap * L_RangeMap = (Epetra_BlockMap *) &Graph_.RangeMap();
    Epetra_BlockMap * U_DomainMap = (Epetra_BlockMap *) &Graph_.DomainMap();
    Epetra_BlockMap * U_RangeMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
    EPETRA_CHK_ERR(L_Graph_->TransformToLocal(L_DomainMap, L_RangeMap));
    EPETRA_CHK_ERR(U_Graph_->TransformToLocal(U_DomainMap, U_RangeMap));

    // At this point L_Graph and U_Graph are filled with the pattern of input graph, 
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...

    int MaxRC = NumMyBlockRows_;
    int *LinkList = new int[MaxRC];
    int *CurrentLevel = new int[MaxRC];
    int **Levels = new int*[MaxRC];
    int *CurrentRow = new int[MaxRC];
    int *LevelsRowU = new int[MaxRC];

    for (i=0; i<NumMyBlockRows_; i++) Levels[i] = 0; // Initialize Levels

    for (i=0; i<NumMyBlockRows_; i++)
    {
      int First, Next, j;
      
      // copy column indices of row into workspace and sort them
      
      int LenL = L_Graph_->NumMyIndices(i);
      int LenU = U_Graph_->NumMyIndices(i);
      int Len = LenL + LenU + 1;
      
      EPETRA_CHK_ERR(L_Graph_->ExtractMyRowCopy(i, LenL, LenL, CurrentRow));      // Get L Indices
      CurrentRow[LenL] = i;                                     // Put in Diagonal
      //EPETRA_CHK_ERR(U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1)); // Get U Indices
      int ierr1 = U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1); // Get U Indices
      if (ierr1!=0) {
	cout << "ierr1 = "<< ierr1 << endl;
	cout << "i = " << i << endl;
	cout << "NumMyBlockRows_ = " << U_Graph_->NumMyBlockRows() << endl;
      }
      
      // Construct linked list for current row
      
      for (j=0; j<Len-1; j++) {
	LinkList[CurrentRow[j]] = CurrentRow[j+1];
	CurrentLevel[CurrentRow[j]] = 0;
      }
      
      LinkList[CurrentRow[Len-1]] = NumMyBlockRows_;
      CurrentLevel[CurrentRow[Len-1]] = 0;
      
      // Merge List with rows in U
      
      First = CurrentRow[0];
      Next = First;
      while (Next < i)
        {
	  int PrevInList = Next;
	  int NextInList = LinkList[Next];
	  int RowU = Next;
	  int LengthRowU;
	  int * IndicesU;
	  // Get Indices for this row of U
	  EPETRA_CHK_ERR(U_Graph_->ExtractMyRowView(RowU, LengthRowU, IndicesU));

	  int ii;
	  
	  // Scan RowU
	  
	  for (ii=0; ii<LengthRowU; /*nop*/)
            {
	      int CurInList = IndicesU[ii];
	      if (CurInList < NextInList)
                {
		  // new fill-in
		  int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
		  if (NewLevel <= LevelFill_)
                    {
		      LinkList[PrevInList]  = CurInList;
		      LinkList[CurInList] = NextInList;
		      PrevInList = CurInList;
		      CurrentLevel[CurInList] = NewLevel;
                    }
		  ii++;
                }
	      else if (CurInList == NextInList)
                {
		  PrevInList = NextInList;
		  NextInList = LinkList[PrevInList];
		  int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
		  CurrentLevel[CurInList] = EPETRA_MIN(CurrentLevel[CurInList], NewLevel);
		  ii++;
                }
	      else // (CurInList > NextInList)
                {
		  PrevInList = NextInList;
		  NextInList = LinkList[PrevInList];
                }
            }
	  Next = LinkList[Next];
        }
      
      // Put pattern into L and U
      
      LenL = 0;

      Next = First;
      
      // Lower

      while (Next < i) {	  
	CurrentRow[LenL++] = Next;
	Next = LinkList[Next];
      }

      EPETRA_CHK_ERR(L_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
      int ierr11 = L_Graph_->InsertMyIndices(i, LenL, CurrentRow);
      if (ierr11 < 0) EPETRA_CHK_ERR(ierr1);

      // Diagonal

      if (Next != i) return(-2); // Fatal:  U has zero diagonal.
      else {
	LevelsRowU[0] = CurrentLevel[Next];
	Next = LinkList[Next];
      }
       

      // Upper

      LenU = 0;

      while (Next < NumMyBlockRows_) // Should be "Next < NumMyBlockRows_"?
        {
	  LevelsRowU[LenU+1] = CurrentLevel[Next];
	  CurrentRow[LenU++] = Next;
	  Next = LinkList[Next];
        }

      EPETRA_CHK_ERR(U_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
      int ierr2 = U_Graph_->InsertMyIndices(i, LenU, CurrentRow);
      if (ierr2<0) EPETRA_CHK_ERR(ierr2);

      // Allocate and fill Level info for this row
      Levels[i] = new int[LenU+1];
      for (int jj=0; jj<LenU+1; jj++) Levels[i][jj] = LevelsRowU[jj];

    }
    
    delete [] LinkList;
    delete [] CurrentLevel;

    for (i=0; i<NumMyBlockRows_; i++) if (Levels[i]!=0) delete [] Levels[i];
    delete [] Levels;
    delete [] CurrentRow;
    delete [] LevelsRowU;
    
  }

  // Complete Fill steps
  Epetra_BlockMap * L_DomainMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
  Epetra_BlockMap * L_RangeMap = (Epetra_BlockMap *) &Graph_.RangeMap();
  Epetra_BlockMap * U_DomainMap = (Epetra_BlockMap *) &Graph_.DomainMap();
  Epetra_BlockMap * U_RangeMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
  EPETRA_CHK_ERR(L_Graph_->TransformToLocal(L_DomainMap, L_RangeMap));
  EPETRA_CHK_ERR(U_Graph_->TransformToLocal(U_DomainMap, U_RangeMap));
      
  // Optimize graph storage
  
  EPETRA_CHK_ERR(L_Graph_->OptimizeStorage());
  EPETRA_CHK_ERR(U_Graph_->OptimizeStorage());

  // Compute global quantities

  NumGlobalBlockDiagonals_ = 0;

  EPETRA_CHK_ERR(L_Graph_->Comm().SumAll(&NumMyBlockDiagonals_, &NumGlobalBlockDiagonals_, 1));

  NumGlobalNonzeros_ = L_Graph_->NumGlobalNonzeros()+U_Graph_->NumGlobalNonzeros();
  NumMyNonzeros_ = L_Graph_->NumMyNonzeros()+U_Graph_->NumMyNonzeros();
  NumGlobalEntries_ = L_Graph_->NumGlobalEntries()+U_Graph_->NumGlobalEntries();
  NumMyEntries_ = L_Graph_->NumMyEntries()+U_Graph_->NumMyEntries();
  return(ierr);
}
//==========================================================================

// Non-member functions

ostream& operator << (ostream& os, const Ifpack_IlukGraph& A)
{
/*  Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
  Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12); */
  int LevelFill = A.LevelFill();
  Epetra_CrsGraph & L = (Epetra_CrsGraph &) A.L_Graph();
  Epetra_CrsGraph & U = (Epetra_CrsGraph &) A.U_Graph();
  os.width(14);
  os <<  "     Level of Fill = "; os << LevelFill;
  os << endl;

  os.width(14);
  os <<  "     Graph of L = "; 
  os << endl;
  os << L; // Let Epetra_CrsGraph handle the rest.

  os.width(14);
  os <<  "     Graph of U = "; 
  os << endl;
  os << U; // Let Epetra_CrsGraph handle the rest.
 
  // Reset os flags

/*  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp); */

  return os;
}
