#include "Ifpack_IlukGraph.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"
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
    NumGlobalDiagonals_(0),
    NumGlobalNonzeros_(0),
    NumMyRows_(Graph.NumMyRows()),
    NumMyCols_(Graph.NumMyCols()),
    NumMyDiagonals_(0),
    NumMyNonzeros_(0)
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
    NumGlobalDiagonals_(Graph.NumGlobalDiagonals_),
    NumGlobalNonzeros_(Graph.NumGlobalNonzeros_),
    NumMyRows_(Graph.NumMyRows_),
    NumMyCols_(Graph.NumMyCols_),
    NumMyDiagonals_(Graph.NumMyDiagonals_),
    NumMyNonzeros_(Graph.NumMyNonzeros_)
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

//==============================================================================
int Ifpack_IlukGraph::ConstructOverlapGraph() {

  OverlapGraph_ = (Epetra_CrsGraph *) &Graph_;
  OverlapRowMap_ = (Epetra_BlockMap *) &Graph_.RowMap();
  
  if (LevelOverlap_==0 || !Graph_.DomainMap().DistributedGlobal()) return(0); // Nothing to do

  Epetra_CrsGraph * OldGraph;
  Epetra_BlockMap * OldRowMap;
  Epetra_BlockMap * DomainMap = (Epetra_BlockMap *) &Graph_.DomainMap();
  for (int level=1; level <= LevelOverlap_; level++) {
    OldGraph = OverlapGraph_; 
    OldRowMap = OverlapRowMap_;

    OverlapImporter_ = (Epetra_Import *) OldGraph->Importer();
    OverlapRowMap_ = new Epetra_BlockMap(OverlapImporter_->TargetMap());

    OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0);
    assert(OverlapGraph_->Import( Graph_, *OverlapImporter_, Insert)==0);
    if (level<LevelOverlap_) 
      assert(OverlapGraph_->TransformToLocal(DomainMap, OverlapRowMap_)==0);
    else {
      // Copy last OverlapImporter because we will use it later
      OverlapImporter_ = new Epetra_Import(*OverlapRowMap_, *DomainMap);
      assert(OverlapGraph_->TransformToLocal()==0);
    }

    if (OldGraph!=&Graph_) delete OldGraph;
    if (OldRowMap!=&Graph_.RowMap()) delete OldRowMap;
  }

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

  
  assert(ConstructOverlapGraph()==0);

  L_Graph_ = new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), 0);
  U_Graph_ = new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), 0);


  // Get Maximun Row length
  int MaxNumIndices = OverlapGraph_->MaxNumIndices();

  L = new int[MaxNumIndices];
  U = new int[MaxNumIndices];
    

  // First we copy the user's graph into L and U, regardless of fill level

  for (i=0; i< NumMyRows_; i++) {


    OverlapGraph_->ExtractMyRowView(i, NumIn, In); // Get Indices

    
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = In[j];

      if (k<NumMyRows_) { // Ignore column elements that are not in the square matrix

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

    if (DiagFound) NumMyDiagonals_++;
    if (NumL) L_Graph_->InsertMyIndices(i, NumL, L);
    if (NumU) U_Graph_->InsertMyIndices(i, NumU, U);
    
  }

  delete [] L;
  delete [] U;

  if (LevelFill_ > 0) {

    // Complete Fill steps
    L_Graph_->TransformToLocal();
    U_Graph_->TransformToLocal();

    // At this point L_Graph and U_Graph are filled with the pattern of input graph, 
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...

    int MaxRC = NumMyRows_;
    int *LinkList = new int[MaxRC];
    int *CurrentLevel = new int[MaxRC];
    int **Levels = new int*[MaxRC];
    int *CurrentRow = new int[MaxRC];
    int *LevelsRowU = new int[MaxRC];

    for (i=0; i<NumMyRows_; i++) Levels[i] = 0; // Initialize Levels

    for (i=0; i<NumMyRows_; i++)
    {
      int First, Next, j;
      
      // copy column indices of row into workspace and sort them
      
      int LenL = L_Graph_->NumMyIndices(i);
      int LenU = U_Graph_->NumMyIndices(i);
      int Len = LenL + LenU + 1;
      
      assert(L_Graph_->ExtractMyRowCopy(i, LenL, LenL, CurrentRow)==0);      // Get L Indices
      CurrentRow[LenL] = i;                                     // Put in Diagonal
      //assert(U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1)==0); // Get U Indices
      int ierr1 = U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1); // Get U Indices
      if (ierr1!=0) {
	cout << "ierr1 = "<< ierr1 << endl;
	cout << "i = " << i << endl;
	cout << "NumMyRows_ = " << U_Graph_->NumMyRows() << endl;
      }
      
      // Construct linked list for current row
      
      for (j=0; j<Len-1; j++) {
	LinkList[CurrentRow[j]] = CurrentRow[j+1];
	CurrentLevel[CurrentRow[j]] = 0;
      }
      
      LinkList[CurrentRow[Len-1]] = NumMyRows_;
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
	  assert(U_Graph_->ExtractMyRowView(RowU, LengthRowU, IndicesU)==0);

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

      assert(L_Graph_->RemoveMyIndices(i)==0); // Delete current set of Indices
      assert(L_Graph_->InsertMyIndices(i, LenL, CurrentRow)>=0);

      // Diagonal

      if (Next != i) return(-2); // Fatal:  U has zero diagonal.
      else {
	LevelsRowU[0] = CurrentLevel[Next];
	Next = LinkList[Next];
      }
       

      // Upper

      LenU = 0;

      while (Next < NumMyRows_) // Should be "Next < NumMyRows_"?
        {
	  LevelsRowU[LenU+1] = CurrentLevel[Next];
	  CurrentRow[LenU++] = Next;
	  Next = LinkList[Next];
        }

      assert(U_Graph_->RemoveMyIndices(i)==0); // Delete current set of Indices
      assert(U_Graph_->InsertMyIndices(i, LenU, CurrentRow)>=0);

      // Allocate and fill Level info for this row
      Levels[i] = new int[LenU+1];
      for (int jj=0; jj<LenU+1; jj++) Levels[i][jj] = LevelsRowU[jj];

    }
    
    delete [] LinkList;
    delete [] CurrentLevel;

    for (i=0; i<NumMyRows_; i++) if (Levels[i]!=0) delete Levels[i];
    delete [] Levels;
    delete [] CurrentRow;
    delete [] LevelsRowU;
    
  }

  // Complete Fill steps
  assert(L_Graph_->TransformToLocal()==0);
  assert(U_Graph_->TransformToLocal()==0);
    
  // Optimize graph storage
  
  assert(L_Graph_->OptimizeStorage()==0);
  assert(U_Graph_->OptimizeStorage()==0);

  // Compute global quantities

  NumGlobalDiagonals_ = 0;

  assert(L_Graph_->Comm().SumAll(&NumMyDiagonals_, &NumGlobalDiagonals_, 1)==0);

  NumGlobalNonzeros_ = L_Graph_->NumGlobalNonzeros()+U_Graph_->NumGlobalNonzeros();
  NumMyNonzeros_ = L_Graph_->NumMyNonzeros()+U_Graph_->NumMyNonzeros();
  return(ierr);
}
#ifdef PETRA_LEVELSCHEDULING
//==============================================================================
int Ifpack_IlukGraph::ComputeLevels(int NumThreads)
{
  L_Graph_->ComputeLevels(NumThreads);
  U_Graph_->ComputeLevels(NumThreads);
  
  return(0);
}
#endif
//==========================================================================

// Non-member functions

ostream& operator << (ostream& os, const Ifpack_IlukGraph& A)
{
  Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
  Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12);
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

  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp);

  return os;
}
