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

#ifndef TIFPACK_ILUK_GRAPH_HPP
#define TIFPACK_ILUK_GRAPH_HPP

#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Import.hpp>

namespace Tifpack {


//! IlukGraph: A class for constructing level filled graphs for use with ILU(k) class preconditioners.

/*! The IlukGraph class enables the construction of matrix graphs using level-fill algorithms.
  The only function required for construction is a getGlobalRowView capability, i.e., the graph that is passed in
 to the constructor must implement the RowGraph interface defined in Tpetra::RowGraph.hpp 
 
 
 <b>Constructing IlukGraph objects</b>
 
 Constructing IlukGraph objects is usually a two step process of passing in a CrsGraph object and 
 an integer indicating the desired level of fill and then calling the ConstructFilledGraph function to complete the
 process.  This allows warning error codes to be returned to the calling routine.
 
 It is worth noting that an IlukGraph object creates two Tpetra::CrsGraph objects containing L and U, 
 the graphs for the lower and upper triangular parts of the ILU(k) graph.
 Thus, it is possible to manually insert and delete graph entries in L and U via the Tpetra::CrsGraph
 InsertIndices and RemoveIndices functions.  However, in this case FillComplete must be
 called before the graph is used for subsequent operations.
 
 */    
template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class IlukGraph : public Teuchos::Describable
{

public:
  
  //! IlukGraph constuctor.
  /*! Creates a IlukGraph object using the input graph and specified level of fill.  
   
   \param In
   Graph_in - An existing CrsGraph.  This object must implement the CrsGraph functions
   that provide graph dimension and pattern information.
   \param In
   LevelFill_in - The level of fill to compute via ILU(k) algorithm.
   \param In
   LevelOverlap_in - The level of overlap between subdomains.
   
   \warning Actual construction occurs in ConstructFilledGraph.  This allows error codes to 
   be passed back to the user.
   */
  IlukGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &Graph_in, int LevelFill_in, int LevelOverlap_in);
  
  //! Copy constructor.
  IlukGraph(const IlukGraph<LocalOrdinal,GlobalOrdinal,Node> & Graph_in);
  
  //! IlukGraph Destructor
  virtual ~IlukGraph();
  
  //!Set parameters using Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
   This method recogizes two parameter names: Level_fill and Level_overlap.
   Both are case insensitive, and in both cases the ParameterEntry must
   have type int.
   */
  void SetParameters(const Teuchos::ParameterList& parameterlist,
                     bool cerr_warning_if_unused=false);
  
  //! This method performs the actual construction of the graph.
  /* 
   \return Integer error code, set to 0 if successful.
   
   */
  int ConstructFilledGraph();
  
  //! Does the actual construction of the overlap matrix graph.
  /* 
   \return Integer error code, set to 0 if successful.
   
   */
  int ConstructOverlapGraph();
  
  //! Returns the level of fill used to construct this graph.
  int getLevelFill() const {return(LevelFill_);}
  
  //! Returns the level of overlap used to construct this graph.
  int getLevelOverlap() const {return(LevelOverlap_);}
  
  //! Returns the number of global matrix rows.
  GlobalOrdinal getNumGlobalRows() const {return(NumGlobalRows_);}
  
  //! Returns the number of global matrix columns.
  GlobalOrdinal getNumGlobalCols() const {return(NumGlobalCols_);}
  
  //! Returns the number of nonzero entries in the global graph.
  GlobalOrdinal getNumGlobalNonzeros() const {return(NumGlobalNonzeros_);}
  
  //! Returns the number of diagonal entries found in the local input graph.
  int NumGlobalDiagonals() const {return(NumGlobalDiagonals_);}
  
  //! Returns the number of local matrix rows.
  LocalOrdinal getNumMyRows() const {return(NumMyRows_);}
  
  //! Returns the number of local matrix columns.
  LocalOrdinal getNumMyCols() const {return(NumMyCols_);}
  
  //! Returns the number of nonzero entries in the local graph.
  size_t getNumMyNonzeros() const {return(NumMyNonzeros_);}
  
  //! Returns the number of diagonal entries found in the local input graph.
  LocalOrdinal getNumMyDiagonals() const {return(NumMyDiagonals_);}
  
  //! Returns the index base for row and column indices for this graph.
  int getIndexBase() const {return(IndexBase_);}
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getL_Graph() {return(*L_Graph_);}
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
   Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & U_Graph() {return(*U_Graph_);}
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
   Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & L_Graph() const {return(*L_Graph_);}
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
   Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & U_Graph() const {return(*U_Graph_);}
  
  //! Returns the importer used to create the overlapped graph.
   Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> * OverlapImporter() const  {return(&*OverlapImporter_);}
  
  //! Returns the the overlapped graph.
   Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> * OverlapGraph() const  {return(&*OverlapGraph_);}
  
  //! Returns the Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> object associated with the domain of this matrix operator.
   const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & DomainMap() const {return(DomainMap_);}
  
  //! Returns the Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> object associated with the range of this matrix operator.
   const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & RangeMap() const{return(RangeMap_);}
  
private:
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> TpetraMapType;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> TpetraCrsGraphType;
  
  Teuchos::RCP<const TpetraCrsGraphType > Graph_;
  Teuchos::RCP<TpetraMapType> DomainMap_;
  Teuchos::RCP<TpetraMapType> RangeMap_;
  Teuchos::RCP<Teuchos::Comm<int> > Comm_;
  Teuchos::RCP<TpetraCrsGraphType > OverlapGraph_;
  Teuchos::RCP<TpetraMapType> OverlapRowMap_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > OverlapImporter_;
  int LevelFill_;
  int LevelOverlap_;
  Teuchos::RCP<TpetraCrsGraphType > L_Graph_;
  Teuchos::RCP<TpetraCrsGraphType > U_Graph_;
  int IndexBase_;
  int NumGlobalRows_;
  int NumGlobalCols_;
  int NumGlobalBlockRows_;
  int NumGlobalBlockCols_;
  int NumGlobalDiagonals_;
  int NumGlobalBlockDiagonals_;
  int NumGlobalNonzeros_;
  int NumGlobalEntries_;
  int NumMyBlockRows_;
  int NumMyBlockCols_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  int NumMyBlockDiagonals_;
  int NumMyNonzeros_;
  int NumMyEntries_;
  
  
};

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::IlukGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph_in, int LevelFill_in, int LevelOverlap_in)
: Graph_(Graph_in),
  DomainMap_(Graph_in.DomainMap()),
  RangeMap_(Graph_in.RangeMap()),
  Comm_(Graph_in.Comm()),
  LevelFill_(LevelFill_in),
  LevelOverlap_(LevelOverlap_in),
  IndexBase_(Graph_in.IndexBase()),
  NumGlobalRows_(Graph_in.NumGlobalRows()),
  NumGlobalCols_(Graph_in.NumGlobalCols()),
  NumGlobalBlockRows_(Graph_in.NumGlobalBlockRows()),
  NumGlobalBlockCols_(Graph_in.NumGlobalBlockCols()),
  NumGlobalDiagonals_(0),
  NumGlobalBlockDiagonals_(0),
  NumGlobalNonzeros_(0),
  NumGlobalEntries_(0),
  NumMyBlockRows_(Graph_in.NumMyBlockRows()),
  NumMyBlockCols_(Graph_in.NumMyBlockCols()),
  NumMyRows_(Graph_in.NumMyRows()),
  NumMyCols_(Graph_in.NumMyCols()),
  NumMyDiagonals_(0),
  NumMyBlockDiagonals_(0),
  NumMyNonzeros_(0),
  NumMyEntries_(0)
{
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::IlukGraph(const IlukGraph<LocalOrdinal,GlobalOrdinal,Node> & Graph_in) 
: Graph_(Graph_in.Graph_),
  DomainMap_(Graph_in.DomainMap()),
  RangeMap_(Graph_in.RangeMap()),
  Comm_(Graph_in.Comm()),
  OverlapGraph_(Graph_in.OverlapGraph_),
  OverlapRowMap_(Graph_in.OverlapRowMap_),
  OverlapImporter_(Graph_in.OverlapImporter_),
  LevelFill_(Graph_in.LevelFill_),
  LevelOverlap_(Graph_in.LevelOverlap_),
  IndexBase_(Graph_in.IndexBase_),
  NumGlobalRows_(Graph_in.NumGlobalRows_),
  NumGlobalCols_(Graph_in.NumGlobalCols_),
  NumGlobalBlockRows_(Graph_in.NumGlobalBlockRows_),
  NumGlobalBlockCols_(Graph_in.NumGlobalBlockCols_),
  NumGlobalDiagonals_(Graph_in.NumGlobalDiagonals_),
  NumGlobalBlockDiagonals_(Graph_in.NumGlobalBlockDiagonals_),
  NumGlobalNonzeros_(Graph_in.NumGlobalNonzeros_),
  NumGlobalEntries_(Graph_in.NumGlobalEntries_),
  NumMyBlockRows_(Graph_in.NumMyBlockRows_),
  NumMyBlockCols_(Graph_in.NumMyBlockCols_),
  NumMyRows_(Graph_in.NumMyRows_),
  NumMyCols_(Graph_in.NumMyCols_),
  NumMyDiagonals_(Graph_in.NumMyDiagonals_),
  NumMyBlockDiagonals_(Graph_in.NumMyBlockDiagonals_),
  NumMyNonzeros_(Graph_in.NumMyNonzeros_),
  NumMyEntries_(Graph_in.NumMyEntries_)
{
  TpetraCrsGraphType & L_Graph_In = Graph_in.L_Graph();
  TpetraCrsGraphType & U_Graph_In = Graph_in.U_Graph();
  L_Graph_ = Teuchos::rcp( new TpetraCrsGraphType(L_Graph_In) );
  U_Graph_ = Teuchos::rcp( new TpetraCrsGraphType(U_Graph_In) );
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::~IlukGraph()
{
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::SetParameters(const Teuchos::ParameterList& parameterlist,
                                                               bool cerr_warning_if_unused)
{
  std::string paramname("LEVEL_FILL");
  if (parameterlist.isParameter(paramname)) {
    if (parameterlist.isType<int>(paramname)) {
      LevelFill_ = parameterlist.get<int>(paramname);
    }
  }

  paramname = "LEVEL_OVERLAP";
  if (parameterlist.isParameter(paramname)) {
    if (parameterlist.isType<int>(paramname)) {
      LevelOverlap_ = parameterlist.get<int>(paramname);
    }
  }
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
int IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::ConstructOverlapGraph() {
  
  OverlapGraph_ = Teuchos::rcp( (TpetraCrsGraphType *) &Graph_, false );
  OverlapRowMap_ = Teuchos::rcp( (TpetraMapType *) &Graph_.RowMap(), false );
  
  if (LevelOverlap_==0 || !Graph_.DomainMap().DistributedGlobal()) return(0); // Nothing to do
  
  Teuchos::RCP<TpetraCrsGraphType> OldGraph;
  Teuchos::RCP<TpetraMapType > OldRowMap;
  TpetraMapType * DomainMap_tmp = (TpetraMapType *) &Graph_.DomainMap();
  TpetraMapType * RangeMap_tmp = (TpetraMapType *) &Graph_.RangeMap();
  for (int level=1; level <= LevelOverlap_; level++) {
  	OldGraph = OverlapGraph_; 
  	OldRowMap = OverlapRowMap_;
  	
  	OverlapImporter_ = OldGraph->getImporter();
  	OverlapRowMap_ = OverlapImporter_->getTargetMap();
  	
  	
  	if (level<LevelOverlap_)
  		OverlapGraph_ = Teuchos::rcp( new TpetraCrsGraphType(OverlapRowMap_, 0) );
  	else
  		// On last iteration, we want to filter out all columns except those that correspond
  		// to rows in the graph.  This assures that our matrix is square
  		OverlapGraph_ = Teuchos::rcp( new TpetraCrsGraphType(OverlapRowMap_, *OverlapRowMap_, 0) );
 
  	EPETRA_CHK_ERR(OverlapGraph_->Import( Graph_, *OverlapImporter_, Insert));
  	if (level<LevelOverlap_) {
  		EPETRA_CHK_ERR(OverlapGraph_->FillComplete(*DomainMap_tmp, *RangeMap_tmp));
  	}
  	else {
  		// Copy last OverlapImporter because we will use it later
  		OverlapImporter_ = Teuchos::rcp( new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(*OverlapRowMap_, *DomainMap_tmp) );
  		EPETRA_CHK_ERR(OverlapGraph_->FillComplete(*DomainMap_tmp, *RangeMap_tmp));
  	}
  }

  NumMyBlockRows_ = OverlapGraph_->NumMyBlockRows();
  NumMyBlockCols_ = OverlapGraph_->NumMyBlockCols();
  NumMyRows_ = OverlapGraph_->NumMyRows();
  NumMyCols_ = OverlapGraph_->NumMyCols();

  return(0);
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
int IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::ConstructFilledGraph() {
  int ierr = 0;
  int i, j;
  int * In=0;
  int NumIn, NumL, NumU;
  bool DiagFound;
  
  
  EPETRA_CHK_ERR(ConstructOverlapGraph());
  
  L_Graph_ = Teuchos::rcp( new Tpetra::CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0) );
  U_Graph_ = Teuchos::rcp( new Tpetra::CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0));
  
  
  // Get Maximun Row length
  int MaxNumIndices = OverlapGraph_->MaxNumIndices();
  
  vector<int> L(MaxNumIndices);
  vector<int> U(MaxNumIndices);
  
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
  	if (NumL) L_Graph_->InsertMyIndices(i, NumL, &L[0]);
  	if (NumU) U_Graph_->InsertMyIndices(i, NumU, &U[0]);
  	
  }
  
  if (LevelFill_ > 0) {
  	
  	// Complete Fill steps
  	TpetraMapType * L_DomainMap = (TpetraMapType *) &OverlapGraph_->RowMap();
  	TpetraMapType * L_RangeMap = (TpetraMapType *) &Graph_.RangeMap();
  	TpetraMapType * U_DomainMap = (TpetraMapType *) &Graph_.DomainMap();
  	TpetraMapType * U_RangeMap = (TpetraMapType *) &OverlapGraph_->RowMap();
  	EPETRA_CHK_ERR(L_Graph_->FillComplete(*L_DomainMap, *L_RangeMap));
  	EPETRA_CHK_ERR(U_Graph_->FillComplete(*U_DomainMap, *U_RangeMap));
  	
  	// At this point L_Graph and U_Graph are filled with the pattern of input graph, 
  	// sorted and have redundant indices (if any) removed.  Indices are zero based.
  	// LevelFill is greater than zero, so continue...
  	
  	int MaxRC = NumMyBlockRows_;
  	vector<vector<int> > Levels(MaxRC);
  	vector<int> LinkList(MaxRC);
  	vector<int> CurrentLevel(MaxRC);
  	vector<int> CurrentRow(MaxRC);
  	vector<int> LevelsRowU(MaxRC);
  	
  	for (i=0; i<NumMyBlockRows_; i++)
  	{
  		int First, Next;
  		
  		// copy column indices of row into workspace and sort them
  		
  		int LenL = L_Graph_->NumMyIndices(i);
  		int LenU = U_Graph_->NumMyIndices(i);
  		int Len = LenL + LenU + 1;
  		
  		EPETRA_CHK_ERR(L_Graph_->ExtractMyRowCopy(i, LenL, LenL, &CurrentRow[0]));      // Get L Indices
  		CurrentRow[LenL] = i;                                     // Put in Diagonal
  		//EPETRA_CHK_ERR(U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1)); // Get U Indices
  		int ierr1 = 0;
  		if (LenU) {
  			// Get U Indices
  			ierr1 = U_Graph_->ExtractMyRowCopy(i, LenU, LenU, &CurrentRow[LenL+1]);
  		}
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
  		int ierr11 = L_Graph_->InsertMyIndices(i, LenL, &CurrentRow[0]);
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
  		int ierr2 = U_Graph_->InsertMyIndices(i, LenU, &CurrentRow[0]);
  		if (ierr2<0) EPETRA_CHK_ERR(ierr2);
  		
  		// Allocate and fill Level info for this row
  		Levels[i] = vector<int>(LenU+1);
  		for (int jj=0; jj<LenU+1; jj++) Levels[i][jj] = LevelsRowU[jj];
  		
  	}
  }    
  
  // Complete Fill steps
  Teuchos::RCP<TpetraMapType> L_DomainMap = OverlapGraph_->RowMap();
  Teuchos::RCP<TpetraMapType> L_RangeMap  = Graph_.RangeMap();
  Teuchos::RCP<TpetraMapType> U_DomainMap = Graph_.DomainMap();
  Teuchos::RCP<TpetraMapType> U_RangeMap  = OverlapGraph_->RowMap();
  L_Graph_->fillComplete(L_DomainMap, L_RangeMap);
  U_Graph_->fillComplete(U_DomainMap, U_RangeMap);
  
  // Optimize graph storage
  
  EPETRA_CHK_ERR(L_Graph_->OptimizeStorage());
  EPETRA_CHK_ERR(U_Graph_->OptimizeStorage());
  
  // Compute global quantities
  
  NumGlobalBlockDiagonals_ = 0;
 
  Teuchos::reduceAll(L_Graph_->Comm(), Teuchos::REDUCE_SUM, 1, &NumMyBlockDiagonals_, &NumGlobalBlockDiagonals_);
 
  NumGlobalNonzeros_ = L_Graph_->NumGlobalNonzeros()+U_Graph_->NumGlobalNonzeros();
  NumMyNonzeros_ = L_Graph_->NumMyNonzeros()+U_Graph_->NumMyNonzeros();
  NumGlobalEntries_ = L_Graph_->NumGlobalEntries()+U_Graph_->NumGlobalEntries();
  NumMyEntries_ = L_Graph_->NumMyEntries()+U_Graph_->NumMyEntries();

  return(ierr);
}

} // namespace Tifpack

#endif /* TIFPACK_ILUK_GRAPH_HPP */
