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

#include <vector>
#include <algorithm>

#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Import.hpp>
#include <Tifpack_CreateOverlapGraph.hpp>
#include <Tifpack_Parameters.hpp>

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
  /** 
   */
  void ConstructFilledGraph();
  
  //! Does the actual construction of the overlap matrix graph.
  /* 
   \return Integer error code, set to 0 if successful.
   
   */
  void ConstructOverlapGraph();
  
  //! Returns the level of fill used to construct this graph.
  int getLevelFill() const {return(LevelFill_);}
  
  //! Returns the level of overlap used to construct this graph.
  int getLevelOverlap() const {return(LevelOverlap_);}
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getL_Graph() {return(*L_Graph_);}
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getU_Graph() {return(*U_Graph_);}
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getL_Graph() const {return(*L_Graph_);}
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getU_Graph() const {return(*U_Graph_);}
  
  //! Returns the the overlapped graph.
  Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> * getOverlapGraph() const  {return(&*OverlapGraph_);}

  //! Returns the global number of diagonals in the ILU(k) graph.
  size_t getNumGlobalDiagonals() const { return NumGlobalDiagonals_; }

private:
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> TpetraMapType;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> TpetraCrsGraphType;
  
  Teuchos::RCP<const TpetraCrsGraphType > Graph_;
  Teuchos::RCP<const TpetraCrsGraphType > OverlapGraph_;
  int LevelFill_;
  int LevelOverlap_;
  Teuchos::RCP<TpetraCrsGraphType > L_Graph_;
  Teuchos::RCP<TpetraCrsGraphType > U_Graph_;
  size_t NumMyDiagonals_;
  size_t NumGlobalDiagonals_;
};

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::IlukGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph_in, int LevelFill_in, int LevelOverlap_in)
: Graph_(Graph_in),
  OverlapGraph_(),
  LevelFill_(LevelFill_in),
  LevelOverlap_(LevelOverlap_in),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0)
{
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::IlukGraph(const IlukGraph<LocalOrdinal,GlobalOrdinal,Node> & Graph_in) 
: Graph_(Graph_in.Graph_),
  OverlapGraph_(Graph_in.OverlapGraph_),
  LevelFill_(Graph_in.LevelFill_),
  LevelOverlap_(Graph_in.LevelOverlap_),
  L_Graph_(),
  U_Graph_(),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0)
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
  GetParameter(parameterlist, "LEVEL_FILL", LevelFill_);
  GetParameter(parameterlist, "LEVEL_OVERLAP", LevelOverlap_);
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::ConstructOverlapGraph() {
  
  if (OverlapGraph_ == Teuchos::null) {
    OverlapGraph_ = CreateOverlapGraph(Graph_, LevelOverlap_);
  }
}

//==============================================================================
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void IlukGraph<LocalOrdinal,GlobalOrdinal,Node>::ConstructFilledGraph() {
  size_t NumIn, NumL, NumU;
  bool DiagFound;
 
  ConstructOverlapGraph();
 
  L_Graph_ = Teuchos::rcp( new TpetraCrsGraphType(OverlapGraph_->getRowMap(), OverlapGraph_->getRowMap(),  0));
  U_Graph_ = Teuchos::rcp( new TpetraCrsGraphType(OverlapGraph_->getRowMap(), OverlapGraph_->getRowMap(),  0));
 
 
  // Get Maximum Row length
  int MaxNumIndices = OverlapGraph_->getNodeMaxNumRowEntries();
 
  Teuchos::Array<int> L(MaxNumIndices);
  Teuchos::Array<int> U(MaxNumIndices);
 
  // First we copy the user's graph into L and U, regardless of fill level
 
  int NumMyRows = OverlapGraph_->getRowMap()->getNodeNumElements();
  NumMyDiagonals_ = 0;

  for (int i=0; i< NumMyRows; i++) {
 
    Teuchos::ArrayRCP<const LocalOrdinal> my_indices = OverlapGraph_->getLocalRowView(i);
 
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    NumIn = my_indices.size();

    for (size_t j=0; j< NumIn; j++) {
      LocalOrdinal k = my_indices[j];
      
      if (k<NumMyRows) { // Ignore column elements that are not in the square matrix
        
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
 
    if (DiagFound) ++NumMyDiagonals_;
    if (NumL) {
      Teuchos::ArrayView<LocalOrdinal> Lview(&L[0], NumL);
      L_Graph_->insertLocalIndices(i, Lview );
    }
    if (NumU) {
      Teuchos::ArrayView<LocalOrdinal> Uview(&U[0], NumU);
      U_Graph_->insertLocalIndices(i, Uview );
    }
  }
 
  if (LevelFill_ > 0) {
    
    // Complete Fill steps
    Teuchos::RCP<const TpetraMapType> L_DomainMap = OverlapGraph_->getRowMap();
    Teuchos::RCP<const TpetraMapType> L_RangeMap = Graph_->getRangeMap();
    Teuchos::RCP<const TpetraMapType> U_DomainMap = Graph_->getDomainMap();
    Teuchos::RCP<const TpetraMapType> U_RangeMap = OverlapGraph_->getRowMap();
    L_Graph_->fillComplete(L_DomainMap, L_RangeMap, Tpetra::DoNotOptimizeStorage);
    U_Graph_->fillComplete(U_DomainMap, U_RangeMap, Tpetra::DoNotOptimizeStorage);
    
    // At this point L_Graph and U_Graph are filled with the pattern of input graph, 
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...
    
    int MaxRC = NumMyRows;
    std::vector<std::vector<int> > Levels(MaxRC);
    std::vector<int> LinkList(MaxRC);
    std::vector<int> CurrentLevel(MaxRC);
    Teuchos::Array<int> CurrentRow(MaxRC+1);
    std::vector<int> LevelsRowU(MaxRC);
 
    for (int i=0; i<NumMyRows; i++)
    {
      int First, Next;

      // copy column indices of row into workspace and sort them

      size_t LenL = L_Graph_->getNumEntriesInLocalRow(i);
      size_t LenU = U_Graph_->getNumEntriesInLocalRow(i);
      size_t Len = LenL + LenU + 1;

      CurrentRow.resize(Len);

      L_Graph_->getLocalRowCopy(i, CurrentRow(), LenL);      // Get L Indices
      CurrentRow[LenL] = i;                                     // Put in Diagonal
      if (LenU > 0) {
        Teuchos::ArrayView<LocalOrdinal> URowView(&CurrentRow[LenL+1], LenU);
        // Get U Indices
        U_Graph_->getLocalRowCopy(i, URowView, LenU);
      }

      // Construct linked list for current row
      
      for (size_t j=0; j<Len-1; j++) {
        LinkList[CurrentRow[j]] = CurrentRow[j+1];
        CurrentLevel[CurrentRow[j]] = 0;
      }
      
      LinkList[CurrentRow[Len-1]] = NumMyRows;
      CurrentLevel[CurrentRow[Len-1]] = 0;
      
      // Merge List with rows in U
      
      First = CurrentRow[0];
      Next = First;
      while (Next < i)
      {
        int PrevInList = Next;
        int NextInList = LinkList[Next];
        int RowU = Next;
        // Get Indices for this row of U
        Teuchos::ArrayRCP<const LocalOrdinal> IndicesU = U_Graph_->getLocalRowView(RowU);
        int LengthRowU = IndicesU.size();
        
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
            CurrentLevel[CurInList] = std::min(CurrentLevel[CurInList], NewLevel);
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
      
      CurrentRow.resize(0);
      
      Next = First;
      
      // Lower
      
      while (Next < i) {    
        CurrentRow.push_back(Next);
        Next = LinkList[Next];
      }

      L_Graph_->removeLocalIndices(i); // Delete current set of Indices
      if (CurrentRow.size() > 0) {
        L_Graph_->insertLocalIndices(i, CurrentRow());
      }
 
      // Diagonal
      
      TEST_FOR_EXCEPTION(Next != i, std::runtime_error,
                         "Tifpack::IlukGraph::ConstructFilledGraph: FATAL: U has zero diagonal")

      LevelsRowU[0] = CurrentLevel[Next];
      Next = LinkList[Next];
      
      // Upper
      
      CurrentRow.resize(0);
      LenU = 0;
      
      while (Next < NumMyRows) {
        LevelsRowU[LenU+1] = CurrentLevel[Next];
        CurrentRow.push_back(Next);
        ++LenU;
        Next = LinkList[Next];
      }
      
      U_Graph_->removeLocalIndices(i); // Delete current set of Indices
      if (LenU > 0) {
        U_Graph_->insertLocalIndices(i, CurrentRow());
      }
 
      // Allocate and fill Level info for this row
      Levels[i] = std::vector<int>(LenU+1);
      for (size_t jj=0; jj<LenU+1; jj++) {
        Levels[i][jj] = LevelsRowU[jj];
      }
    }
  }    
  
  // Complete Fill steps
  Teuchos::RCP<const TpetraMapType> L_DomainMap = OverlapGraph_->getRowMap();
  Teuchos::RCP<const TpetraMapType> L_RangeMap  = Graph_->getRangeMap();
  Teuchos::RCP<const TpetraMapType> U_DomainMap = Graph_->getDomainMap();
  Teuchos::RCP<const TpetraMapType> U_RangeMap  = OverlapGraph_->getRowMap();
  L_Graph_->fillComplete(L_DomainMap, L_RangeMap);//DoOptimizeStorage is default here...
  U_Graph_->fillComplete(U_DomainMap, U_RangeMap);//DoOptimizeStorage is default here...

  Teuchos::reduceAll<int,size_t>(*(L_DomainMap->getComm()), Teuchos::REDUCE_SUM, 1, &NumMyDiagonals_, &NumGlobalDiagonals_);
}

} // namespace Tifpack

#endif /* TIFPACK_ILUK_GRAPH_HPP */
