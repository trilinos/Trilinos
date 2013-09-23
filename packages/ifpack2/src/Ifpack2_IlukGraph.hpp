/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_ILUK_GRAPH_HPP
#define IFPACK2_ILUK_GRAPH_HPP

#include <vector>
#include <algorithm>

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Import.hpp>
#include <Ifpack2_CreateOverlapGraph.hpp>
#include <Ifpack2_Parameters.hpp>

namespace Ifpack2 {


//! A class for constructing level filled graphs for use with ILU(k) class preconditioners.

/*! Ifpack2::IlukGraph enables the construction of matrix graphs using level-fill algorithms.
  The only function required for construction is a getGlobalRowView capability, i.e., the graph that is passed in
 to the constructor must implement the RowGraph interface defined in Tpetra_RowGraph.hpp 
 
 
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
template<class GraphType>
class IlukGraph : public Teuchos::Describable {
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
   
   \warning Actual construction occurs in constructFilledGraph.
   */
  IlukGraph(const Teuchos::RCP<const GraphType> &Graph_in, int LevelFill_in, int LevelOverlap_in);
  
  //! Copy constructor.
  IlukGraph(const IlukGraph<GraphType> & Graph_in);
  
  //! IlukGraph Destructor
  virtual ~IlukGraph();
  
  //!Set parameters using Teuchos::ParameterList object.
  /**
   This method recogizes two parameter names: Level_fill and Level_overlap.
   Both are case insensitive, and in both cases the ParameterEntry must
   have type int.
   */
  void setParameters(const Teuchos::ParameterList& parameterlist);
  
  //! This method performs the actual construction of the graph.
  /** 
   */
  void constructFilledGraph();
  
  //! Does the actual construction of the overlap matrix graph.
  /**
   */
  void constructOverlapGraph();
  
  //! Returns the level of fill used to construct this graph.
  int getLevelFill() const {return(LevelFill_);}
  
  //! Returns the level of overlap used to construct this graph.
  int getLevelOverlap() const {return(LevelOverlap_);}
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Teuchos::RCP<GraphType> 
  getL_Graph () const {
    return L_Graph_;
  }
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Teuchos::RCP<GraphType> 
  getU_Graph () const {
    return U_Graph_;
  }
  
  //! Returns the the overlapped graph.
  Teuchos::RCP<GraphType> 
  getOverlapGraph () const {
    return OverlapGraph_;
  }

  //! Returns the global number of diagonals in the ILU(k) graph.
  size_t getNumGlobalDiagonals() const { return NumGlobalDiagonals_; }

private:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::map_type map_type;
  
  Teuchos::RCP<const GraphType > Graph_;
  Teuchos::RCP<const GraphType > OverlapGraph_;
  int LevelFill_;
  int LevelOverlap_;
  Teuchos::RCP<GraphType> L_Graph_;
  Teuchos::RCP<GraphType> U_Graph_;
  size_t NumMyDiagonals_;
  size_t NumGlobalDiagonals_;
};

//==============================================================================
template<class GraphType>
IlukGraph<GraphType>::IlukGraph(const Teuchos::RCP<const GraphType>& Graph_in, int LevelFill_in, int LevelOverlap_in)
: Graph_(Graph_in),
  OverlapGraph_(),
  LevelFill_(LevelFill_in),
  LevelOverlap_(LevelOverlap_in),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0)
{
}

//==============================================================================
template<class GraphType>
IlukGraph<GraphType>::IlukGraph(const IlukGraph<GraphType> & Graph_in) 
: Graph_(Graph_in.Graph_),
  OverlapGraph_(Graph_in.OverlapGraph_),
  LevelFill_(Graph_in.LevelFill_),
  LevelOverlap_(Graph_in.LevelOverlap_),
  L_Graph_(),
  U_Graph_(),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0)
{
  GraphType & L_Graph_In = Graph_in.L_Graph();
  GraphType & U_Graph_In = Graph_in.U_Graph();
  L_Graph_ = Teuchos::rcp( new GraphType(L_Graph_In) );
  U_Graph_ = Teuchos::rcp( new GraphType(U_Graph_In) );
}

//==============================================================================
template<class GraphType>
IlukGraph<GraphType>::~IlukGraph()
{
}

//==============================================================================
template<class GraphType>
void IlukGraph<GraphType>::setParameters(const Teuchos::ParameterList& parameterlist)
{
  getParameter(parameterlist, "iluk level-of-fill", LevelFill_);
  getParameter(parameterlist, "iluk level-of-overlap", LevelOverlap_);
}

//==============================================================================
template<class GraphType>
void IlukGraph<GraphType>::constructOverlapGraph() {
  
  if (OverlapGraph_ == Teuchos::null) {
    OverlapGraph_ = CreateOverlapGraph<GraphType>(Graph_, LevelOverlap_);
  }
}

//==============================================================================
template<class GraphType>
void IlukGraph<GraphType>::constructFilledGraph() {
  size_t NumIn, NumL, NumU;
  bool DiagFound;
 
  constructOverlapGraph();
 
  L_Graph_ = Teuchos::rcp( new GraphType(OverlapGraph_->getRowMap(), OverlapGraph_->getRowMap(),  0));
  U_Graph_ = Teuchos::rcp( new GraphType(OverlapGraph_->getRowMap(), OverlapGraph_->getRowMap(),  0));
 
 
  // Get Maximum Row length
  int MaxNumIndices = OverlapGraph_->getNodeMaxNumRowEntries();
 
  Teuchos::Array<local_ordinal_type> L(MaxNumIndices);
  Teuchos::Array<local_ordinal_type> U(MaxNumIndices);
 
  // First we copy the user's graph into L and U, regardless of fill level
 
  int NumMyRows = OverlapGraph_->getRowMap()->getNodeNumElements();
  NumMyDiagonals_ = 0;

  for (int i=0; i< NumMyRows; i++) {
 
    Teuchos::ArrayView<const local_ordinal_type> my_indices;
    OverlapGraph_->getLocalRowView(i,my_indices);
 
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    NumIn = my_indices.size();

    for (size_t j=0; j< NumIn; j++) {
      local_ordinal_type k = my_indices[j];
      
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
      Teuchos::ArrayView<local_ordinal_type> Lview(&L[0], NumL);
      L_Graph_->insertLocalIndices(i, Lview );
    }
    if (NumU) {
      Teuchos::ArrayView<local_ordinal_type> Uview(&U[0], NumU);
      U_Graph_->insertLocalIndices(i, Uview );
    }
  }
 
  if (LevelFill_ > 0) {
    
    // Complete Fill steps
    Teuchos::RCP<const map_type> L_DomainMap = OverlapGraph_->getRowMap();
    Teuchos::RCP<const map_type> L_RangeMap = Graph_->getRangeMap();
    Teuchos::RCP<const map_type> U_DomainMap = Graph_->getDomainMap();
    Teuchos::RCP<const map_type> U_RangeMap = OverlapGraph_->getRowMap();
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
    params->set("Optimize Storage",false);
    L_Graph_->fillComplete(L_DomainMap, L_RangeMap, params);
    U_Graph_->fillComplete(U_DomainMap, U_RangeMap, params);
    L_Graph_->resumeFill();
    U_Graph_->resumeFill();
    
    // At this point L_Graph and U_Graph are filled with the pattern of input graph, 
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...
    
    int MaxRC = NumMyRows;
    std::vector<std::vector<int> > Levels(MaxRC);
    std::vector<int> LinkList(MaxRC);
    std::vector<int> CurrentLevel(MaxRC);
    Teuchos::Array<local_ordinal_type> CurrentRow(MaxRC+1);
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
        Teuchos::ArrayView<local_ordinal_type> URowView(&CurrentRow[LenL+1], LenU);
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
        Teuchos::ArrayView<const local_ordinal_type> IndicesU;
        U_Graph_->getLocalRowView(RowU,IndicesU);
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
      
      TEUCHOS_TEST_FOR_EXCEPTION(Next != i, std::runtime_error,
                         "Ifpack2::IlukGraph::constructFilledGraph: FATAL: U has zero diagonal")

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
  Teuchos::RCP<const map_type> L_DomainMap = OverlapGraph_->getRowMap();
  Teuchos::RCP<const map_type> L_RangeMap  = Graph_->getRangeMap();
  Teuchos::RCP<const map_type> U_DomainMap = Graph_->getDomainMap();
  Teuchos::RCP<const map_type> U_RangeMap  = OverlapGraph_->getRowMap();
  L_Graph_->fillComplete(L_DomainMap, L_RangeMap);//DoOptimizeStorage is default here...
  U_Graph_->fillComplete(U_DomainMap, U_RangeMap);//DoOptimizeStorage is default here...

  Teuchos::reduceAll<int,size_t>(*(L_DomainMap->getComm()), Teuchos::REDUCE_SUM, 1, &NumMyDiagonals_, &NumGlobalDiagonals_);
}

} // namespace Ifpack2

#endif /* IFPACK2_ILUK_GRAPH_HPP */
