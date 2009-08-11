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

#ifndef TIFPACK_CONSTRUCTLEVELFILLGRAPH_HPP
#define TIFPACK_CONSTRUCTLEVELFILLGRAPH_HPP

#include "Tifpack_ConfigDefs.hpp"
#include <Teuchos_Describable.hpp>
#include "Tpetra_RowGraph.hpp" 
#include "Tpetra_CrsGraph.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
	class ParameterList;
}

namespace Tifpack {
	
	
	//! ConstructLevelFillGraph: A function for constructing level filled graphs for use with ILU(k) class preconditioners.
	
	/*! The IlukGraph class enable the construction matrix graphs using level-fill algorithms.  The only function
	 required for construction is an ExtractRowView capability, i.e., the matrix that is passed in to the constructor
	 must implement the CrsGraph interface defined in CrsMatrix.hpp 
	 
	 
	 <b>Constructing IlukGraph objects</b>
	 
	 Constructing IlukGraph objects is usually a two step process of passing in a CrsGraph object and 
	 an integer indicating the desired level of fill and then calling the ConstructFilledGraph function to complete the
	 process.  This allows warning error codes to be returned to the calling routine.
	 
	 It is worth noting that an IlukGraph object has two Tpetra::CrsGraph objects containing L and U, 
	 the graphs for the lower and upper triangular parts of the ILU(k) graph.
	 Thus, it is possible to manually insert and delete graph entries in L and U via the Tpetra_CrsGraph
	 InsertIndices and RemoveIndices functions.  However, in this case FillComplete must be
	 called before the graph is used for subsequent operations.
	 
	 \param userGraph (In) - An existing CrsGraph.  This object must implement the RowGraph functions that provide graph dimension and pattern information.
	 \param levelFill (In) - The level of fill to compute via ILU(k) algorithm.
	 \param graphL (Out) - Lower triangular graph with the level k fill pattern.
	 \param graphU (Out) - Upper triangular graph with the level k fill pattern.
	 */    
	template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal>
	void ConstructLevelFillGraph(const RowGraph<LocalOrdinal,GlobalOrdinal> & userGraph, 
								 Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal>> graphL, Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal>> graphU) {
	
	LocalOrdinal i, j;
	LocalOrdinal numIn, numL, numU;
	LocalOrdinal numMyDiagonals = 0;
	bool diagFound;
	
	
	
	graphL = Teuchos::rcp(new CrsGraph(userGraph.getRowMap(), userGraph.getRowMap(),  0));
	graphU = Teuchos::rcp(new CrsGraph(userGraph.getRowMap(), userGraph.getRowMap(),  0));
	
	
	// Get Maximun Row length
	int MaxNumIndices = userGraph.getMaxNumIndices();
	
	Teuchos::ArrayView<LocalOrdinal> userGraphRow;
	vector<LocalOrdinal> L(MaxNumIndices);
	vector<LocalOrdinal> U(MaxNumIndices);
	LocalOrdinal * pL = &L[0];
	LocalOrdinal * pU = &U[0];
	
	LocalOrdinal numMyRows = userGraph.getNumMyRows();
	
	// First we copy the user's graph into L and U, regardless of fill level
	
	for (i=0; i< numMyRows; ++i) {
		
		
		userGraph.extractMyRowConstView(i, userGraphRow); // Get Indices
		
		
		// Split into L and U (we don't assume that indices are ordered).
		
		numL = 0; 
		numU = 0; 
		diagFound = false;
		numIn = userGraphRow.size();
		
		for (j=0; j< numIn; ++j) {
			int k = userGraphRow[j];
			
			if (k<numMyRows) { // Ignore column elements that are not in the square matrix
				if (k==i) {DiagFound = true;}	
				else if (k < i) {pL[numL++] = k;}
				else {pU[numU++] = k;}
			}
		}
		
		// Check in things for this row of L and U
		
		if (diagFound) numMyDiagonals++;
		if (NumL>0) {graphL->insertMyIndices(i, ArrayView(pL, numL));}
		if (NumU>0) {graphU->insertMyIndices(i, ArrayView(pU, numU));}
	}
	
	//	if (LevelFill_ > 0) {
	//		
	//		// Complete Fill steps
	//		graphL->fillComplete(userGraph.getRowMap(), userGraph.getRangeMap());
	//		graphU->fillComplete(userGraph.getDomainMap(), userGraph.getRowMap());
	//		
	//		// At this point L_Graph and U_Graph are filled with the pattern of input graph, 
	//		// sorted and have redundant indices (if any) removed.  Indices are zero based.
	//		// LevelFill is greater than zero, so continue...
	//		
	//		vector<vector<LocalOrdinal> > levels(numMyRows);
	//		vector<LocalOrdinal> linkList(numMyRows);
	//		vector<LocalOrdinal> currentLevel(numMyRows);
	//		vector<LocalOrdinal> currentRow(numMyRows);
	//		vector<LocalOrdinal> levelsRowU(numMyRows);
	//		
	//		for (i=0; i<numMyRows; ++i) {
	//			int First, Next;
	//			
	//			// copy column indices of row into workspace and sort them
	//			
	//			int LenL = L_Graph_->NumMyIndices(i);
	//			int LenU = U_Graph_->NumMyIndices(i);
	//			int Len = LenL + LenU + 1;
	//			
	//			EPETRA_CHK_ERR(L_Graph_->ExtractMyRowCopy(i, LenL, LenL, &CurrentRow[0]));      // Get L Indices
	//			CurrentRow[LenL] = i;                                     // Put in Diagonal
	//			int ierr1 = 0;
	//			if (LenU) {
	//				// Get U Indices
	//				ierr1 = U_Graph_->ExtractMyRowCopy(i, LenU, LenU, &CurrentRow[LenL+1]);
	//			}
	//			if (ierr1!=0) {
	//				cout << "ierr1 = "<< ierr1 << endl;
	//				cout << "i = " << i << endl;
	//				cout << "NumMyBlockRows_ = " << graphU->NumMyBlockRows() << endl;
	//			}
	//				
	//			// Construct linked list for current row
	//				
	//			for (j=0; j<Len-1; j++) {
	//				LinkList[CurrentRow[j]] = CurrentRow[j+1];
	//				CurrentLevel[CurrentRow[j]] = 0;
	//			}
	//				
	//				LinkList[CurrentRow[Len-1]] = NumMyBlockRows_;
	//				CurrentLevel[CurrentRow[Len-1]] = 0;
	//				
	//				// Merge List with rows in U
	//				
	//				First = CurrentRow[0];
	//				Next = First;
	//				while (Next < i)
	//				{
	//					int PrevInList = Next;
	//					int NextInList = LinkList[Next];
	//					int RowU = Next;
	//					int LengthRowU;
	//					int * IndicesU;
	//					// Get Indices for this row of U
	//					EPETRA_CHK_ERR(U_Graph_->ExtractMyRowView(RowU, LengthRowU, IndicesU));
	//					
	//					int ii;
	//					
	//					// Scan RowU
	//					
	//					for (ii=0; ii<LengthRowU; /*nop*/) {
	//						int CurInList = IndicesU[ii];
	//						if (CurInList < NextInList) {
	//							// new fill-in
	//							int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
	//							if (NewLevel <= LevelFill_)
	//							{
	//								LinkList[PrevInList]  = CurInList;
	//								LinkList[CurInList] = NextInList;
	//								PrevInList = CurInList;
	//								CurrentLevel[CurInList] = NewLevel;
	//							}
	//							ii++;
	//						}
	//						else if (CurInList == NextInList)
	//						{
	//							PrevInList = NextInList;
	//							NextInList = LinkList[PrevInList];
	//							int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
	//							CurrentLevel[CurInList] = EPETRA_MIN(CurrentLevel[CurInList], NewLevel);
	//							ii++;
	//						}
	//						else // (CurInList > NextInList)
	//						{
	//							PrevInList = NextInList;
	//							NextInList = LinkList[PrevInList];
	//						}
	//					}
	//					Next = LinkList[Next];
	//				}
	//				
	//				// Put pattern into L and U
	//				
	//				LenL = 0;
	//				
	//				Next = First;
	//				
	//				// Lower
	//				
	//				while (Next < i) {	  
	//					CurrentRow[LenL++] = Next;
	//					Next = LinkList[Next];
	//				}
	//				
	//				EPETRA_CHK_ERR(L_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
	//				int ierr11 = L_Graph_->InsertMyIndices(i, LenL, &CurrentRow[0]);
	//				if (ierr11 < 0) EPETRA_CHK_ERR(ierr1);
	//				
	//				// Diagonal
	//				
	//				if (Next != i) return(-2); // Fatal:  U has zero diagonal.
	//				else {
	//					LevelsRowU[0] = CurrentLevel[Next];
	//					Next = LinkList[Next];
	//				}
	//				
	//				// Upper
	//				
	//				LenU = 0;
	//				
	//				while (Next < NumMyBlockRows_) // Should be "Next < NumMyBlockRows_"?
	//				{
	//					LevelsRowU[LenU+1] = CurrentLevel[Next];
	//					CurrentRow[LenU++] = Next;
	//					Next = LinkList[Next];
	//				}
	//				
	//				EPETRA_CHK_ERR(U_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
	//				int ierr2 = U_Graph_->InsertMyIndices(i, LenU, &CurrentRow[0]);
	//				if (ierr2<0) EPETRA_CHK_ERR(ierr2);
	//				
	//				// Allocate and fill Level info for this row
	//				Levels[i] = vector<int>(LenU+1);
	//				for (int jj=0; jj<LenU+1; jj++) Levels[i][jj] = LevelsRowU[jj];
	//				
	//				}
	//				}    
	
	// Complete Fill steps
	graphL->fillComplete(userGraph.getRowMap(), userGraph.getRangeMap());
	graphU->fillComplete(userGraph.getDomainMap(), userGraph.getRowMap());
	
	
	// Optimize graph storage
	
	graphL->optimizeStorage();
	graphU->optimizeStorage();
	
	// Compute global quantities
	
	GlobalOrdinal numGlobalDiagonals = 0;
	
	graphL->getComm().sumAll(&numMyDiagonals, &numGlobalDiagonals, 1));
	
	return;
}

} // namespace Tifpack

#endif /* TIFPACK_CONSTRUCTLEVELFILLGRAPH_HPP */
