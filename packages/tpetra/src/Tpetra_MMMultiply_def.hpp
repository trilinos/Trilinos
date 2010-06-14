//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef TPETRA_MMMULTIPLY_DEF_HPP
#define TPETRA_MMMULTIPLY_DEF_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_map.hpp"
#include <algorithm>
#include <vector>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#ifdef DOXYGEN_USE_ONLY
//#include "Tpetra_MMMultiply_decl.hpp"
#endif

/*! \file Tpetra_MMMultiply_def.hpp 

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */
namespace Tpetra {




	template<class GlobalOrdinal, class Scalar>
	void mmInsertValue(std::vector<GlobalOrdinal>& indexArray, std::vector<Scalar>& valueArray, GlobalOrdinal indexToInsert, Scalar valueToInsert){
		typename std::vector<GlobalOrdinal>::iterator it=find(indexArray.begin(), indexArray.end(), indexToInsert);
		if(it != indexArray.end()){
			valueArray[it-indexArray.begin()] += valueToInsert;
		}
		else{
			it=indexArray.begin();
			while(it !=indexArray.end() && (*it) < indexToInsert)  ++it;
			if(it!=indexArray.end()){
				indexArray.insert(it, indexToInsert);
				valueArray.insert(valueArray.begin() + (it-indexArray.begin()), valueToInsert);
			}
			else{
				indexArray.push_back(indexToInsert);
				valueArray.push_back(valueToInsert);
			}
		}
	}


	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	MatrixMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MatrixMatrixMultiply( 
		const Teuchos::RCP<const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A,
		const Teuchos::RCP<const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > B, 
		Teuchos::RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > C)
	:matrixA(A), matrixB(B), matrixC(C){}

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	int MatrixMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::multiply(){
		TEST_FOR_EXCEPTION(matrixA->getGlobalNumRows() != matrixB->getGlobalNumCols(), std::runtime_error,
			"Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
			"a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
			"\n\n" << Teuchos::typeName(*this) << "::multiply(): Matrix A must have the same number of rows as Matrix B has columns.");

		TEST_FOR_EXCEPTION(matrixA->getGlobalNumRows() != matrixC->getGlobalNumRows(), std::runtime_error,
			"Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
			"a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
			"\n\n" << Teuchos::typeName(*this) << "::multiply(): Matrix A must have the same number of rows as Matrix C.");

		if(matrixC->hasColMap()){
			TEST_FOR_EXCEPTION(matrixB->getGlobalNumCols() != matrixC->getGlobalNumCols(), std::runtime_error,
				"Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
				"a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
				"\n\n" << Teuchos::typeName(*this) << "::multiply(): Matrix B must have the same number of rows as Matrix C.");
		}

		TEST_FOR_EXCEPTION(!matrixA->isFillComplete(), std::runtime_error,
			"Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
			"a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
			"\n\n" << Teuchos::typeName(*this) << "::multiply(): Matrix A is not fill complete.");

		TEST_FOR_EXCEPTION(!matrixB->isFillComplete(), std::runtime_error,
			"Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
			"a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
			"\n\n" << Teuchos::typeName(*this) << "::multiply(): Matrix B is not fill complete.");



		//const LocalOrdinal LO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
		//const GlobalOrdinal GO0 = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
		Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tB;
		RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(matrixB);
		transposer.createTranspose(DoOptimizeStorage, tB, matrixA->getRowMap());
		Teuchos::ArrayRCP<const LocalOrdinal> aCurRowIndices;
		Teuchos::ArrayRCP<const Scalar> aCurRowValues;
		Teuchos::ArrayRCP<const LocalOrdinal> bCurRowIndices;
		Teuchos::ArrayRCP<const Scalar> bCurRowValues;
		std::vector<GlobalOrdinal> cCurRowIndices;
		std::vector<Scalar> cCurRowValues;
		typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator findResult;

		Teuchos::ArrayRCP<const GlobalOrdinal> importedBRowIndicies;
		Teuchos::ArrayRCP<Teuchos::ArrayRCP<const GlobalOrdinal> > importedBIndicies;
		Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > importedBValues;

		getAndSendBRows(tB, importedBRowIndicies, importedBIndicies, importedBValues);
/*
		for(LocalOrdinal i = LO0; (size_t)i<matrixA->getNodeNumRows(); ++i){
		//	std::cout << "Row: " << i << "\n";
			cCurRowIndices = std::vector<GlobalOrdinal>(0);
			cCurRowValues = std::vector<Scalar>(0);
			matrixA->getLocalRowView(i, aCurRowIndices, aCurRowValues);
			for(GlobalOrdinal j= GO0; (global_size_t)j<tB->getGlobalNumRows(); ++j){
				//std::cout << "\tCol: " << j << "\n";
				if(tB->getRowMap()->isNodeGlobalElement(j)){
					tB->getLocalRowView(tB->getRowMap()->getLocalElement(j), bCurRowIndices, bCurRowValues);
				}
				else{

				}
					
				for(LocalOrdinal k =LO0; k<aCurRowIndices.size(); ++k){
					findResult = std::find(bCurRowIndices.begin(), bCurRowIndices.end(), aCurRowIndices[k]);
					if(findResult != bCurRowIndices.end()){
						LocalOrdinal actualBValueIndex = findResult-bCurRowIndices.begin();
				//		std::cout << "\t\t" << aCurRowValues[k] << " * " << bCurRowValues[actualBValueIndex] << "\n";
						mmInsertValue(cCurRowIndices, cCurRowValues, (GlobalOrdinal)j, aCurRowValues[k]*bCurRowValues[actualBValueIndex]);
					}
				}
			}
			std::cout << "Going to insert: \n";
			std::cout << "Indicies: \n";
			for(unsigned int l=0; l < cCurRowIndices.size(); ++l){
				std::cout << cCurRowIndices[l] << "\n";
			}
			std::cout << "With values: " << Teuchos::ArrayView<Scalar>(cCurRowValues) << "\n";
			std::cout << "Into row: " << i << "\n";
			matrixC->insertGlobalValues((GlobalOrdinal)i, Teuchos::ArrayView<GlobalOrdinal>(cCurRowIndices), Teuchos::ArrayView<Scalar>(cCurRowValues));
		}*/
		matrixC->fillComplete();
		return 0;
	}

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	MatrixMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~MatrixMatrixMultiply(){}

	template<class LocalOrdinal, class GlobalOrdinal, class Node>
	Teuchos::ArrayRCP<const GlobalOrdinal> convertLocalIndicesToGlobal(Teuchos::ArrayRCP<const LocalOrdinal> localIndicies, Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map){
		Teuchos::ArrayRCP<GlobalOrdinal> toReturn(localIndicies.size());
		const GlobalOrdinal GO0 = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
		typename Teuchos::ArrayRCP<const LocalOrdinal>::const_iterator it;
		it = localIndicies.begin();
		for(GlobalOrdinal i= GO0; i<toReturn.size(); ++i, ++it){
			toReturn[i] = map->getGlobalElement(*it);
		}
		return toReturn;
	}

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void MatrixMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getAndSendBRows(
		Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &tB,
		Teuchos::ArrayRCP<const GlobalOrdinal> &importedBRowIndicies,
		Teuchos::ArrayRCP<Teuchos::ArrayRCP<const GlobalOrdinal> > &importedBIndicies,
		Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &importedBValues){


		const LocalOrdinal LO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
		const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
		//Calculate rows that need to be sent and recieved

		//A map of the rows that this node owns and to whom they need to be sent
		Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> > tBGraph = tB->getGraph();
		Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> > mAGraph = matrixA->getGraph();
		//Maps a row id to a list of NodeIDs that the row needs to be sent To
		typedef Teuchos::map<GlobalOrdinal, Teuchos::ArrayView<int> > DestinationMap;
		DestinationMap mySendTos;
		Teuchos::ArrayRCP<const GlobalOrdinal> curRowIndicies;
		for(LocalOrdinal i = LO0; (size_t)i< tBGraph->getNodeNumRows(); i++){
			curRowIndicies = convertLocalIndicesToGlobal(tBGraph->getLocalRowView(i), tBGraph->getRowMap());
			if((global_size_t)curRowIndicies.size() != GST0){
				mySendTos.insert(std::pair<GlobalOrdinal, const Teuchos::ArrayView<int> >(i, Teuchos::ArrayView<int>()));
				const Teuchos::ArrayView<int> nodesToSendTo;
				mAGraph->getColMap()->getRemoteIndexList(curRowIndicies(), mySendTos[i]);
			}
		}

		DestinationMap myRecieves;
		for(LocalOrdinal i = LO0; (size_t)i< mAGraph->getNodeNumRows(); i++){
			curRowIndicies = convertLocalIndicesToGlobal(mAGraph->getLocalRowView(i), mAGraph->getRowMap());
			if((global_size_t)curRowIndicies.size() != GST0){
				myRecieves.insert(std::pair<GlobalOrdinal, const Teuchos::ArrayView<int> >(i, Teuchos::ArrayView<int>()));
				const Teuchos::ArrayView<int> nodesToSendTo;
				tBGraph->getColMap()->getRemoteIndexList(curRowIndicies(), mySendTos[i]);
			}
		}

		if(matrixA->getComm()->getRank()==0){
			std::cout << "Sent to map for 0 \n";
			for(typename DestinationMap::iterator it = mySendTos.begin(); it!=mySendTos.end(); it++){
				std::cout << "Row: " << it->first << "ids to send to: " << it->second << "\n";

			}

		}
			
	}


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MATRIXMATRIXMULTIPLY_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class MatrixMatrixMultiply< SCALAR , LO , GO , NODE >;

}
#endif // TPETRA_MMMULTIPLY_DEF_HPP
