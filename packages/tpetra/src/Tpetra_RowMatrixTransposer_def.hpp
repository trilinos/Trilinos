//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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

#include "Tpetra_Map.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra{

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrixTransposer(const Teuchos::RCP<const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > origMatrix)
  : origMatrix_(origMatrix), comm_(origMatrix->getComm()) {}
	
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RowMatrixTransposer(){}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createTranspose (const OptimizeOption optimizeTranspose, 
						 Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &transposeMatrix, Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap){
	if(!transposeMatrix.is_null()){
		TEST_FOR_EXCEPTION(transposeMatrix->isFillComplete(), std::runtime_error,
          Teuchos::typeName(*this) << "::createTranspose(): The matrix to transpose must not be fill complete before calling this function.");
	}
	optimizeTranspose_ = optimizeTranspose;
	int myRank = comm_->getRank(); 
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const size_t LST0 = Teuchos::OrdinalTraits<size_t>::zero();
	if(origMatrix_->isLocallyIndexed()){
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap;
		
		global_size_t origGlobalNumCols = (global_size_t)origMatrix_->getDomainMap()->getGlobalNumElements();
		global_size_t origGlobalNumRows = (global_size_t)origMatrix_->getRangeMap()->getGlobalNumElements();
		if(transposeRowMap.is_null()){
			//tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getNodeNumCols(), origMatrix_->getIndexBase(), comm_));
			tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumCols, origMatrix_->getIndexBase(), comm_));
		}
		else{
			tRowMap = transposeRowMap;
		}
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tColMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumRows, origMatrix_->getIndexBase(), comm_));

		//Teuchos::ArrayRCP<size_t> origColLengths = Teuchos::ArrayRCP<size_t>((LocalOrdinal)origMatrix_->getNodeNumCols(), (LocalOrdinal)0);
		//Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > sourceMap = Teuchos::rcp(Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumCols, origGlobalNumCols, GST0, comm_));
		//Vector<global_size_t> origColLengths = Vector<global_size_t>(sourceMap);
		Teuchos::ArrayRCP<global_size_t> origColLengths = Teuchos::ArrayRCP<global_size_t>(origGlobalNumCols, GST0);

		for(size_t i = LST0; i < origMatrix_->getNodeNumRows(); ++i){
			Teuchos::ArrayRCP<const LocalOrdinal> indicies;
			Teuchos::ArrayRCP<const Scalar> values;
			origMatrix_->getLocalRowView(i, indicies, values);
			GlobalOrdinal minColIndex = origMatrix_->getColMap()->getMinGlobalIndex();
			for(LocalOrdinal j=0; j<indicies.size(); ++j){
				origColLengths[indicies[j]+minColIndex] = origColLengths[indicies[j]+minColIndex] +GST1;
			}
		}

		Teuchos::reduceAll<typename Teuchos::ArrayRCP<global_size_t>, global_size_t>(*comm_, Teuchos::REDUCE_SUM, origColLengths.size(), origColLengths.getRawPtr(), origColLengths.getRawPtr());
		if(myRank == 0){
			std::cout << "original col lengths " << origColLengths() << "\n";
		}
		/*global_size_t targetNumElements;
		if(myRank ==0){
			targetNumElements = origGlobalNumCols;
		}
		else{
			targetNumElements = GST0;
		}
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap = Teuchos::rcp(Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumCols, targetNumElements, GST0, comm_));
		Import importer<LocalOrdinal, GlobalOrdinal, Node>(sourceMap, targetMap);
		Vector target(targetMap);*/



		size_t maxOrigColLength = origColLengths[0];
		for(size_t i = 1; i<(size_t)origColLengths.size(); i++){
			if(origColLengths[i]>maxOrigColLength){
				maxOrigColLength = origColLengths[i];
			}
		}
		transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, tColMap, maxOrigColLength));
		for(size_t i=0; i<transposeMatrix_->getNodeNumRows(); ++i){
			Teuchos::ArrayRCP<LocalOrdinal> rowToInsertIndicies = Teuchos::ArrayRCP<LocalOrdinal>(origColLengths[i], (LocalOrdinal)0);
			Teuchos::ArrayRCP<Scalar> rowToInsertVals=Teuchos::ArrayRCP<Scalar>(origColLengths[i], (Scalar)0);
			LocalOrdinal toInsertCounter =0;
			for(size_t j = 0; j < origMatrix_->getNodeNumRows(); ++j){
				Teuchos::ArrayRCP<const LocalOrdinal> indicies;
				Teuchos::ArrayRCP<const Scalar> values;
				origMatrix_->getLocalRowView(j, indicies, values);
				for(size_t k=0; k<(size_t)indicies.size(); k++){
					if((size_t)indicies[k] ==i){
						rowToInsertIndicies[toInsertCounter] = j;
						rowToInsertVals[toInsertCounter] = values[k];
						toInsertCounter++;
						break;
					}
				}
				if(toInsertCounter >= (LocalOrdinal)origColLengths[i]){
					break;
				}
			}
			transposeMatrix_->insertLocalValues(i, rowToInsertIndicies(), rowToInsertVals());
		}
	}
	else{
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap;
		if(transposeRowMap.is_null()){
			tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getGlobalNumCols(), origMatrix_->getIndexBase(), origMatrix_->getComm()));
		}
		else{
			tRowMap = transposeRowMap;
		}
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tColMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getGlobalNumRows(), origMatrix_->getIndexBase(), origMatrix_->getComm()));
		Teuchos::ArrayRCP<global_size_t> origColLengths = Teuchos::ArrayRCP<global_size_t>((GlobalOrdinal)origMatrix_->getGlobalNumCols(), (global_size_t)0);

		for(global_size_t i = 0; i < origMatrix_->getGlobalNumRows(); ++i){
			Teuchos::ArrayRCP<const GlobalOrdinal> indicies;
			Teuchos::ArrayRCP<const Scalar> values;
			origMatrix_->getGlobalRowView(i, indicies, values);
			for(GlobalOrdinal j=0; j<indicies.size(); ++j){
				origColLengths[indicies[j]] = origColLengths[indicies[j]] +1;
			}
		}
		global_size_t maxOrigColLength = origColLengths[0];
		for(global_size_t i = 1; i<(global_size_t)origColLengths.size(); i++){
			if(origColLengths[i]>maxOrigColLength){
				maxOrigColLength = origColLengths[i];
			}
		}
		transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, tColMap, maxOrigColLength));
		for(global_size_t i=0; i<transposeMatrix_->getGlobalNumRows(); ++i){
			Teuchos::ArrayRCP<GlobalOrdinal> rowToInsertIndicies = Teuchos::ArrayRCP<GlobalOrdinal>(origColLengths[i], (GlobalOrdinal)0);
			Teuchos::ArrayRCP<Scalar> rowToInsertVals=Teuchos::ArrayRCP<Scalar>(origColLengths[i], (Scalar)0);
			GlobalOrdinal toInsertCounter =0;
			for(global_size_t j = 0; j < origMatrix_->getGlobalNumRows(); ++j){
				Teuchos::ArrayRCP<const GlobalOrdinal> indicies;
				Teuchos::ArrayRCP<const Scalar> values;
				origMatrix_->getGlobalRowView(j, indicies, values);
				for(global_size_t k=0; k<(global_size_t)indicies.size(); k++){
					if((global_size_t)indicies[k] ==i){
						rowToInsertIndicies[toInsertCounter] = j;
						rowToInsertVals[toInsertCounter] = values[k];
						toInsertCounter++;
						break;
					}
				}
				if(toInsertCounter >= (GlobalOrdinal)origColLengths[i]){
					break;
				}
			}
			transposeMatrix_->insertGlobalValues(i, rowToInsertIndicies(), rowToInsertVals());
		}
	}
	if(optimizeTranspose == DoOptimizeStorage){
		transposeMatrix_->fillComplete(DoOptimizeStorage);
	}
	else{
		transposeMatrix_->fillComplete(DoNotOptimizeStorage);
	}
	transposeMatrix = transposeMatrix_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(LO,GO,NODE) \
  \
  template class RowMatrixTransposer< LO , GO , NODE >;


}

