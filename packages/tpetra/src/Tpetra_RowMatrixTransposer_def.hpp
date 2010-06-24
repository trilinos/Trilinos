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
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra{

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrixTransposer(const Teuchos::RCP<const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > origMatrix)
  : origMatrix_(origMatrix), comm_(origMatrix->getComm()), indexBase_(origMatrix_->getIndexBase()) {}
	
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RowMatrixTransposer(){}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createTranspose (const OptimizeOption optimizeTranspose, 
						 Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &transposeMatrix, Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap){
	optimizeTranspose_ = optimizeTranspose;
   	const size_t LST0 = Teuchos::OrdinalTraits<size_t>::zero();
    //const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
   	const LocalOrdinal LO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
   	const Scalar S0 = Teuchos::ScalarTraits<Scalar>::zero();
	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap;
	if(transposeRowMap.is_null()){
		tRowMap = origMatrix_->getDomainMap();
	}
	else{
		tRowMap = transposeRowMap;
	}

	size_t numMyRows = origMatrix_->getNodeNumRows();
	size_t numMyCols = origMatrix_->getNodeNumCols();
	Teuchos::ArrayRCP<size_t> transNumNz = Teuchos::ArrayRCP<size_t>(numMyCols, LST0);
	Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > transIndicies = Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> >(numMyCols);
	Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > transValues = Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >(numMyCols);

	size_t numIndicies;
	size_t maxNumEntries = LST0;
	size_t numEntries;
	for(LocalOrdinal i = LST0; (size_t)i<numMyRows;++i){
		numEntries = origMatrix_->getNumEntriesInLocalRow(i);
		if(numEntries > maxNumEntries){
			maxNumEntries = numEntries;
		}
	}
	Teuchos::ArrayRCP<const LocalOrdinal> indicies;
	Teuchos::ArrayRCP<const Scalar> values;
	for(LocalOrdinal i=LO0; (size_t)i<numMyRows;++i){
		origMatrix_->getLocalRowView(i,indicies,values);
		for(LocalOrdinal j=LO0; j<indicies.size(); ++j){
			++transNumNz[indicies[j]];
		}
	}

	for(LocalOrdinal i=LO0; (size_t)i<numMyCols; i++){
		numIndicies= transNumNz[i];
		if(numIndicies>0){
			transIndicies[i] = Teuchos::ArrayRCP<LocalOrdinal>(numIndicies, LST0);
			transValues[i] = Teuchos::ArrayRCP<Scalar>(numIndicies, S0);
		}
	}

	Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > globalTransIndicies = transIndicies;
	for(LocalOrdinal i=LO0;(size_t)i<numMyCols;++i) transNumNz[i]=0;
	for(LocalOrdinal i=LO0;(size_t)i<numMyRows;++i){
		origMatrix_->getLocalRowView(i,indicies,values);
		GlobalOrdinal ii = origMatrix_->getRowMap()->getGlobalElement(i);
		for(LocalOrdinal j=LO0; j<indicies.size();++j){
			LocalOrdinal transRow= indicies[j];
			size_t loc= transNumNz[transRow];
			globalTransIndicies[transRow][loc] = ii;
			transValues[transRow][loc] = values[j];
			++transNumNz[transRow];
		}
	}

	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transMap = origMatrix_->getColMap();
	CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> tempTransA1(transMap, transNumNz);
	Teuchos::ArrayView<const GlobalOrdinal> transMyGlobalEquations = transMap->getNodeElementList();

   	const LocalOrdinal GO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
	for(GlobalOrdinal i= GO0; (size_t)i<numMyCols; ++i){
		tempTransA1.insertGlobalValues(transMyGlobalEquations[i], transIndicies[i](), transValues[i]());
	}

	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap = origMatrix_->getDomainMap();
	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > rangeMap = origMatrix_->getRangeMap();
	tempTransA1.fillComplete(rangeMap, domainMap, DoNotOptimizeStorage);
	transposeMatrix_ = Teuchos::rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, 0));

	Export<LocalOrdinal, GlobalOrdinal, Node> transExporter(transMap, tRowMap);

	transposeMatrix_->doExport(tempTransA1, transExporter, ADD);
	transposeMatrix_->fillComplete(rangeMap, domainMap, optimizeTranspose_ );
	transposeMatrix = transposeMatrix_;

	
	/*const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transMap = origMatrix_->getColMap();
	transposeMatrix_ = Teuchos::rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, maxNumEntries));
	Teuchos::ArrayView<const GlobalOrdinal> transMyGlobalEquations = transMap->getNodeElementList();
	for(global_size_t i = GST0; i<numMyCols; ++i){
		transposeMatrix_->insertGlobalValues(transMyGlobalEquations[i], transIndicies[i](), transValues[i]());
	}

	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domain_map = origMatrix_->getDomainMap();
	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > range_map = origMatrix_->getRangeMap();

   	const LocalOrdinal GO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
	for(GlobalOrdinal i=GO0; i<=transposeMatrix_->getRowMap()->getMaxAllGlobalIndex(); ++i){
		if(transposeMatrix_->getRowMap()->isNodeGlobalElement(i)){
			global_size_t currentRowLength = transposeMatrix_->getNumEntriesInGlobalRow(i);
			Scalar* valueArrayPtr = new Scalar;
			GlobalOrdinal* indexArrayPtr = new GlobalOrdinal;
			Teuchos::ArrayView<Scalar> currentRowValues(valueArrayPtr, currentRowLength);
			Teuchos::ArrayView<GlobalOrdinal> currentRowIndices(indexArrayPtr, currentRowLength);
			transposeMatrix_->getGlobalRowCopy(i, currentRowIndices, currentRowValues, currentRowLength);
			std::cout << "Node: " << transposeMatrix_->getComm()->getRank() << " Global Row: " << i << " Indicies: " << currentRowIndices << " Values: " << currentRowValues << "\n";
			delete valueArrayPtr;
			delete indexArrayPtr;
		}
	}

	std::cout << "right before fill complete\n";
	Export transExporter(origMatrix_->getColMap(), tRowMap());
	transposeMatrix_->doExport(transposeMatrix_,
	transposeMatrix_->fillComplete(origMatrix_->getDomainMap(),origMatrix_->getRangeMap(), optimizeTranspose_);
	std::cout << "right after fill complete\n";
	transposeMatrix = transposeMatrix_;*/
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

