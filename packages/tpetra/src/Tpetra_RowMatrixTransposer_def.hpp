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
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
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
	
	const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transMap = origMatrix_->getColMap();
	transposeMatrix_ = Teuchos::rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, maxNumEntries));
	Teuchos::ArrayView<const GlobalOrdinal> transMyGlobalEquations = transMap->getNodeElementList();
	for(global_size_t i = GST0; i<numMyCols; ++i){
		transposeMatrix_->insertGlobalValues(transMyGlobalEquations[i], transIndicies[i](), transValues[i]());
	}

	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domain_map = origMatrix_->getDomainMap();
	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > range_map = origMatrix_->getRangeMap();

	transposeMatrix_->fillComplete(domain_map,range_map, optimizeTranspose_);
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

