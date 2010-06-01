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
  : origMatrix_(origMatrix)
{
}
//=============================================================================
	/*
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrixTransposer(const RowMatrixtransposer& Source)
  :OrigMatrix_(Source.OrigMatrix_)
{
  transposeMatrix_ = Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*Source.transposeMatrix_));
  if (MakeDataContiguous_) transposeMatrix_->MakeDataContiguous();
}*/
//=========================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RowMatrixTransposer(){
}

//=========================================================================
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createTranspose (const OptimizeOption optimizeTranspose, 
						 Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &transposeMatrix){
	if(!transposeMatrix.is_null()){
		TEST_FOR_EXCEPTION(transposeMatrix->isFillComplete(), std::runtime_error,
          Teuchos::typeName(*this) << "::createTranspose(): The matrix to transpose must not be fill complete before calling this function.");
	}
	optimizeTranspose_ = optimizeTranspose;
	if(origMatrix_->isLocallyIndexed()){
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getNodeNumCols(), origMatrix_->getIndexBase(), origMatrix_->getComm()));
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tColMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getNodeNumRows(), origMatrix_->getIndexBase(), origMatrix_->getComm()));
		Teuchos::ArrayRCP<size_t> origColLengths = Teuchos::ArrayRCP<size_t>((LocalOrdinal)origMatrix_->getNodeNumCols(), (LocalOrdinal)0);

		for(size_t i = 0; i < origMatrix_->getNodeNumRows(); ++i){
			Teuchos::ArrayRCP<const LocalOrdinal> indicies;
			Teuchos::ArrayRCP<const Scalar> values;
			origMatrix_->getLocalRowView(i, indicies, values);
			for(LocalOrdinal j=0; j<indicies.size(); ++j){
				origColLengths[indicies[j]] = origColLengths[indicies[j]] +1;
			}
		}
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
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origMatrix_->getGlobalNumCols(), origMatrix_->getIndexBase(), origMatrix_->getComm()));
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
	return 0;
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

