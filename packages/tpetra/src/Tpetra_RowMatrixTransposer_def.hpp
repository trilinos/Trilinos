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
	for(size_t i = LST0; i<numMyRows;++i){
		numEntries = origMatrix_->getNumEntriesInLocalRow(numMyRows);
		if(numEntries > maxNumEntries){
			maxNumEntries = numEntries;
		}
	}
	Teuchos::ArrayRCP<const LocalOrdinal> indicies;// = Teuchos::ArrayRCP<LocalOrdinal>(maxNumEntries);
	Teuchos::ArrayRCP<const Scalar> values;// = Teuchos::ArrayRCP<Scalar>(maxNumEntries);
	for(LocalOrdinal i=LO0; (size_t)i<numMyRows;++i){
		origMatrix_->getLocalRowView(i,indicies,values);
		for(LocalOrdinal j=LO0; j<indicies.size(); ++j){
			++transNumNz[indicies[j]];
		}
	}
	for(LocalOrdinal i=LO0; (size_t)i<numMyCols; i++){
		numIndicies= transNumNz[i];
		if(numIndicies>0){
			transIndicies[i] = Teuchos::ArrayRCP<LocalOrdinal>(numIndicies);
			transValues[i] = Teuchos::ArrayRCP<Scalar>(numIndicies);
		}
	}

	Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > globalTransIndicies = transIndicies;
	for(LocalOrdinal i=LO0;(size_t)i<numMyCols;++i) transNumNz[i]=0;
	for(LocalOrdinal i=LO0;(size_t)i<numMyRows;++i){
		origMatrix_->getLocalRowView(i,indicies,values);
		GlobalOrdinal ii = origMatrix_->getRowMap()->getGlobalElement(i);
		for(LocalOrdinal j=LO0;(size_t)j<numIndicies;++j){
			LocalOrdinal transRow= indicies[j];
			size_t loc= transNumNz[transRow];
			globalTransIndicies[transRow][loc] = ii;
			transValues[transRow][loc] = values[j];
			++transNumNz[transRow];
		}
	}

	const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transMap = origMatrix_->getColMap();
	transposeMatrix_ = Teuchos::rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(transMap, transNumNz));
	Teuchos::ArrayView<const GlobalOrdinal> transMyGlobalEquations = transMap->getNodeElementList();
	for(global_size_t i = GST0; i<numMyCols; ++i){
		transposeMatrix_->insertGlobalValues(transMyGlobalEquations[i], transIndicies[i](), transValues[i]());
	}

	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domain_map = origMatrix_->getDomainMap();
	Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > range_map = origMatrix_->getRangeMap();

	transposeMatrix_->fillComplete(domain_map,range_map, optimizeTranspose_);
	transposeMatrix = transposeMatrix_;



/*	if(!transposeMatrix.is_null()){
		TEST_FOR_EXCEPTION(transposeMatrix->isFillComplete(), std::runtime_error,
          Teuchos::typeName(*this) << "::createTranspose(): The matrix to transpose must not be fill complete before calling this function.");
	}
	optimizeTranspose_ = optimizeTranspose;
	int myRank = comm_->getRank(); 
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const size_t LST0 = Teuchos::OrdinalTraits<size_t>::zero();
    //const Scalar SST0 = Teuchos::ScalarTraits<Scalar>::zero();
	global_size_t origGlobalNumCols = (global_size_t)origMatrix_->getDomainMap()->getGlobalNumElements();
	global_size_t origGlobalNumRows = (global_size_t)origMatrix_->getRangeMap()->getGlobalNumElements();
	if(origMatrix_->isLocallyIndexed()){
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap;
		
		if(transposeRowMap.is_null()){
			tRowMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumCols, indexBase_, comm_));
		}
		else{
			tRowMap = transposeRowMap;
		}
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tColMap = Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> >(new Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumRows, indexBase_, comm_));
		

//		Teuchos::ArrayRCP<global_size_t> origColLengths = Teuchos::ArrayRCP<global_size_t>(origGlobalNumCols, GST0);
		Teuchos::RCP<Teuchos::SerialComm<GlobalOrdinal> > sComm = Teuchos::rcp(new Teuchos::SerialComm<GlobalOrdinal>());
		Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > colLengthsMap = Teuchos::rcp(new Map<LocalOrdinal, GlobalOrdinal, Node>(origGlobalNumCols, indexBase_, sComm));
		Vector<global_size_t, LocalOrdinal, GlobalOrdinal, Node> origColLengthsVector(colLengthsMap);

		GlobalOrdinal minColIndex = origMatrix_->getColMap()->getMinGlobalIndex();
		Teuchos::ArrayRCP<global_size_t> origColLengths = Teuchos::ArrayRCP<global_size_t>(origGlobalNumCols, GST0);
		for(size_t i = LST0; i < origMatrix_->getNodeNumRows(); ++i){
			Teuchos::ArrayRCP<const LocalOrdinal> indicies;
			Teuchos::ArrayRCP<const Scalar> values;
			origMatrix_->getLocalRowView(i, indicies, values);
			for(LocalOrdinal j=0; j<indicies.size(); ++j){
				origColLengthsVector.sumIntoLocalValue(indicies[j]+minColIndex, GST1);
			}
		}
		std::cout << "gonna print now on " << myRank <<"\n";
		Teuchos::ArrayView<global_size_t> a = Teuchos::ArrayView<global_size_t>(new global_size_t, origGlobalNumCols);
		origColLengthsVector.get1dCopy(a);
		std::cout << "array " << a << "\n";

		Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal,Node> > sourceMap = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(origGlobalNumCols, indexBase_, comm_));
		global_size_t targetNumElements;
		if(myRank==0) targetNumElements=origGlobalNumCols;
		else targetNumElements=0;

		Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal,Node> > targetMap = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(origGlobalNumCols, indexBase_, comm_));

		
		Teuchos::ArrayRCP<global_size_t> masterOrigColLengths = Teuchos::ArrayRCP<global_size_t>(origGlobalNumCols, GST0);
		Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, origColLengths.size(), origColLengths.getRawPtr(), masterOrigColLengths.getRawPtr());
		size_t maxOrigColLength = masterOrigColLengths[0];
		for(size_t i = 1; i<(size_t)origColLengths.size(); i++){
			if(masterOrigColLengths[i]>maxOrigColLength){
				maxOrigColLength = masterOrigColLengths[i];
			}
		}

		transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tRowMap, maxOrigColLength));
		
		for(global_size_t i=0; i<transposeMatrix_->getGlobalNumRows(); ++i){
			Teuchos::ArrayRCP<GlobalOrdinal> rowToInsertIndicies = Teuchos::ArrayRCP<GlobalOrdinal>(masterOrigColLengths[i], GST0);
			Teuchos::ArrayRCP<Scalar> rowToInsertVals=Teuchos::ArrayRCP<Scalar>(masterOrigColLengths[i], SST0);
			GlobalOrdinal toInsertCounter =0;
			for(global_size_t j = 0; j < origMatrix_->getNodeNumRows(); ++j){
				Teuchos::ArrayRCP<const GlobalOrdinal> indicies;
				Teuchos::ArrayRCP<const Scalar> values;
				origMatrix_->getLocalRowView(j, indicies, values);
				for(size_t k=LST0; k<(size_t)indicies.size(); k++){
					if((size_t)indicies[k] ==i){
						rowToInsertIndicies[toInsertCounter] = j+minColIndex;
						rowToInsertVals[toInsertCounter] = values[k];
						toInsertCounter++;
						break;
					}
				}
				if(toInsertCounter >= (GlobalOrdinal)masterOrigColLengths[i]){
					break;
				}
			}
			transposeMatrix_->insertGlobalValues(i, rowToInsertIndicies(), rowToInsertVals());
		}
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

	Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tDomainMap = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(origGlobalNumCols,origMatrix_->getIndexBase(), comm_));
	Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > tRangeMap = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(origGlobalNumRows,origMatrix_->getIndexBase(), comm_));
	if(optimizeTranspose == DoOptimizeStorage){
		transposeMatrix_->fillComplete(DoOptimizeStorage);
		//transposeMatrix_->fillComplete(tDomainMap, tRangeMap, DoOptimizeStorage);
	}
	else{
		transposeMatrix_->fillComplete(DoNotOptimizeStorage);
		//transposeMatrix_->fillComplete(tDomainMap, tRangeMap, DoNotOptimizeStorage);
	}
	transposeMatrix = transposeMatrix_;
	*/
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

