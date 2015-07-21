// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_ReorderADValues_Evaluator_impl_hpp__
#define __Panzer_ReorderADValues_Evaluator_impl_hpp__


#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Phalanx_DataLayout.hpp"

template<typename EvalT,typename TRAITS>
panzer::ReorderADValues_Evaluator<EvalT, TRAITS>::
ReorderADValues_Evaluator(const std::string & outPrefix,
                          const std::vector<std::string> & inFieldNames,
                          const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                          const std::string & elementBlock,
                          const UniqueGlobalIndexerBase & indexerSrc,
                          const UniqueGlobalIndexerBase & indexerDest)
{ 
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
    outFields_.push_back(PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]));

    // tell the field manager that we depend on this field
    this->addDependentField(inFields_[eq]);
    this->addEvaluatedField(outFields_[eq]);
  }

  this->setName(outPrefix+" Reorder AD Values");
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
panzer::ReorderADValues_Evaluator<EvalT, TRAITS>::
ReorderADValues_Evaluator(const std::string & outPrefix,
                          const std::vector<std::string> & inFieldNames,
                          const std::vector<std::string> & inDOFs,
                          const std::vector<std::string> & outDOFs,
                          const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                          const std::string & elementBlock,
                          const UniqueGlobalIndexerBase & indexerSrc,
                          const UniqueGlobalIndexerBase & indexerDest)
{ 
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());
  TEUCHOS_ASSERT(inDOFs.size()==outDOFs.size());

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
    outFields_.push_back(PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]));

    // tell the field manager that we depend on this field
    this->addDependentField(inFields_[eq]);
    this->addEvaluatedField(outFields_[eq]);
  }

  this->setName("Reorder AD Values");
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
void panzer::ReorderADValues_Evaluator<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData d, 
		      PHX::FieldManager<TRAITS>& fm)
{
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<inFields_.size();++fd) {
    // fill field data object
    this->utils.setFieldData(inFields_[fd],fm);
    this->utils.setFieldData(outFields_[fd],fm);
  }
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
void panzer::ReorderADValues_Evaluator<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  // just copy fields if there is no AD data
  //for(std::size_t i = 0; i < inFields_.size(); ++i)
  //  for(typename PHX::MDField<ScalarT>::size_type j = 0; j < inFields_[i].size(); ++j)
  //    outFields_[i][j] = inFields_[i][j];
  for(std::size_t i = 0; i < inFields_.size(); ++i)
   outFields_[i].deep_copy(inFields_[i]);
}

// **********************************************************************
// Jacobian
// **********************************************************************

template<typename TRAITS>
panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
ReorderADValues_Evaluator(const std::string & outPrefix,
                          const std::vector<std::string> & inFieldNames,
                          const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                          const std::string & elementBlock,
                          const UniqueGlobalIndexerBase & indexerSrc,
                          const UniqueGlobalIndexerBase & indexerDest)
{ 
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
    outFields_.push_back(PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]));

    // tell the field manager that we depend on this field
    this->addDependentField(inFields_[eq]);
    this->addEvaluatedField(outFields_[eq]);
  }

  buildSrcToDestMap(elementBlock,
                    indexerSrc,
                    indexerDest);

  this->setName(outPrefix+" Reorder AD Values");
}

// **********************************************************************

template<typename TRAITS>
panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
ReorderADValues_Evaluator(const std::string & outPrefix,
                          const std::vector<std::string> & inFieldNames,
                          const std::vector<std::string> & inDOFs,
                          const std::vector<std::string> & outDOFs,
                          const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                          const std::string & elementBlock,
                          const UniqueGlobalIndexerBase & indexerSrc,
                          const UniqueGlobalIndexerBase & indexerDest)
{ 
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());
  TEUCHOS_ASSERT(inDOFs.size()==outDOFs.size());

  // build the vector of fields that this is dependent on
  std::map<int,int> fieldNumberMaps;
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
    outFields_.push_back(PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]));

    // tell the field manager that we depend on this field
    this->addDependentField(inFields_[eq]);
    this->addEvaluatedField(outFields_[eq]);
  }

  // build a int-int map that associates fields
  for(std::size_t i=0;i<inDOFs.size();i++) {
    int srcFieldNum = indexerSrc.getFieldNum(inDOFs[i]);
    int dstFieldNum = indexerDest.getFieldNum(outDOFs[i]);
    TEUCHOS_ASSERT(srcFieldNum>=0);
    TEUCHOS_ASSERT(dstFieldNum>=0);

    fieldNumberMaps[srcFieldNum] = dstFieldNum;
  }

  buildSrcToDestMap(elementBlock,
                    fieldNumberMaps,
                    indexerSrc,
                    indexerDest);

  this->setName("Reorder AD Values");
}

// **********************************************************************
template<typename TRAITS>
void panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData d, 
		      PHX::FieldManager<TRAITS>& fm)
{
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<inFields_.size();++fd) {
    // fill field data object
    this->utils.setFieldData(inFields_[fd],fm);
    this->utils.setFieldData(outFields_[fd],fm);
  }
}

// **********************************************************************
template<typename TRAITS>
void panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  // for AD data do a reordering

  // TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"ERROR: panzer::ReorderADValues_Evaluator: This is currently broken for the Kokkkos Transition!  Contact Drekar team to fix!");

  for(std::size_t fieldIndex = 0; fieldIndex < inFields_.size(); ++fieldIndex) {

    const PHX::MDField<ScalarT>& inField = inFields_[fieldIndex];                                                                                                    
    const PHX::MDField<ScalarT>& outField = outFields_[fieldIndex];

    if(inField.size()>0) {

      switch (inFields_[fieldIndex].rank()) {
      case (1):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i) {
	  outField(i).val() = inField(i).val();
	  for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
	    outField(i).fastAccessDx(dx) = inField(i).fastAccessDx(dstFromSrcMap_[dx]);
	}
	break;
      case (2):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i)
	  for (typename PHX::MDField<ScalarT>::size_type j = 0; j < inField.dimension(1); ++j) {
	    outField(i,j).val() = inField(i,j).val();
	    for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
	      outField(i,j).fastAccessDx(dx) = inField(i,j).fastAccessDx(dstFromSrcMap_[dx]);
	  }
	break;
      case (3):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i)
	  for (typename PHX::MDField<ScalarT>::size_type j = 0; j < inField.dimension(1); ++j)
	    for (typename PHX::MDField<ScalarT>::size_type k = 0; k < inField.dimension(2); ++k) {
	      outField(i,j,k).val() = inField(i,j,k).val();
	      for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
		outField(i,j,k).fastAccessDx(dx) = inField(i,j,k).fastAccessDx(dstFromSrcMap_[dx]);
	    }
	break;
      case (4):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i)
	  for (typename PHX::MDField<ScalarT>::size_type j = 0; j < inField.dimension(1); ++j)
	    for (typename PHX::MDField<ScalarT>::size_type k = 0; k < inField.dimension(2); ++k)
	      for (typename PHX::MDField<ScalarT>::size_type l = 0; l < inField.dimension(3); ++l) {
		outField(i,j,k,l).val() = inField(i,j,k,l).val();
		for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
		  outField(i,j,k,l).fastAccessDx(dx) = inField(i,j,k,l).fastAccessDx(dstFromSrcMap_[dx]);
	      }
	break;
      case (5):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i)
	  for (typename PHX::MDField<ScalarT>::size_type j = 0; j < inField.dimension(1); ++j)
	    for (typename PHX::MDField<ScalarT>::size_type k = 0; k < inField.dimension(2); ++k)
	      for (typename PHX::MDField<ScalarT>::size_type l = 0; l < inField.dimension(3); ++l)
		for (typename PHX::MDField<ScalarT>::size_type m = 0; m < inField.dimension(4); ++m) {
		  outField(i,j,k,l,m).val() = inField(i,j,k,l,m).val();
		  for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
		    outField(i,j,k,l,m).fastAccessDx(dx) = inField(i,j,k,l,m).fastAccessDx(dstFromSrcMap_[dx]);
		}
	break;
      case (6):
	for (typename PHX::MDField<ScalarT>::size_type i = 0; i < inField.dimension(0); ++i)
	  for (typename PHX::MDField<ScalarT>::size_type j = 0; j < inField.dimension(1); ++j)
	    for (typename PHX::MDField<ScalarT>::size_type k = 0; k < inField.dimension(2); ++k)
	      for (typename PHX::MDField<ScalarT>::size_type l = 0; l < inField.dimension(3); ++l)
		for (typename PHX::MDField<ScalarT>::size_type m = 0; m < inField.dimension(4); ++m)
		  for (typename PHX::MDField<ScalarT>::size_type n = 0; n < inField.dimension(5); ++n) {
		    outField(i,j,k,l,m,n).val() = inField(i,j,k,l,m,n).val();
		    for (typename PHX::MDField<ScalarT>::size_type dx = 0; dx < Teuchos::as<typename PHX::MDField<ScalarT>::size_type>(dstFromSrcMap_.size()); ++dx)
		      outField(i,j,k,l,m,n).fastAccessDx(dx) = inField(i,j,k,l,m,n).fastAccessDx(dstFromSrcMap_[dx]);
		  }
	break;
      }

    }

  }

//Irina TOFIX
/*
  for(std::size_t i = 0; i < inFields_.size(); ++i) {

    for(typename PHX::MDField<ScalarT>::size_type j = 0; j < inFields_[i].size(); ++j) {
      // allocated scalar fields
      outFields_[i][j] = ScalarT(dstFromSrcMap_.size(), inFields_[i][j].val());

      ScalarT & outField = outFields_[i][j];
      const ScalarT & inField = inFields_[i][j];

      // the jacobian must be initialized, otherwise its just a value copy
      if(inField.size()>0) {
        // loop over the sensitivity indices: all DOFs on a cell
        outField.resize(dstFromSrcMap_.size());

        // copy jacobian entries correctly reordered
        for(std::size_t k=0;k<dstFromSrcMap_.size();k++) 
          outField.fastAccessDx(k) = inField.fastAccessDx(dstFromSrcMap_[k]);
      }
 
      outField.val() = inField.val();
    }
  }
*/
}

// **********************************************************************
template<typename TRAITS>
void panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
buildSrcToDestMap(const std::string & elementBlock,
                  const UniqueGlobalIndexerBase & indexerSrc,
                  const UniqueGlobalIndexerBase & indexerDest)
{
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);

  TEUCHOS_ASSERT(indexerSrc.getComm()!=Teuchos::null);
  TEUCHOS_ASSERT(indexerDest.getComm()!=Teuchos::null);

  const std::vector<int> & dstFieldsNum = indexerDest.getBlockFieldNumbers(elementBlock);

  // build a map between destination field numbers and source field numbers
  std::map<int,int> fieldNumberMaps;
  for(std::size_t i=0;i<dstFieldsNum.size();i++) {
    std::string fieldName = indexerDest.getFieldString(dstFieldsNum[i]);

    int srcFieldNum = indexerSrc.getFieldNum(fieldName);
    if(srcFieldNum>=0)
      fieldNumberMaps[srcFieldNum] = dstFieldsNum[i];
    else
      out << "Warning: Reorder AD Values can't find field \"" << fieldName << "\"" << std::endl;
  }

  buildSrcToDestMap(elementBlock,fieldNumberMaps,indexerSrc,indexerDest);
}

// **********************************************************************
template<typename TRAITS>
void panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
buildSrcToDestMap(const std::string & elementBlock,
                  const std::map<int,int> & fieldNumberMaps,
                  const UniqueGlobalIndexerBase & indexerSrc,
                  const UniqueGlobalIndexerBase & indexerDest)
{
  int maxDest = -1;
  std::map<int,int> offsetMap; // map from source to destination offsets
  for(std::map<int,int>::const_iterator itr=fieldNumberMaps.begin();
      itr!=fieldNumberMaps.end();++itr) {
    int srcField = itr->first;
    int dstField = itr->second;

    const std::vector<int> & srcOffsets = indexerSrc.getGIDFieldOffsets(elementBlock,srcField);
    const std::vector<int> & dstOffsets = indexerDest.getGIDFieldOffsets(elementBlock,dstField);

    // field should be the same size
    TEUCHOS_ASSERT(srcOffsets.size()==dstOffsets.size());
    for(std::size_t i=0;i<srcOffsets.size();i++) {
      offsetMap[srcOffsets[i]] = dstOffsets[i];

      // provides a size for allocating an array below: we will be able
      // to index into dstFromSrcMap_ in a simple way
      maxDest = dstOffsets[i]>maxDest ? dstOffsets[i] : maxDest;
    }
  }

  // Build map
  TEUCHOS_ASSERT(maxDest>0);
  dstFromSrcMap_ = std::vector<int>(maxDest+1,-1);
  for(std::map<int,int>::const_iterator itr=offsetMap.begin();
      itr!=offsetMap.end();++itr) {
    dstFromSrcMap_[itr->second] = itr->first;
  }
}

// **********************************************************************

#endif

