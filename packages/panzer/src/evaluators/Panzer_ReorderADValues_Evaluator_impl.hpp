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

template<typename EvalT,typename Traits>
panzer::ReorderADValues_Evaluator<EvalT, Traits>::
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
template<typename EvalT,typename Traits>
void panzer::ReorderADValues_Evaluator<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<inFields_.size();++fd) {
    // fill field data object
    this->utils.setFieldData(inFields_[fd],fm);
    this->utils.setFieldData(outFields_[fd],fm);
  }
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ReorderADValues_Evaluator<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  // just copy fields if there is no AD data
  for(std::size_t i = 0; i < inFields_.size(); ++i)
    for(std::size_t j = 0; j < inFields_[i].size(); ++j)
      outFields_[i][j] = inFields_[i][j];
}

// **********************************************************************
// Jacobian
// **********************************************************************

template<typename Traits>
panzer::ReorderADValues_Evaluator<panzer::Traits::Jacobian, Traits>::
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
template<typename Traits>
void panzer::ReorderADValues_Evaluator<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<inFields_.size();++fd) {
    // fill field data object
    this->utils.setFieldData(inFields_[fd],fm);
    this->utils.setFieldData(outFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer::ReorderADValues_Evaluator<panzer::Traits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  // for AD data do a reordering
  for(std::size_t i = 0; i < inFields_.size(); ++i) {

    for(std::size_t j = 0; j < inFields_[i].size(); ++j) {
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
}

// **********************************************************************
template<typename Traits>
void panzer::ReorderADValues_Evaluator<panzer::Traits::Jacobian, Traits>::
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

