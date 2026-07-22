// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ReorderADValues_Evaluator_impl_hpp__
#define __Panzer_ReorderADValues_Evaluator_impl_hpp__


#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Panzer_GlobalIndexer.hpp"

#include "Phalanx_DataLayout.hpp"

template<typename EvalT,typename TRAITS>
panzer::ReorderADValues_Evaluator<EvalT, TRAITS>::
ReorderADValues_Evaluator(const std::string & outPrefix,
                          const std::vector<std::string> & inFieldNames,
                          const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                          const std::string & /* elementBlock */,
                          const GlobalIndexer & /* indexerSrc */,
                          const GlobalIndexer & /* indexerDest */)
{
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<const ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
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
                          const std::string & /* elementBlock */,
                          const GlobalIndexer & /* indexerSrc */,
                          const GlobalIndexer & /* indexerDest */)
{
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());
  TEUCHOS_ASSERT(inDOFs.size()==outDOFs.size());

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<const ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
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
evaluateFields(typename TRAITS::EvalData /* workset */)
{
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
                          const GlobalIndexer & indexerSrc,
                          const GlobalIndexer & indexerDest)
{
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());

  // build the vector of fields that this is dependent on
  inFields_.resize(inFieldNames.size());
  outFields_.resize(inFieldNames.size());
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_[eq] = PHX::MDField<const ScalarT>(inFieldNames[eq],fieldLayouts[eq]);
    outFields_[eq] = PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]);

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
                          const GlobalIndexer & indexerSrc,
                          const GlobalIndexer & indexerDest)
{
  TEUCHOS_ASSERT(inFieldNames.size()==fieldLayouts.size());
  TEUCHOS_ASSERT(inDOFs.size()==outDOFs.size());

  // build the vector of fields that this is dependent on
  std::map<int,int> fieldNumberMaps;
  for (std::size_t eq = 0; eq < inFieldNames.size(); ++eq) {
    inFields_.push_back(PHX::MDField<const ScalarT>(inFieldNames[eq],fieldLayouts[eq]));
    outFields_.push_back(PHX::MDField<ScalarT>(outPrefix+inFieldNames[eq],fieldLayouts[eq]));

    // tell the field manager that we depend on this field
    this->addDependentField(inFields_[eq]);
    this->addEvaluatedField(outFields_[eq]);
    // Don't share so we can avoid zeroing out off blck Jacobian entries
    this->addUnsharedField(outFields_[eq].fieldTag().clone());
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
evaluateFields(typename TRAITS::EvalData workset)
{
  // for AD data do a reordering
  for(std::size_t fieldIndex = 0; fieldIndex < inFields_.size(); ++fieldIndex) {

    const auto & inField_v = inFields_[fieldIndex].get_view();
    const auto & outField_v = outFields_[fieldIndex].get_view();
    const auto & dstFromSrcMap_v = dstFromSrcMapView_;

    if(inField_v.size()>0) {

      const auto rank = inField_v.rank();
      if (rank==1) {
        Kokkos::parallel_for("ReorderADValues: Jacobian rank 2",workset.num_cells,KOKKOS_LAMBDA(const int& i){
          outField_v(i).val() = inField_v(i).val();
          for (size_t dx = 0; dx < dstFromSrcMap_v.size(); ++dx)
            outField_v(i).fastAccessDx(dx) = inField_v(i).fastAccessDx(dstFromSrcMap_v(dx));
        });
      }
      else if (rank==2) {
        Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<2>> policy({0,0},{static_cast<int64_t>(workset.num_cells),static_cast<int64_t>(inField_v.extent(1))});
        Kokkos::parallel_for("ReorderADValues: Jacobian rank 2",policy,KOKKOS_LAMBDA(const int& i, const int& j){
          outField_v(i,j).val() = inField_v(i,j).val();
          for (size_t dx = 0; dx < dstFromSrcMap_v.size(); ++dx)
            outField_v(i,j).fastAccessDx(dx) = inField_v(i,j).fastAccessDx(dstFromSrcMap_v(dx));
        });
      }
      else if (rank==3) {
        Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{static_cast<int64_t>(workset.num_cells),
              static_cast<int64_t>(inField_v.extent(1)),static_cast<int64_t>(inField_v.extent(2))});
        Kokkos::parallel_for("ReorderADValues: Jacobian rank 2",policy,KOKKOS_LAMBDA(const int& i, const int& j, const int& k){
          outField_v(i,j,k).val() = inField_v(i,j,k).val();
          for (size_t dx = 0; dx < dstFromSrcMap_v.size(); ++dx)
            outField_v(i,j,k).fastAccessDx(dx) = inField_v(i,j,k).fastAccessDx(dstFromSrcMap_v(dx));
        });
      }
      else if (rank==4) {
        Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<4>> policy({0,0,0,0},{static_cast<int64_t>(workset.num_cells),
              static_cast<int64_t>(inField_v.extent(1)),static_cast<int64_t>(inField_v.extent(2)),
              static_cast<int64_t>(inField_v.extent(3))});
        Kokkos::parallel_for("ReorderADValues: Jacobian rank 2",policy,KOKKOS_LAMBDA(const int& i, const int& j, const int& k, const int& l){
          outField_v(i,j,k,l).val() = inField_v(i,j,k,l).val();
          for (size_t dx = 0; dx < dstFromSrcMap_v.size(); ++dx)
            outField_v(i,j,k,l).fastAccessDx(dx) = inField_v(i,j,k,l).fastAccessDx(dstFromSrcMap_v(dx));
        });
      }
      else if (rank==5) {
        Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<5>> policy({0,0,0,0,0},{static_cast<int64_t>(workset.num_cells),
              static_cast<int64_t>(inField_v.extent(1)),static_cast<int64_t>(inField_v.extent(2)),
              static_cast<int64_t>(inField_v.extent(3)),static_cast<int64_t>(inField_v.extent(4))});
        Kokkos::parallel_for("ReorderADValues: Jacobian rank 2",policy,KOKKOS_LAMBDA(const int& i, const int& j, const int& k, const int& l, const int& m){
          outField_v(i,j,k,l,m).val() = inField_v(i,j,k,l,m).val();
            for (size_t dx = 0; dx < dstFromSrcMap_v.size(); ++dx)
            outField_v(i,j,k,l,m).fastAccessDx(dx) = inField_v(i,j,k,l,m).fastAccessDx(dstFromSrcMap_v(dx));
        });
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR AD Reorder, rank size " << rank << " not supported!");
      }
    }
  }
}

// **********************************************************************
template<typename TRAITS>
void panzer::ReorderADValues_Evaluator<typename TRAITS::Jacobian, TRAITS>::
buildSrcToDestMap(const std::string & elementBlock,
                  const GlobalIndexer & indexerSrc,
                  const GlobalIndexer & indexerDest)
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
                  const GlobalIndexer & indexerSrc,
                  const GlobalIndexer & indexerDest)
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
  std::vector<int> dstFromSrcMap(maxDest+1,-1);
  for(std::map<int,int>::const_iterator itr=offsetMap.begin();
      itr!=offsetMap.end();++itr) {
    dstFromSrcMap[itr->second] = itr->first;
  }

  dstFromSrcMapView_ = Kokkos::View<int*>("dstFromSrcMapView_",dstFromSrcMap.size());
  auto dfsm_host = Kokkos::create_mirror_view(Kokkos::HostSpace(),dstFromSrcMapView_);
  for (size_t i=0; i < dstFromSrcMapView_.size(); ++i)
    dfsm_host(i) = dstFromSrcMap[i];

  Kokkos::deep_copy(dstFromSrcMapView_,dfsm_host);
}

// **********************************************************************

#endif
