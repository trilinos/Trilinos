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

#ifndef PANZER_SCATTER_RESIDUAL_TPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_TPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                       const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , globalDataKey_("Residual Scatter Container")
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();
   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // loop over each field to be scattered
      for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            r_array[lid] += (scatterFields_[fieldIndex])(worksetCellIndex,basis);
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                       const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , globalDataKey_("Residual Scatter Container")
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   TEUCHOS_ASSERT(false);

   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();
   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // loop over each field to be scattered
      for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            r_array[lid] += (scatterFields_[fieldIndex])(worksetCellIndex,basis).val();
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
                       const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  scratch_offsets_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  fieldIds_.resize(scatterFields_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  std::string blockId = workset_0.block_id;

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    int fieldNum = fieldIds_[fd];
    const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
    scratch_offsets_[fd] = Kokkos::View<int*,PHX::Device>("offsets",offsets.size());
    for(std::size_t i=0;i<offsets.size();i++)
      scratch_offsets_[fd](i) = offsets[i];
  }

  scratch_lids_ = Kokkos::View<LO**,PHX::Device>("lids",scatterFields_[0].dimension_0(),
                                                 globalIndexer_->getElementBlockGIDCount(blockId));
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
  }
}


// **********************************************************************
namespace panzer {
namespace {

template <typename ScalarT,typename LO,typename GO,typename NodeT,typename LocalMatrixT>
class ScatterResidual_Jacobian_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;

  bool fillResidual;
  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;
  LocalMatrixT jac; // Kokkos jacobian type

  Kokkos::View<const LO**,PHX::Device> lids;    // local indices for unknowns
  Kokkos::View<const int*,PHX::Device> offsets; // how to get a particular field
  FieldType field;


  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    LO cLIDs[256];
    typename Sacado::ScalarType<ScalarT>::type vals[256];
    int numIds = lids.dimension_1();

    for(int i=0;i<numIds;i++)
      cLIDs[i] = lids(cell,i);

    // loop over the basis functions (currently they are nodes)
    for(std::size_t basis=0; basis < offsets.dimension_0(); basis++) {
       typename FieldType::array_type::reference_type scatterField = field(cell,basis);
       int offset = offsets(basis);
       LO lid    = lids(cell,offset);

       // Sum residual
       if(fillResidual)
         Kokkos::atomic_add(&r_data(lid,0), scatterField.val());

       // loop over the sensitivity indices: all DOFs on a cell
       for(int sensIndex=0;sensIndex<numIds;++sensIndex)
          vals[sensIndex] = scatterField.fastAccessDx(sensIndex);

       // Sum Jacobian
       jac.sumIntoValues(lid, cLIDs,numIds, vals, true);
    } // end basis
  }
};
}
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

#if 0
   std::vector<GO> GIDs;
   std::vector<LO> cLIDs, rLIDs;
   std::vector<double> jacRow;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();
   Teuchos::RCP<typename LOC::CrsMatrixType> Jac = tpetraContainer_->get_A();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets" may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];


      rLIDs = globalIndexer_->getElementLIDs(cellLocalId);
      cLIDs = rLIDs;

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over the basis functions (currently they are nodes)
         for(std::size_t rowBasisNum = 0; rowBasisNum < elmtOffset.size(); rowBasisNum++) {
            const ScalarT scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,rowBasisNum);
            int rowOffset = elmtOffset[rowBasisNum];
            LO row = rLIDs[rowOffset];

            // Sum residual
            if(r!=Teuchos::null)
               r->sumIntoLocalValue(row,scatterField.val());

            // loop over the sensitivity indices: all DOFs on a cell
            jacRow.resize(scatterField.size());

            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);

            // Sum Jacobian
            Jac->sumIntoLocalValues(row, cLIDs, jacRow);

         } // end rowBasisNum
      } // end fieldIndex
   }
#else
   typedef typename LOC::CrsMatrixType::local_matrix_type LocalMatrixT;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();
   Teuchos::RCP<typename LOC::CrsMatrixType> Jac = tpetraContainer_->get_A();

   globalIndexer_->getElementLIDs(workset.cell_local_ids_k,scratch_lids_);

   ScatterResidual_Jacobian_Functor<ScalarT,LO,GO,NodeT,LocalMatrixT> functor;
   functor.fillResidual = (r!=Teuchos::null);
   if(functor.fillResidual)
     functor.r_data = r->template getLocalView<PHX::Device>();
   functor.jac = Jac->getLocalMatrix();
   functor.lids = scratch_lids_;

   // for each field, do a parallel for loop
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     functor.offsets = scratch_offsets_[fieldIndex];
     functor.field = scatterFields_[fieldIndex];

     Kokkos::parallel_for(workset.num_cells,functor);
   }
#endif
}

// **********************************************************************

#endif
