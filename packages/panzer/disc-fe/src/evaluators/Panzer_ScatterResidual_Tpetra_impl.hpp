// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_RESIDUAL_TPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_TPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

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
ScatterResidual_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
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

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  fieldIds_.resize(scatterFields_.size());
  const Workset & workset_0 = (*d.worksets_)[0];
  std::string blockId = this->wda(workset_0).block_id;


  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldIds_[fd]);
    scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
    Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
  }
  scratch_lids_ = PHX::View<LO**>("lids",scatterFields_[0].extent(0),
                                                 globalIndexer_->getElementBlockGIDCount(blockId));

}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
  }
}


// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
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
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // this is the list of parameters and their names that this scatter has to account for
  std::vector<std::string> activeParameters =
    rcp_dynamic_cast<ParameterList_GlobalEvaluationData>(d.gedc->getDataObject("PARAMETER_NAMES"))->getActiveParameters();

  dfdp_vectors_.clear();
  for(std::size_t i=0;i<activeParameters.size();i++) {
    RCP<typename LOC::VectorType> vec =
      rcp_dynamic_cast<LOC>(d.gedc->getDataObject(activeParameters[i]),true)->get_f();
    Teuchos::ArrayRCP<double> vec_array = vec->get1dViewNonConst();
    dfdp_vectors_.push_back(vec_array);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // loop over each field to be scattered
      for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            ScalarT value = (scatterFields_[fieldIndex])(worksetCellIndex,basis);
            for(int d=0;d<value.size();d++)
              dfdp_vectors_[d][lid] += value.fastAccessDx(d);
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                       const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
   , my_derivative_size_(0)
   , other_derivative_size_(0)
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
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  fieldIds_.resize(scatterFields_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  std::string blockId = this->wda(workset_0).block_id;

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    int fieldNum = fieldIds_[fd];
    const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
    scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
    Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
  }

  my_derivative_size_ = globalIndexer_->getElementBlockGIDCount(blockId);
  if (Teuchos::nonnull(workset_0.other)) {
    auto otherBlockId = workset_0.other->block_id;
    other_derivative_size_ = globalIndexer_->getElementBlockGIDCount(otherBlockId);
  }
  scratch_lids_ = Kokkos::View<LO**, Kokkos::LayoutRight, PHX::Device>(
    "lids", scatterFields_[0].extent(0), my_derivative_size_ + other_derivative_size_ );
  scratch_vals_ = Kokkos::View<typename Sacado::ScalarType<ScalarT>::type**, Kokkos::LayoutRight, PHX::Device>(
    "vals", scatterFields_[0].extent(0), my_derivative_size_ + other_derivative_size_ );
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
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

  Kokkos::View<const LO**, Kokkos::LayoutRight, PHX::Device> lids; // local indices for unknowns.
  Kokkos::View<typename Sacado::ScalarType<ScalarT>::type**, Kokkos::LayoutRight, PHX::Device> vals;
  PHX::View<const int*> offsets; // how to get a particular field
  FieldType field;


  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    int numIds = lids.extent(1);

    // loop over the basis functions (currently they are nodes)
    for(std::size_t basis=0; basis < offsets.extent(0); basis++) {
       typename FieldType::array_type::reference_type scatterField = field(cell,basis);
       int offset = offsets(basis);
       LO lid    = lids(cell,offset);

       // Sum residual
       if(fillResidual)
         Kokkos::atomic_add(&r_data(lid,0), scatterField.val());

       // loop over the sensitivity indices: all DOFs on a cell
       for(int sensIndex=0;sensIndex<numIds;++sensIndex)
          vals(cell,sensIndex) = scatterField.fastAccessDx(sensIndex);

       // Sum Jacobian
       jac.sumIntoValues(lid, &lids(cell,0), numIds, &vals(cell,0), true, true);
    } // end basis
  }
};

template <typename ScalarT,typename LO,typename GO,typename NodeT>
class ScatterResidual_Residual_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;

  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;

  PHX::View<const LO**> lids;    // local indices for unknowns
  PHX::View<const int*> offsets; // how to get a particular field
  FieldType field;


  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {

    // loop over the basis functions (currently they are nodes)
    for(std::size_t basis=0; basis < offsets.extent(0); basis++) {
       int offset = offsets(basis);
       LO lid    = lids(cell,offset);
       Kokkos::atomic_add(&r_data(lid,0), field(cell,basis));

   } // end basis
  }
};

}
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // for convenience pull out some objects from workset
  std::string blockId = this->wda(workset).block_id;

  Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();

  globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);

  ScatterResidual_Residual_Functor<ScalarT,LO,GO,NodeT> functor;
  functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
  functor.lids = scratch_lids_;

  // for each field, do a parallel for loop
  for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
    functor.offsets = scratch_offsets_[fieldIndex];
    functor.field = scatterFields_[fieldIndex];

    Kokkos::parallel_for(workset.num_cells,functor);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   typedef typename LOC::CrsMatrixType::local_matrix_device_type LocalMatrixT;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f();
   Teuchos::RCP<typename LOC::CrsMatrixType> Jac = tpetraContainer_->get_A();

   // Cache scratch lids. For interface bc problems the derivative
   // dimension extent spans two cells. Use subviews to get the self
   // lids and the other lids.
   if (Teuchos::nonnull(workset.other)) {
     auto my_scratch_lids = Kokkos::subview(scratch_lids_,Kokkos::ALL,std::make_pair(0,my_derivative_size_));
     globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,my_scratch_lids);
     auto other_scratch_lids = Kokkos::subview(scratch_lids_,Kokkos::ALL,std::make_pair(my_derivative_size_,my_derivative_size_ + other_derivative_size_));
     globalIndexer_->getElementLIDs(workset.other->cell_local_ids_k,other_scratch_lids);
   }
   else {
     globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);
   }

   ScatterResidual_Jacobian_Functor<ScalarT,LO,GO,NodeT,LocalMatrixT> functor;
   functor.fillResidual = (r!=Teuchos::null);
   if(functor.fillResidual)
     functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
   functor.jac = Jac->getLocalMatrixDevice();
   functor.lids = scratch_lids_;
   functor.vals = scratch_vals_;

   // for each field, do a parallel for loop
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     functor.offsets = scratch_offsets_[fieldIndex];
     functor.field = scatterFields_[fieldIndex];

     Kokkos::parallel_for(workset.num_cells,functor);
   }

}

// **********************************************************************

#endif
