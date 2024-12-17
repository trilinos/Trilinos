// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_TPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_TPETRA_IMPL_HPP

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

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
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

  // determine if we are scattering an initial condition
  scatterIC_ = p.isParameter("Scatter Initial Condition") ? p.get<bool>("Scatter Initial Condition") : false;

  Teuchos::RCP<PHX::DataLayout> dl = (!scatterIC_) ? 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional :
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional ;
  if (!scatterIC_) {
    side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
    local_side_id_ = p.get<int>("Local Side ID");
    scratch_basisIds_.resize(names.size());
  }
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  scratch_offsets_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.isParameter("Check Apply BC") ? p.get<bool>("Check Apply BC") : false;
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
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

    if (!scatterIC_) {
      const std::pair<std::vector<int>,std::vector<int> > & indicePair 
        = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldIds_[fd], side_subcell_dim_, local_side_id_);
      const std::vector<int> & offsets = indicePair.first;
      const std::vector<int> & basisIdMap = indicePair.second;

      scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
      Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
 
      scratch_basisIds_[fd] = PHX::View<int*>("basisIds",basisIdMap.size()); 
      Kokkos::deep_copy(scratch_basisIds_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(basisIdMap.data(), basisIdMap.size()));
 
    } else {
      const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldIds_[fd]);
      scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
      Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
    }
 }

  scratch_lids_ = PHX::View<LO**>("lids",scatterFields_[0].extent(0),
                                  globalIndexer_->getElementBlockGIDCount(blockId));
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
namespace panzer {
namespace {

template <typename ScalarT,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Residual_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> ScalarFieldType;
  typedef PHX::MDField<const bool,Cell,NODE> BoolFieldType;

  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;
  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> dirichlet_counter;

  PHX::View<const LO**> lids;    // local indices for unknowns
  PHX::View<const int*> offsets; // how to get a particular field
  PHX::View<const int*> basisIds;
  ScalarFieldType field;
  BoolFieldType applyBC;

  bool checkApplyBC;

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {

    // loop over the basis functions (currently they are nodes)
    for(std::size_t basis=0; basis < offsets.extent(0); basis++) {
       int offset = offsets(basis);
       LO lid    = lids(cell,offset);
       if (lid<0) continue; // not on this processor

       int basisId = basisIds(basis);
       if (checkApplyBC)
         if(!applyBC(cell,basisId)) continue;

       r_data(lid,0) = field(cell,basisId);

       // record that you set a dirichlet condition
       dirichlet_counter(lid,0) = 1.0;

   } // end basis
  }
};

template <typename ScalarT,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidualIC_Residual_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;

  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;
  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> dirichlet_counter;

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
       if (lid<0) continue; // not on this processor

       r_data(lid,0) = field(cell,basis);

       // record that you set a dirichlet condition
       dirichlet_counter(lid,0) = 1.0;

   } // end basis
  }
};
}
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   std::vector<GO> GIDs;
   std::vector<LO> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;

  globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);

  Teuchos::RCP<typename LOC::VectorType> r = (!scatterIC_) ? 
    tpetraContainer_->get_f() :
    tpetraContainer_->get_x(); 

  if (scatterIC_) {
    ScatterDirichletResidualIC_Residual_Functor<ScalarT,LO,GO,NodeT> functor;
    functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.lids = scratch_lids_;
    functor.dirichlet_counter = dirichletCounter_->getLocalViewDevice(Tpetra::Access::ReadWrite);

      // for each field, do a parallel for loop
    for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      functor.offsets = scratch_offsets_[fieldIndex];
      functor.field = scatterFields_[fieldIndex];

      Kokkos::parallel_for(workset.num_cells,functor);
    }
  } else {
    ScatterDirichletResidual_Residual_Functor<ScalarT,LO,GO,NodeT> functor;
    functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.lids = scratch_lids_;
    functor.dirichlet_counter = dirichletCounter_->getLocalViewDevice(Tpetra::Access::ReadWrite);

      // for each field, do a parallel for loop
    for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      functor.offsets = scratch_offsets_[fieldIndex];
      functor.field = scatterFields_[fieldIndex];
      if (checkApplyBC_) functor.applyBC = applyBC_[fieldIndex];
      functor.checkApplyBC = checkApplyBC_;
      functor.basisIds = scratch_basisIds_[fieldIndex];

      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }

}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************


template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
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

  // determine if we are scattering an initial condition
  scatterIC_ = p.isParameter("Scatter Initial Condition") ? p.get<bool>("Scatter Initial Condition") : false;

  Teuchos::RCP<PHX::DataLayout> dl = (!scatterIC_) ? 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional :
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional ;
  if (!scatterIC_) {
    side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
    local_side_id_ = p.get<int>("Local Side ID");
    scratch_basisIds_.resize(names.size());
  }
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  scratch_offsets_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.isParameter("Check Apply BC") ? p.get<bool>("Check Apply BC") : false;
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
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

    if (!scatterIC_) {
      const std::pair<std::vector<int>,std::vector<int> > & indicePair 
        = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldIds_[fd], side_subcell_dim_, local_side_id_);
      const std::vector<int> & offsets = indicePair.first;
      const std::vector<int> & basisIdMap = indicePair.second;

      scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
      Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
 
      scratch_basisIds_[fd] = PHX::View<int*>("basisIds",basisIdMap.size()); 
      Kokkos::deep_copy(scratch_basisIds_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(basisIdMap.data(), basisIdMap.size()));
 
    } else {
      const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldIds_[fd]);
      scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
      Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
    }
 }

  scratch_lids_ = PHX::View<LO**>("lids",scatterFields_[0].extent(0),
                                  globalIndexer_->getElementBlockGIDCount(blockId));

}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // this is the list of parameters and their names that this scatter has to account for
  std::vector<std::string> activeParameters =
    rcp_dynamic_cast<ParameterList_GlobalEvaluationData>(d.gedc->getDataObject("PARAMETER_NAMES"))->getActiveParameters();

  dfdpFieldsVoV_.initialize("ScatterResidual_Tpetra<Tangent>::dfdpFieldsVoV_",activeParameters.size());

  for(std::size_t i=0;i<activeParameters.size();i++) {
    RCP<typename LOC::VectorType> vec =
      rcp_dynamic_cast<LOC>(d.gedc->getDataObject(activeParameters[i]),true)->get_f();
    auto dfdp_view = vec->getLocalViewDevice(Tpetra::Access::ReadWrite);

    dfdpFieldsVoV_.addView(dfdp_view,i);
  }

  dfdpFieldsVoV_.syncHostToDevice();

}

// **********************************************************************
namespace panzer {
namespace {

template <typename ScalarT,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Tangent_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> ScalarFieldType;
  typedef PHX::MDField<const bool,Cell,NODE> BoolFieldType;

  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;
  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> dirichlet_counter;

  Kokkos::View<Kokkos::View<double**,Kokkos::LayoutLeft,PHX::Device>*>  dfdp_fields; // tangent fields
  double num_params;

  PHX::View<const LO**> lids;    // local indices for unknowns
  PHX::View<const int*> offsets; // how to get a particular field
  PHX::View<const int*> basisIds;
  ScalarFieldType field;
  BoolFieldType applyBC;

  bool checkApplyBC;

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {

    // loop over the basis functions (currently they are nodes)
    for(std::size_t basis=0; basis < offsets.extent(0); basis++) {
       int offset = offsets(basis);
       LO lid    = lids(cell,offset);
       if (lid<0) continue; // not on this processor

       int basisId = basisIds(basis);
       if (checkApplyBC)
         if(!applyBC(cell,basisId)) continue;

       r_data(lid,0) = field(cell,basisId).val();

       // loop over the tangents
       for(int i_param=0; i_param<num_params; i_param++)
         dfdp_fields(i_param)(lid,0) = field(cell,basisId).fastAccessDx(i_param);

       // record that you set a dirichlet condition
       dirichlet_counter(lid,0) = 1.0;

   } // end basis
  }
};

template <typename ScalarT,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidualIC_Tangent_Functor {
public:
  typedef typename PHX::Device execution_space;
  typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;

  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> r_data;
  Kokkos::View<double**, Kokkos::LayoutLeft,PHX::Device> dirichlet_counter;

  Kokkos::View<Kokkos::View<double**,Kokkos::LayoutLeft,PHX::Device>*>  dfdp_fields; // tangent fields
  double num_params;

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
       if (lid<0) continue; // not on this processor

       r_data(lid,0) = field(cell,basis).val();

       // loop over the tangents
       for(int i_param=0; i_param<num_params; i_param++)
          dfdp_fields(i_param)(lid,0) = field(cell,basis).fastAccessDx(i_param);

       // record that you set a dirichlet condition
       dirichlet_counter(lid,0) = 1.0;

   } // end basis
  }
};
}
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   std::vector<GO> GIDs;
   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
     
  globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);

   Teuchos::RCP<typename LOC::VectorType> r = (!scatterIC_) ?
     tpetraContainer_->get_f() :
     tpetraContainer_->get_x();

  if (scatterIC_) {
    ScatterDirichletResidualIC_Tangent_Functor<ScalarT,LO,GO,NodeT> functor;
    functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.lids = scratch_lids_;
    functor.dirichlet_counter = dirichletCounter_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.dfdp_fields = dfdpFieldsVoV_.getViewDevice();

      // for each field, do a parallel for loop
    for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      functor.offsets = scratch_offsets_[fieldIndex];
      functor.field = scatterFields_[fieldIndex];
      functor.num_params = Kokkos::dimension_scalar(scatterFields_[fieldIndex].get_view())-1;

      Kokkos::parallel_for(workset.num_cells,functor);
    }
  } else {
    ScatterDirichletResidual_Tangent_Functor<ScalarT,LO,GO,NodeT> functor;
    functor.r_data = r->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.lids = scratch_lids_;
    functor.dirichlet_counter = dirichletCounter_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    functor.dfdp_fields = dfdpFieldsVoV_.getViewDevice();

      // for each field, do a parallel for loop
    for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      functor.offsets = scratch_offsets_[fieldIndex];
      functor.field = scatterFields_[fieldIndex];
      if (checkApplyBC_) functor.applyBC = applyBC_[fieldIndex];
      functor.checkApplyBC = checkApplyBC_;
      functor.basisIds = scratch_basisIds_[fieldIndex];
      functor.num_params = Kokkos::dimension_scalar(scatterFields_[fieldIndex].get_view())-1;

      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }

}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
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
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
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

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].extent(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   std::vector<GO> GIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f(); 
   Teuchos::RCP<typename LOC::CrsMatrixType> Jac = tpetraContainer_->get_A();

   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();
   Teuchos::ArrayRCP<double> dc_array = dirichletCounter_->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   auto LIDs = globalIndexer_->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);
   // loop over each field to be scattered
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     auto scatterFields_h = Kokkos::create_mirror_view(scatterFields_[fieldIndex].get_static_view());
     Kokkos::deep_copy(scatterFields_h, scatterFields_[fieldIndex].get_static_view());
     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

       globalIndexer_->getElementGIDs(cellLocalId,GIDs); 
   
       // this call "should" get the right ordering according to the Intrepid2 basis
       const std::pair<std::vector<int>,std::vector<int> > & indicePair 
	 = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
       const std::vector<int> & elmtOffset = indicePair.first;
       const std::vector<int> & basisIdMap = indicePair.second;
       
       // loop over basis functions
       for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
	 int offset = elmtOffset[basis];
	 int lid = LIDs_h(cellLocalId, offset);
	 if(lid<0) // not on this processor
	   continue;

	 int basisId = basisIdMap[basis];

	 if (checkApplyBC_)
	   if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
	     continue;

	 // zero out matrix row
	 {
               std::size_t sz = Jac->getNumEntriesInLocalRow(lid);
               std::size_t numEntries = 0;
	       typename LOC::CrsMatrixType::nonconst_local_inds_host_view_type rowIndices("indices", sz);
	       typename LOC::CrsMatrixType::nonconst_values_host_view_type rowValues("values", sz);

               // Jac->getLocalRowView(lid,numEntries,rowValues,rowIndices);
               Jac->getLocalRowCopy(lid,rowIndices,rowValues,numEntries);

               for(std::size_t i=0;i<numEntries;i++)
		 rowValues(i) = 0.0;

               Jac->replaceLocalValues(lid,rowIndices,rowValues);
            }
 
            GO gid = GIDs[offset];
            const ScalarT scatterField = scatterFields_h(worksetCellIndex,basisId);
    
            r_array[lid] = scatterField.val();
            dc_array[lid] = 1.0; // mark row as dirichlet
    
            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);
            TEUCHOS_ASSERT(jacRow.size()==GIDs.size());
    
            Jac->replaceGlobalValues(gid, GIDs, jacRow);
       }
     }
   }
}

// **********************************************************************

#endif
