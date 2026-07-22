// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SCATTER_CELL_AVG_QUANTITY_IMPL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_QUANTITY_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace panzer_stk {

template<typename EvalT, typename Traits>
ScatterCellAvgQuantity<EvalT, Traits>::
ScatterCellAvgQuantity(
  const Teuchos::ParameterList& p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::Point;

  std::string scatterName = p.get<std::string>("Scatter Name");

  // if it's there pull the ouptut scaling
  if (p.isParameter("Variable Scale Factors Map"))
    varScaleFactors_ = p.get<Teuchos::RCP<std::map<std::string,double>>>("Variable Scale Factors Map");

  const std::vector<std::string> & names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::IntegrationRule> intRule =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // build dependent fields
  scatterFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] =
      PHX::MDField<const ScalarT,Cell,Point>(names[fd],intRule->dl_scalar);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Cell Fields");
}

template<typename EvalT, typename Traits>
void
ScatterCellAvgQuantity<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* d */,
  PHX::FieldManager<Traits>& /* fm */)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<double>(stk::topology::ELEMENT_RANK, fieldName);
  }
}

template<typename EvalT, typename Traits>
void
ScatterCellAvgQuantity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
   panzer::MDFieldArrayFactory af("",true);

   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;
   std::string blockId = this->wda(workset).block_id;

   for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
     auto field = scatterFields_[fieldIndex].get_static_view();
     auto average = af.buildStaticArray<double,panzer::Cell,panzer::NODE>("",field.extent(0),1).get_static_view();
     // write to double field
     Kokkos::parallel_for("ScatterCellAvgQuantity",field.extent(0), KOKKOS_LAMBDA(int i) {
       for(unsigned j=0; j<field.extent(1);j++)
	 average(i,0) += Sacado::scalarValue(field(i,j));
       average(i,0) /= field.extent(1);
     });
     Kokkos::fence();

     std::string varname = scatterFields_[fieldIndex].fieldTag().name();
     double scalef = 1.0;

     if (!varScaleFactors_.is_null())
     {
       std::map<std::string,double> *tmp_sfs = varScaleFactors_.get();
       if(tmp_sfs->find(varname) != tmp_sfs->end())
         scalef = (*tmp_sfs)[varname];
     }

     mesh_->setCellFieldData(varname,blockId,localCellIds,average,scalef);
   }
}

} // end panzer_stk

#endif
