// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SCATTER_CELL_QUANTITY_IMPL_HPP
#define PANZER_STK_SCATTER_CELL_QUANTITY_IMPL_HPP

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
ScatterCellQuantity<EvalT, Traits>::
ScatterCellQuantity(
  const Teuchos::ParameterList& p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::Point;

  std::string scatterName = p.get<std::string>("Scatter Name");
  int worksetSize = p.get<int>("Workset Size");

  // if it's there pull the ouptut scaling
  if (p.isParameter("Variable Scale Factors Map"))
    varScaleFactors_ = p.get<Teuchos::RCP<std::map<std::string,double>>>("Variable Scale Factors Map");

  const std::vector<std::string> & names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));

  // build dependent fields
  scatterFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] = PHX::MDField<const ScalarT,Cell>(names[fd],dl_cell);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Cell Quantity Fields");
}

template<typename EvalT, typename Traits>
void
ScatterCellQuantity<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* d */,
  PHX::FieldManager<Traits>& /* fm */)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd)
    std::string fieldName = scatterFields_[fd].fieldTag().name();
}

template<typename EvalT, typename Traits>
void
ScatterCellQuantity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
   panzer::MDFieldArrayFactory af("",true);

   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;
   std::string blockId = this->wda(workset).block_id;

   for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
      PHX::MDField<const ScalarT,panzer::Cell> & field = scatterFields_[fieldIndex];
      // std::vector<double> value(field.extent(0),0.0);
      PHX::MDField<double,panzer::Cell,panzer::NODE> value = af.buildStaticArray<double,panzer::Cell,panzer::NODE>("",field.extent(0),1);

      // write to double field
      for(unsigned i=0; i<field.extent(0);i++)
         value(i,0) = Sacado::scalarValue(field(i));

      std::string varname = field.fieldTag().name();
      double scalef = 1.0;

      if (!varScaleFactors_.is_null())
      {
        std::map<std::string,double> *tmp_sfs = varScaleFactors_.get();
        if(tmp_sfs->find(varname) != tmp_sfs->end())
          scalef = (*tmp_sfs)[varname];
      }

      mesh_->setCellFieldData(field.fieldTag().name(),blockId,localCellIds,value.get_view(),scalef);
   }
}

} // end panzer_stk

#endif
