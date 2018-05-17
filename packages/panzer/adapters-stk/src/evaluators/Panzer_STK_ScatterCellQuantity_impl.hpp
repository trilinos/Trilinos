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
  typename Traits::SetupData  /* d */,
  PHX::FieldManager<Traits>&  fm)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    // setup the field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
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
         value(i,0) = Sacado::ScalarValue<ScalarT>::eval(field(i));

      mesh_->setCellFieldData(field.fieldTag().name(),blockId,localCellIds,value);
   }
}

} // end panzer_stk

#endif
