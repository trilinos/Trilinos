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

#ifndef PANZER_STK_SCATTER_CELL_AVG_VECTOR_IMPL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_VECTOR_IMPL_HPP

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
ScatterCellAvgVector<EvalT, Traits>::
ScatterCellAvgVector(
  const Teuchos::ParameterList& p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

  std::string scatterName = p.get<std::string>("Scatter Name");
 
  const std::vector<std::string> & names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::IntegrationRule> intRule = 
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // build dependent fields
  scatterFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) 
  {
    scatterFields_[fd] = PHX::MDField<const ScalarT,Cell,Point,Dim>(names[fd],intRule->dl_vector);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Cell Vectors");
}


template<typename EvalT, typename Traits>
void
ScatterCellAvgVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* d */,
  PHX::FieldManager<Traits>&  fm)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) 
  {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(stk::topology::ELEMENT_RANK, fieldName);

    // setup the field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}


template<typename EvalT, typename Traits>
void
ScatterCellAvgVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  panzer::MDFieldArrayFactory af("",true);

  // for convenience pull out some objects from workset
  const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;
  std::string blockId = this->wda(workset).block_id;
  std::string d_mod[3] = {"X","Y","Z"};
   
  // loop over the number of vector fields requested for exodus output
  for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) 
  {
    PHX::MDField<const ScalarT,panzer::Cell,panzer::Point,panzer::Dim> & field = scatterFields_[fieldIndex];
    std::string fieldName = field.fieldTag().name();
    int numCells = field.extent(0);
    int numPoints = field.extent(1);  
    int numDims = field.extent(2);
    
    for (int dim = 0; dim < numDims; dim++)
    {  
      // std::vector<double> average(numCells,0.0);
      PHX::MDField<double,panzer::Cell,panzer::NODE> average = af.buildStaticArray<double,panzer::Cell,panzer::NODE>("",numCells,1);
      
      // write to double field
      for(int i = 0; i < numCells; i++)  // loop over cells
      {
         average(i,0) = 0.0;
         for(int j = 0; j < numPoints; j++)  // loop over IPs 
            average(i,0) += Sacado::ScalarValue<ScalarT>::eval(field(i,j,dim));
         
         average(i,0) /= numPoints;
      }

      mesh_->setCellFieldData(fieldName+d_mod[dim],blockId,localCellIds,average);
      
    }
  }
  
}

} // end panzer_stk

#endif
