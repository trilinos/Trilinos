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

#ifndef PANZER_STK_SCATTER_VECTOR_FIELDS_IMPL_HPP
#define PANZER_STK_SCATTER_VECTOR_FIELDS_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_Traits.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

PHX_EVALUATOR_CTOR(ScatterVectorFields,p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
   TEUCHOS_ASSERT(false);
}

template <typename EvalT,typename TraitsT>
ScatterVectorFields<EvalT,TraitsT>::
ScatterVectorFields(const std::string & scatterName,
              const Teuchos::RCP<STK_Interface> mesh,
              const Teuchos::RCP<const panzer::PointRule> & pointRule,
              const std::vector<std::string> & names)
   : mesh_(mesh)
{
  using panzer::Cell;
  using panzer::IP;
  using panzer::Dim;

  spatialDimension_ = pointRule->spatial_dimension;

  pointField_ = PHX::MDField<ScalarT,Cell,IP,Dim>(pointRule->getName()+"_point_coords",pointRule->dl_vector);

  // build dependent fields
  names_ = names;
  scatterFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] = 
      PHX::MDField<ScalarT,Cell,IP,Dim>(names_[fd]+"_"+pointRule->getName(),pointRule->dl_vector);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Vector Fields");
}

PHX_POST_REGISTRATION_SETUP(ScatterVectorFields,d,fm)
{
  this->utils.setFieldData(pointField_,fm);
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) {
    // setup the field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

PHX_EVALUATE_FIELDS(ScatterVectorFields,workset)
{
   TEUCHOS_ASSERT(false);
}

template < >
void ScatterVectorFields<panzer::Traits::Residual,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData workset)
{
   std::vector<std::string> dimStrings(3);
   dimStrings[0] = "X";
   dimStrings[1] = "Y";
   dimStrings[2] = "Z";

   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   std::string blockId = workset.block_id;

/*
   std::cout << "Point Field " << std::endl;
   for(int i=0; i<pointField_.dimension(0);i++) { // cell
      std::cout << "   ";
      for(int j=0; j<pointField_.dimension(1);j++) {// IP
         std::cout << pointField_(i,j,0) << " " << pointField_(i,j,1) << std::endl;
      }
   }    

   for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
      PHX::MDField<ScalarT,panzer::Cell,panzer::IP,panzer::Dim> & field = scatterFields_[fieldIndex];
      std::cout << "Field \"" << field.fieldTag().name() << "\"" << std::endl;
      for(int i=0; i<field.dimension(0);i++) { // cell
         std::cout << "   ";
         for(int j=0; j<field.dimension(1);j++) {// IP
            std::cout << field(i,j,0) << " " << field(i,j,1) << std::endl;
         }
      }    
   }
*/

   for(int d=0;d<spatialDimension_;d++) {
      for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
         PHX::MDField<ScalarT,panzer::Cell,panzer::IP,panzer::Dim> & field = scatterFields_[fieldIndex];
         std::string fieldName = names_[fieldIndex]+dimStrings[d];

         std::vector<double> average(field.dimension(0),0.0);
         // write to double field
         for(int i=0; i<field.dimension(0);i++) { // cell
            for(int j=0; j<field.dimension(1);j++) // IP
               average[i] += field(i,j,d);         // insert dimension
            average[i] /= field.dimension(1);
         }
      
         // add in vector value at d^th point
         mesh_->setCellFieldData(fieldName,blockId,localCellIds,average);
      }
   }
}

} // end panzer_stk

#endif
