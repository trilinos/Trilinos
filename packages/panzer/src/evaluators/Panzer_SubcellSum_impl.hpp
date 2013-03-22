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

#ifndef __Panzer_SubcellSum_impl_hpp__
#define __Panzer_SubcellSum_impl_hpp__

#include "Panzer_PureBasis.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(SubcellSum,p) 
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const std::string inName = p.get<std::string>("Field Name");
  const std::string outName = p.get<std::string>("Sum Name");
  Teuchos::RCP<const PureBasis> basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  multiplier = p.get<double>("Multiplier");

  inField = PHX::MDField<ScalarT,Cell,BASIS>( inName, basis->functional);
  outField = PHX::MDField<ScalarT,Cell>( outName, basis->cell_data);

  this->addDependentField(inField);
  this->addEvaluatedField(outField);

  // build a field pattern object so that looking up closure indices is easy
  fieldPattern_ = Teuchos::rcp(new IntrepidFieldPattern(basis->getIntrepidBasis()));
    
  std::string n = "SubcellSum: " + outField.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(SubcellSum,sd,fm)
{
  this->utils.setFieldData(inField,fm);
  this->utils.setFieldData(outField,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(SubcellSum,workset)
{ 
  std::vector<int> indices;
 
  // figure out which indices to sum (this can be made more efficient by 
  // simply saving the indices and only updating if the subcell dimension
  // and index changes)
  fieldPattern_->getSubcellClosureIndices(workset.subcell_dim,workset.subcell_index,indices);

  for(std::size_t c=0;c<workset.num_cells;c++) {
    outField(c) = 0.0; // initialize field 

    // sum over all relevant indices for this subcell
    for(std::size_t i=0;i<indices.size();i++)
      outField(c) += inField(c,indices[i]);
 
    // scale by what ever the user wants
    outField(c) *= multiplier;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
SubcellSum<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Sum Name", "?");
  p->set<std::string>("Field Name", "?");
  p->set<double>("Multiplier",1.0);

  Teuchos::RCP<const panzer::PureBasis> basis;
  p->set("Basis", basis);

  return p;
}

//**********************************************************************

}

#endif
