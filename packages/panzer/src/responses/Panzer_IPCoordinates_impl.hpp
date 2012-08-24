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

#ifndef PANZER_EVALUATOR_IP_COORDINATES_IMPL_HPP
#define PANZER_EVALUATOR_IP_COORDINATES_IMPL_HPP

#include "Panzer_Workset_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_Assert.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
IPCoordinates<EvalT, Traits>::
IPCoordinates(int in_ir_order,
	      const Teuchos::RCP<std::string>& in_block_id,
	      const Teuchos::RCP<std::vector<ScalarT> >& in_coords) :
  first_evaluation(true),
  ir_order(in_ir_order),
  block_id(in_block_id),
  coords(in_coords)
{ 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::MDALayout;
  using PHX::MDField;

  TEUCHOS_ASSERT(nonnull(coords));

  RCP<MDALayout<Cell,IP,Dim> > dl = rcp(new MDALayout<Cell,IP,Dim>(0,0,0));
  dummy_field = PHX::MDField<ScalarT, Cell,IP,Dim>("IP Coordinates", dl);
  
  this->addEvaluatedField(dummy_field);
 
  this->setName("IPCoordinates Evaluator");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dummy_field,fm);

  ir_index = panzer::getIntegrationRuleIndex(ir_order,(*sd.worksets_)[0]);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::preEvaluate(typename Traits::PreEvalData data)
{
  
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  if (first_evaluation) {

    *block_id = workset.block_id;

    Intrepid::FieldContainer<double>& workset_coords = (workset.int_rules[ir_index])->ip_coordinates;

    if (tmp_coords.size() != Teuchos::as<std::size_t>(workset_coords.dimension(2)))
      tmp_coords.resize(workset_coords.dimension(2));

    // This ordering is for the DataTransferKit.  It blocks all x
    // coordinates for a set of points, then all y coordinates and if
    // required all z coordinates.
    for (int dim = 0; dim < workset_coords.dimension(2); ++dim) {
      for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
	for (int ip = 0; ip < workset_coords.dimension(1); ++ip) {
	  tmp_coords[dim].push_back(workset_coords(static_cast<int>(cell),ip,dim));
	}
      }
    }

  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::postEvaluate(typename Traits::PostEvalData data)
{
  if (first_evaluation) {
    coords->clear();
    for (std::size_t dim = 0; dim < tmp_coords.size(); ++dim) {
      for (typename std::vector<ScalarT>::const_iterator x=tmp_coords[dim].begin(); x != tmp_coords[dim].end(); ++ x)
	coords->push_back(*x);
    }
    tmp_coords.clear();
    first_evaluation = false;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
const PHX::MDField<typename IPCoordinates<EvalT, Traits>::ScalarT,Cell,IP,Dim> 
IPCoordinates<EvalT, Traits>::getEvaluatedField() const
{
  return dummy_field;
}

//**********************************************************************

}

#endif
