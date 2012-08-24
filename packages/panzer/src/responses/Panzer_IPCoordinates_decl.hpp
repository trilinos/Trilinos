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

#ifndef PANZER_EVALUATOR_IP_COORDINATES_DECL_HPP
#define PANZER_EVALUATOR_IP_COORDINATES_DECL_HPP

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {
  
  /** Evaluates the coordinates on fist evaluation
   *
   *  NOTE: This assumes a static mesh and does not recompute the
   *  coordinates on each response evaluation.  It is simple to change
   *  in the future if needed.
   *
   *  NOTE: This ordering of points is blocked by dimension.  All x
   *  points values, then all y point values and then if required all
   *  z point values.  This is ordering is done to support
   *  DataTransferKit requirements.
   */
  template<typename EvalT, typename Traits>
  class IPCoordinates : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
  
  public:
    
    typedef typename EvalT::ScalarT ScalarT;
    
    IPCoordinates(int ir_order,
		  const Teuchos::RCP<std::string>& block_id,
		  const Teuchos::RCP<std::vector<ScalarT> >& coords);
    
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& vm);
    
    void preEvaluate(typename Traits::PreEvalData data);
    void evaluateFields(typename Traits::EvalData ud);
    void postEvaluate(typename Traits::PostEvalData data);

    const PHX::MDField<ScalarT,Cell,IP,Dim> getEvaluatedField() const;

  private:
    
    PHX::MDField<ScalarT,Cell,IP,Dim> dummy_field;
 
    bool first_evaluation;

    //! integration rule order
    int ir_order;

    //! integration rule index
    int ir_index;

    //! Block name from ResponseData object
    Teuchos::RCP<std::string> block_id;

    //! Coordinate vector from ResponseData object
    Teuchos::RCP<std::vector<ScalarT> > coords;

    //! Temporary coordinates blocked by dimension
    std::vector<std::vector<ScalarT> > tmp_coords;
  };
  
}

#endif
