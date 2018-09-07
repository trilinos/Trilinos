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

#ifndef PANZER_EVALUATOR_PARAMETER_DECL_HPP
#define PANZER_EVALUATOR_PARAMETER_DECL_HPP

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_ParameterLibrary.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
  template <typename EvalT> class ScalarParameterEntry;

//! Constant parameter from sacado parameter library
  template<typename EvalT, typename TRAITS>
  class Parameter : 
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS> {
  
  public:
    
    Parameter(const std::string parameter_name,
	      const std::string field_name,
	      const Teuchos::RCP<PHX::DataLayout>& data_layout,
	      panzer::ParamLib& param_lib);
    
    void evaluateFields(typename TRAITS::EvalData ud);
    
  private:
    
    typedef typename EvalT::ScalarT ScalarT;
    
    PHX::MDField<ScalarT, Cell, Point> target_field;
    
    std::size_t cell_data_size;
    
    ScalarT initial_value;
    
    Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > param;
  };
  
}

#endif
