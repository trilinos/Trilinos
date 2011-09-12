// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_EVALUATOR_UTILITIES_H
#define PHX_EVALUATOR_UTILITIES_H

#include <vector>

#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_FieldManager.hpp"

namespace PHX {

  /*! @brief Utilities to hide templating in concrete Evaluators.
   
  */
  template<typename EvalT, typename Traits> 
  struct EvaluatorUtilities {
    
    template <typename DataT>
    void setFieldData(PHX::Field<DataT>& f, PHX::FieldManager<Traits>& fm) 
    {
      fm.template getFieldData<DataT,EvalT>(f);
    }

    template <typename DataT,
	      typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	      typename Tag4, typename Tag5, typename Tag6, typename Tag7>
    void setFieldData(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,
		      Tag6,Tag7>& f, PHX::FieldManager<Traits>& fm) 
    {
      fm.template getFieldData<DataT,EvalT>(f);
    }

  };
}

#endif
