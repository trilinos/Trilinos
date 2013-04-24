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

#ifndef __Panzer_SubcellSum_decl_hpp__
#define __Panzer_SubcellSum_decl_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_FieldPattern.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {
    
/** This performs a sum over all the fields limited to the subcell
  * specified in the workset. It is useful for computing high-order
  * surface integrals as responses. 
  *
  * The field specified with "Sum Name" will be dimensioned as the number
  * of cells in the workset. The "Field Name" object is dimension as the number
  * of cells by the number of basis functions specified by the "Basis" object.
  * The "Evaluate On Closure" indicates if the subcells are to use the closure
  * index (i.e. all subcells of lesser dimension contained within a subcell) or
  * simply sum on those fields on the subcell proper.

  \verbatim
    <ParameterList>
      <Parameter name="Sum Name" type="string" value="<Name to give to the summed field>"/>
      <Parameter name="Field Name" type="string" value="<Name of field to sum>"/>
      <Parameter name="Basis" type="RCP<const PureBasis>" value="<user specified PureBasis object>"/>
      <Parameter name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
      <Parameter name="Evaluate On Closure" type="bool" value="false"/>
    </ParameterList>
  \endverbatim
  */
PHX_EVALUATOR_CLASS(SubcellSum)
  
  PHX::MDField<ScalarT,Cell> outField;  // result
    
  PHX::MDField<ScalarT,Cell,BASIS> inField; // function to be integrated

  double multiplier;

public:

  const PHX::FieldTag & getFieldTag() const 
  { return outField.fieldTag(); }

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
 
  // This is used to lookup closure indices (local Ids that live on a subcell)
  Teuchos::RCP<const panzer::FieldPattern> fieldPattern_;
  
  // evalaute on the "closure" of the indicated sub-cells
  bool evaluateOnClosure_;

PHX_EVALUATOR_CLASS_END

}

#endif
