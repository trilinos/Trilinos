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

#ifndef PANZER_EVALUATOR_NORMALS_DECL_HPP
#define PANZER_EVALUATOR_NORMALS_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
/** Compute normals on a particular side of an element.
  * By default the normals are normalized. A second option
  * would be for the normals to be unormalized values.
  
    <ParameterList name="Name" type="string" value="<Name to give to the normals field>"/>
    <ParameterList name="Side Id" type="int" value="<side id to use for computing normals>"/>
    <ParameterList name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
    <ParameterList name="Normalize" type="bool" value="true"/>
  
  * The Name used to define the normals field is specified by "Name"
  * and the data layout is defined by the dl_vector field in the IntegrationRule.
  * The side ID must be legitimate for this topology and will be used to
  * construct an outward facing normal. The normals will be normalized by
  * default.  However, if Normalize=false then the resulting normals will have
  * the determinant of the side jacobian built in thus they correspond to a
  * differential on the side.
  */
PHX_EVALUATOR_CLASS(Normals)

  int side_id;
  int quad_order, quad_index;

  std::size_t num_qp, num_dim;

  PHX::MDField<ScalarT,Cell,Point,Dim> normals;
  bool normalize;

public:
  // for testing purposes
  const PHX::FieldTag & getFieldTag() const 
  { return normals.fieldTag(); }

PHX_EVALUATOR_CLASS_END

}

#endif
