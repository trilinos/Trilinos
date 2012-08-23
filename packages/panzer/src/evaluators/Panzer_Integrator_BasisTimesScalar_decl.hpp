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

#ifndef PANZER_EVALUATOR_BASISTIMESSCALAR_DECL_HPP
#define PANZER_EVALUATOR_BASISTIMESSCALAR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {

/** Computes the integral
  * \f[
    \int s(x) \phi(x) \, d x
    \f]
  * where \f$\phi\f$is the test function and \f$s\f$ is
  * the scalar quantity. The parameter list passed into the construction
  * is formatted as follows
    \verbatim
    <ParameterList>
       <Parameter name="Residual Name" type="string" value=(required)/>
       <Parameter name="Value Name" type="string" value="(required)"/>
       <Parameter name="Basis" type="RCP<BasisIRLayout>" value=(required)/>
       <Parameter name="IR" type="RCP<IntegrationRule>" value="(required)"/>
       <Parameter name="Multiplier" type="double" value="(required)"/>
       <Parameter name="Field Multipliers" type="RCP<const std::vector>" value=Null (default)/>
    </ParameterList>
    \endverbatim
  */
PHX_EVALUATOR_CLASS(Integrator_BasisTimesScalar)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
    
  PHX::MDField<ScalarT,Cell,IP> scalar;

  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;

  std::size_t num_qp;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  Intrepid::FieldContainer<ScalarT> tmp;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

PHX_EVALUATOR_CLASS_END

}

#endif
