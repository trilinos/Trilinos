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

#ifndef __PoissonExample_EquationSet_Energy_hpp__
#define __PoissonExample_EquationSet_Energy_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace Example {

/** The equation set serves two roles. The first is to let the panzer library
  * know which fields this equation set defines and their names. It registers
  * the evaluators required for a particular equation set. The level of the
  * granularity is largely up to a user. For instance this could be the momentum
  * or continuity equation in Navier-Stokes, or it could simply be the Navier-Stokes
  * equations. 
  *
  * Generally, this inherits from the panzer::EquationSet_DefaultImpl which takes
  * care of adding the gather (extract basis coefficients from solution vector) and 
  * scatter (using element matrices and vectors distribute and sum their values
  * to a global linear system) evaluators. These use data members that must be set by
  * the user.
  */
template <typename EvalT>
class PoissonEquationSet : public panzer::EquationSet_DefaultImpl<EvalT> {
public:    

   /** In the constructor you set all the fields provided by this
     * equation set. 
     */
   PoissonEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		      const int& default_integration_order,
                      const panzer::CellData& cell_data,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      const bool build_transient_support);
    
   /** The specific evaluators are registered with the field manager argument.
     */
   void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					      const panzer::FieldLibrary& field_library,
                                              const Teuchos::ParameterList& user_data) const;

};

}

#include "Example_PoissonEquationSet_impl.hpp"

#endif
