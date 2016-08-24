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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_PureBasis.hpp"

// Evaluators
#include "Panzer_Constant.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::BCStrategy_ConstantDirichlet<EvalT>::
BCStrategy_ConstantDirichlet(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) 
  : panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc,global_data)
// the constructor called here stores bc and global_data
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Constant");
}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_ConstantDirichlet<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& user_data)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  this->addDOF(this->m_bc.equationSetName());              // DOF Name

  // add in the targert
  this->addTarget("Constant_"+this->m_bc.equationSetName(),  // Target Name
                  this->m_bc.equationSetName(),              // DOF Name
                  "Residual_"+this->m_bc.identifier());      // Residual Name
 
  // find the basis for this dof 
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();
  
  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
       dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    if (dof_it->first == this->m_bc.equationSetName())
      this->basis = dof_it->second;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
		     "Error the name \"" << this->m_bc.equationSetName()
		     << "\" is not a valid DOF for the boundary condition:\n"
		     << this->m_bc << "\n");
}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_ConstantDirichlet<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
			   const Teuchos::ParameterList& models,
			   const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
 
  // THIS LINE GETS PARAMETER "VALUE" FROM THE INPUT DECK..
  //  p.set("Value", this->m_bc.params()->template get<double>("Value"));

  // provide a constant target value to map into residual
  {
    ParameterList p("BC Constant Dirichlet");
    p.set("Name", "Constant_" + this->m_bc.equationSetName());
    p.set("Data Layout", basis->functional);
    p.set("Value", 2.9);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }
}
