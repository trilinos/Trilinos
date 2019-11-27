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
#ifndef __Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian_impl_hpp__
#define __Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// blocked Tpetra scatter dirichlet file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & /* indexer */,
                                      const Teuchos::ParameterList& p)
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  this->addEvaluatedField(*scatterHolder_);

  this->setName(scatterName+" Scatter Dirichlet Residual BlockedTpetra (Hessian)");
}
  
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* vm */)
{
}

template <typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData /* d */)
{
}
  
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData /* workset */)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "ScatterDirichletResidual_BlockedTpetra<Hessian> is not yet implemented"); // just in case
}

}

// **************************************************************
#endif

#endif
