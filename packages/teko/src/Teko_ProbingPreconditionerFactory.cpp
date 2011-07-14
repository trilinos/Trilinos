/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "Teko_ProbingPreconditionerFactory.hpp"

#ifdef Teko_ENABLE_Isorropia

#include "Teko_EpetraOperatorWrapper.hpp"

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_CrsMatrix.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;
using Teuchos::RCP;

namespace Teko {

/*****************************************************/


ProbingPreconditionerFactory::ProbingPreconditionerFactory()
{
   prober = rcp(new Isorropia::Epetra::Prober);
}

LinearOp ProbingPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   // make an epetra operator to be probed
   RCP<Epetra_Operator> epetraLo = rcp(new Teko::Epetra::EpetraOperatorWrapper(lo));

   // build color scheme
   prober->color();

   // probe operator: take me to your leader
   RCP<Epetra_CrsMatrix> retOp = prober->probe(*epetraLo);
   Teko::LinearOp probedOp = Thyra::epetraLinearOp(retOp);

   return Teko::buildInverse(*invFactory_,probedOp);
}

void ProbingPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   const std::string inverse_type = "Inverse Type";
   const std::string probing_graph_operator = "Probing Graph Operator";
   const std::string probing_graph = "Probing Graph";
   const std::string user_graph = "User Will Set Probing Graph";

   // get string specifying default inverse
   std::string invStr ="Amesos"; 
   if(pl.isParameter(inverse_type))
      invStr = pl.get<std::string>(inverse_type);

   if(pl.isParameter(probing_graph_operator))
      setGraphOperator(pl.get<Teko::LinearOp>(probing_graph_operator));
   else if(pl.isParameter(probing_graph))
      setGraph(pl.get<RCP<const Epetra_CrsGraph> >(probing_graph));
   else if(pl.isParameter(user_graph) && pl.get<bool>("User Will Set Probing Graph")){
     //noop
   }	  
   else {
      Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
      rh->preRequest<RCP<const Epetra_CrsGraph> >(Teko::RequestMesg("Probing Graph"));
      setGraph(rh->request<RCP<const Epetra_CrsGraph> >(Teko::RequestMesg("Probing Graph")));
   }

   setInverseFactory(invLib->getInverseFactory(invStr));
}

void ProbingPreconditionerFactory::setGraphOperator(const Teko::LinearOp & graphOp)
{
   RCP<const Epetra_CrsMatrix> crsMatrix = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*graphOp));
   setGraph(Teuchos::rcpFromRef(crsMatrix->Graph()));
}

void ProbingPreconditionerFactory::setGraph(const Teuchos::RCP<const Epetra_CrsGraph> & graph)
{
   prober->setGraph(graph);
}

void ProbingPreconditionerFactory::setProberList(const Teuchos::ParameterList & list)
{
   prober->setList(list);
}

} // end namespace Teko

#endif
