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

#include "Panzer_config.hpp"

#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseFunctional_Aggregator_Builder::
build<panzer::Traits::Residual>() const
{ 
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}

#ifdef HAVE_STOKHOS
template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseFunctional_Aggregator_Builder::
build<panzer::Traits::SGResidual>() const
{ 
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::SGResidual,panzer::Traits>); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}
#endif

}

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_ResponseFunctional_Aggregator_impl.hpp"
#ifdef HAVE_STOKHOS
   #include "Panzer_ResponseFunctional_AggregatorSG_impl.hpp"
#endif 

template class panzer::ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>;

#ifdef HAVE_STOKHOS
template class panzer::ResponseFunctional_Aggregator<panzer::Traits::SGResidual,panzer::Traits>;
#endif

#endif
