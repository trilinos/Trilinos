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

#ifndef   __myClosureModelFactory_TemplateBuilder_hpp__
#define   __myClosureModelFactory_TemplateBuilder_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// Sacado
#include "Sacado_mpl_apply.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Files for this specific example.
#include "myClosureModelFactory.hpp"

/**
 * 	\brief A means of building the closure model factory.
 *
 * 	This class allows us to build our closure model factory.
 */
class MyClosureModelFactory_TemplateBuilder
{
  public:
    
    /**
     * 	\brief Build the closure model factory.
     *
     * 	\returns A pointer to our built closure model factory.
     */
    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase>
    build() const 
    {
      using panzer::ClosureModelFactoryBase;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_static_cast;
      RCP<MyClosureModelFactory<EvalT>> closureFactory =
        rcp(new MyClosureModelFactory<EvalT>);
      return rcp_static_cast<ClosureModelFactoryBase>(closureFactory);
    } // end of build()

}; // end of class MyClosureModelFactory_TemplateBuilder

#endif // __myClosureModelFactory_TemplateBuilder_hpp__
