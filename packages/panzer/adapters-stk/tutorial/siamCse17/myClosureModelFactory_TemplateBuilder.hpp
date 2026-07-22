// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
