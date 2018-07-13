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

#ifndef   __myBCStrategyFactory_hpp__
#define   __myBCStrategyFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Traits.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Files specifically for this example.
#include "myBCStrategy.hpp"

// Create an object that can build our boundary condition strategy.
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(MyBCStrategy, MyBCStrategy)

/**
 *  \brief Our boundary condition strategy factory.
 *
 *  This class is used to build any boundary conditions we'd like to support.
 */
class MyBCStrategyFactory
  :
  public panzer::BCStrategyFactory
{
  public:

    /**
     * 	\brief Build our boundary condition strategy.
     *
     * 	Create all of the various boundary conditions that we'll support.  For
     * 	our specific example, there's only one (MyBCStrategy), but you could
     * 	have as many as you like.
     *
     * 	\param[in] bc          A boundary condition corresponding to one of the
     * 	                       BCs specified in the input XML file.
     * 	\param[in] global_data The global data object, which stores the
     * 	                       parameter library and the default output stream.
     *
     * 	\returns An object that manages all the different boundary conditions
     * 	         we've constructed.
     */
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits>>
    buildBCStrategy(
      const panzer::BC&                       bc,
      const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {
      using panzer::BCStrategy_TemplateManager;
      using panzer::Traits;
      using std::logic_error;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // Create the BCStrategy_TemplateManager and add our BCStrategy to it.
      RCP<BCStrategy_TemplateManager<Traits>> bcs_tm =
        rcp(new BCStrategy_TemplateManager<Traits>);
      bool found(false);
      PANZER_BUILD_BCSTRATEGY_OBJECTS("MyBCStrategy", MyBCStrategy)
      TEUCHOS_TEST_FOR_EXCEPTION(not found, logic_error, "Error:  The "       \
        "BCStrategy called \"" << bc.strategy() << "\" is not a valid "       \
        "identifier in the BCStrategyFactory.  Either add a valid "           \
        "implementation to your factory or fix your input file.  The "        \
        "relevant boundary condition is:\n\n" << bc << std::endl)
      return bcs_tm;
    } // end of buildBCStrategy()

}; // end of class MyBCStrategyFactory

#endif // __myBCStrategyFactory_hpp__
