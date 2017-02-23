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

#ifndef   __myEquationSetFactory_hpp__
#define   __myEquationSetFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Files specifically for this example.
#include "myEquationSet.hpp"

// Description.                                                                  // JMG:  Fill this out.                          
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(MyEquationSet, MyEquationSet)

/**
 * 	\brief Description.                                                          // JMG:  Fill this out.                          
 *                                                                               //                                               
 * 	Detailed description.                                                        //                                               
 */
class
MyEquationSetFactory
  :
  public panzer::EquationSetFactory
{
  public:

    /**
     *  \brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     *  Detailed description.                                                    //                                               
     *                                                                           //                                               
     *  \param[?] params                    Description.                         //                                               
     *  \param[?] default_integration_order Description.                         //                                               
     *  \param[?] cell_data                 Description.                         //                                               
     *  \param[?] global_data               Description.                         //                                               
     *  \param[?] build_transient_support   Description.                         //                                               
     *                                                                           //                                               
     *  \returns Something.                                                      //                                               
     */
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits>>
    buildEquationSet(
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const int&                                  default_integration_order,
      const panzer::CellData&                     cell_data,
      const Teuchos::RCP<panzer::GlobalData>&     global_data,
      const bool                                  build_transient_support)
      const
    {
      using panzer::EquationSet_TemplateManager;
      using panzer::Traits;
      using std::logic_error;
      using std::string;
      using Teuchos::RCP;
      using Teuchos::rcp;
      RCP<EquationSet_TemplateManager<Traits>> eq_set =
        rcp(new EquationSet_TemplateManager<Traits>);
      bool found(false);
      PANZER_BUILD_EQSET_OBJECTS("MyEquationSet", MyEquationSet)
      if (not found)
      {
        string msg("Error - the \"Equation Set\" with \"Type\" = \"" +
          params->get<string>("Type") + "\" is not a valid equation set "     \
          "identifier.  Please supply the correct factory.\n");
        TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, msg);
      } // end if (not found)

      return eq_set;
    } // end of buildEquationSet()

}; // end of class MyEquationSetFactory

#endif //  __myEquationSetFactory_hpp__
