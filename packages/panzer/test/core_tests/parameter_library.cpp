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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ParameterLibraryAcceptor_DefaultImpl.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(parameter_library, scalar_parameter_entry)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Residual
    { 
      const RCP<ScalarParameterEntry<panzer::Traits::Residual> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Residual>);
      
      const ScalarParameterEntry<panzer::Traits::Residual>::ScalarT value(5.0);

      entry->setValue(value);
      
      const panzer::Traits::Residual::ScalarT tol = 
	10.0*std::numeric_limits<panzer::Traits::Residual::ScalarT>::epsilon();
      const double double_tol = 10.0*std::numeric_limits<double>::epsilon();
      
      TEST_FLOATING_EQUALITY(entry->getValue(), value, tol);
      TEST_FLOATING_EQUALITY(entry->getRealValue(), 5.0, double_tol);

      const double double_value = 6.0;

      entry->setRealValue(double_value);
      
      TEST_FLOATING_EQUALITY(entry->getRealValue(), double_value, double_tol);
    }

    // Jacobian
    {
      RCP<ScalarParameterEntry<panzer::Traits::Jacobian> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Jacobian>);
      
      typedef Sacado::ScalarType<panzer::Traits::Jacobian::ScalarT>::type ValueType;

      ScalarParameterEntry<panzer::Traits::Jacobian>::ScalarT param(5.0);
      param.resize(1);
      ValueType d_val = ValueType(1.0);
      param.fastAccessDx(0) = d_val;

      entry->setValue(param);
      
      ValueType tol(10.0*std::numeric_limits<ValueType>::epsilon());
      const double double_tol = 10.0*std::numeric_limits<double>::epsilon();
      
      TEST_FLOATING_EQUALITY(entry->getValue().val(), param.val(), tol);
      TEST_FLOATING_EQUALITY(entry->getValue().fastAccessDx(0), d_val, tol);
      TEST_FLOATING_EQUALITY(entry->getRealValue(), 5.0, double_tol);

      const double double_value = 6.0;

      entry->setRealValue(double_value);
      
      TEST_FLOATING_EQUALITY(entry->getRealValue(), double_value, double_tol);
    }

  }

  TEUCHOS_UNIT_TEST(parameter_library, utilities)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using panzer::ParamLib;
    
    RCP<ParamLib> pl = rcp(new ParamLib);
    
    std::string name = "viscosity";

    panzer::createAndRegisterScalarParameter<panzer::Traits::Residual>(name,*pl);
    
    panzer::createAndRegisterScalarParameter<panzer::Traits::Jacobian>(name,*pl);
    
    const double value = 5.0;
    
    pl->setRealValue<panzer::Traits::Jacobian>(name,value);

    double test_value = pl->getRealValue<panzer::Traits::Jacobian>(name);

    const double tol = 10.0*std::numeric_limits<double>::epsilon();

    TEST_FLOATING_EQUALITY(test_value, value, tol);

  }

  class TestAcceptor : public panzer::ParameterLibraryAcceptor_DefaultImpl {

  public:

    TestAcceptor() {}

    TestAcceptor(const Teuchos::RCP<panzer::ParamLib>& pl) :
      panzer::ParameterLibraryAcceptor_DefaultImpl(pl)
    { }

  };


  TEUCHOS_UNIT_TEST(parameter_library, acceptor)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using panzer::ParamLib;

    RCP<ParamLib> p1 = rcp(new ParamLib);
    RCP<ParamLib> p2 = rcp(new ParamLib);

    TestAcceptor t1(p1);
    TestAcceptor t2;

    t2.setParameterLibrary(p2);

    TEST_EQUALITY(p1, t1.getParameterLibrary());
    TEST_EQUALITY(p2, t2.getParameterLibrary());
    TEST_INEQUALITY(p1, t2.getParameterLibrary());
    TEST_INEQUALITY(p2, t1.getParameterLibrary());

  }

}
