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
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_MaterialModelEntry.hpp"
#include <sstream>

namespace panzer {

  TEUCHOS_UNIT_TEST(material_model_entry, no_params)
  {
    std::string factory_name_m1 = "one";
    std::string factory_name_m2 = "two";
    
    panzer::MaterialModelEntry m1(factory_name_m1);
    panzer::MaterialModelEntry m2(factory_name_m2);

    panzer::MaterialModelEntry copy_m1(m1);

    panzer::MaterialModelEntry copy_m2;
    copy_m2 = m2;

    TEST_EQUALITY(m1.factoryName(), factory_name_m1);
    TEST_EQUALITY(m2.factoryName(), factory_name_m2);

    TEST_EQUALITY(copy_m1.factoryName(), factory_name_m1);
    TEST_EQUALITY(copy_m2.factoryName(), factory_name_m2);

    TEST_ASSERT(m1 == copy_m1);
    TEST_ASSERT(m2 == copy_m2);
    TEST_ASSERT(m1 != m2);
    TEST_ASSERT(copy_m1 != copy_m2);

    std::stringstream s;
    s << m1;
    TEST_EQUALITY(s.str(),"Material Model Entry: one"); 
  }

  TEUCHOS_UNIT_TEST(material_model_entry, with_params)
  {
    std::string factory_name = "one";
    Teuchos::ParameterList p1;
    p1.set<double>("value", 1.0);
    Teuchos::ParameterList p2;
    p2.set<int>("value", 1);
    

    panzer::MaterialModelEntry m1(factory_name,p1);
    panzer::MaterialModelEntry m2(factory_name,p2);

    panzer::MaterialModelEntry copy_m1(m1);

    panzer::MaterialModelEntry copy_m2;
    copy_m2 = m2;

    TEST_EQUALITY(m1.factoryName(), factory_name);
    TEST_EQUALITY(m2.factoryName(), factory_name);

    TEST_EQUALITY(copy_m1.factoryName(), factory_name);
    TEST_EQUALITY(copy_m2.factoryName(), factory_name);

    TEST_EQUALITY(m1.params(), p1);
    TEST_EQUALITY(m2.params(), p2);
    TEST_INEQUALITY(m1.params(), p2);

    TEST_EQUALITY(copy_m1.params(), p1);
    TEST_EQUALITY(copy_m2.params(), p2);

    TEST_ASSERT(m1 == copy_m1);
    TEST_ASSERT(m2 == copy_m2);
    TEST_ASSERT(m1 != m2);
    TEST_ASSERT(copy_m1 != copy_m2);

  }


}
