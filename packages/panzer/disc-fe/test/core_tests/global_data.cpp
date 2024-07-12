// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Panzer_GlobalData.hpp"

namespace panzer {

  class TestObject : public panzer::GlobalDataAcceptorDefaultImpl {

  public:

    TestObject() {}

  };

  class TestObject2 : public panzer::GlobalDataAcceptorDefaultImpl {

  public:

    TestObject2(const Teuchos::RCP<panzer::GlobalData>& gd) : 
      panzer::GlobalDataAcceptorDefaultImpl(gd)
    {
      
    }

  };
  
  TEUCHOS_UNIT_TEST(global_data, builder)
  {
    
    Teuchos::RCP<panzer::GlobalData> gd;

    TEST_ASSERT(is_null(gd));

    {
      gd = panzer::createGlobalData(false);

      TEST_ASSERT(nonnull(gd));
      TEST_ASSERT(nonnull(gd->pl));
      TEST_ASSERT(is_null(gd->os));
    }

    {
      gd = panzer::createGlobalData(true);

      TEST_ASSERT(nonnull(gd));
      TEST_ASSERT(nonnull(gd->pl));
      TEST_ASSERT(nonnull(gd->os));
    }
  }

  TEUCHOS_UNIT_TEST(global_data, accessor_default_impl)
  {
    
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    
    TestObject t;
    t.setGlobalData(gd);

    TestObject2 t2(gd);

    TEST_EQUALITY(gd.get(), t2.getGlobalData().get());

  }

}
