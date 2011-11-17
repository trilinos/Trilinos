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

    TestObject2(const Teuchos::RCP<const panzer::GlobalData>& gd) : 
      panzer::GlobalDataAcceptorDefaultImpl(gd)
    {
      
    }

  };
  
  TEUCHOS_UNIT_TEST(global_data_accessor, default_impl)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    TestObject t;
    t.setGlobalData(gd);

    TestObject2 t2(gd);

    TEST_EQUALITY(gd.get(), t2.getGlobalData().get());

  }

}
