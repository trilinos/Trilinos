#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Teuchos_ParameterList.hpp>
#include "Panzer_Constant.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include <sstream>

namespace panzer {

  TEUCHOS_UNIT_TEST(evaluators, constant)
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using panzer::Dim;
    using panzer::Cell;
    using panzer::Point;
    
    RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));

    ParameterList p("Constant Test");
    p.set("Value", 4.0);
    p.set("Name", "Viscosity");
    p.set("Data Layout", dl);

    panzer::Constant<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::Constant<panzer::Traits::Jacobian,panzer::Traits> e_J(p);

  }

}
