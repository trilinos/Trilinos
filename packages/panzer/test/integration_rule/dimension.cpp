#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"

namespace panzer {

  using shards::Array;
  using shards::NaturalOrder;

  TEUCHOS_UNIT_TEST(dimension, default)
  {
    Array<double,NaturalOrder,Dim,IP,BASIS,NODE,Point,Cell,Dummy> a;
    Array<double,NaturalOrder,Dim,IP,BASIS,NODE,Point,Cell,Dummy> b;
    
    a = b;
  }

}
