
#ifndef PANZER_INTEGRATION_RULE_HPP
#define PANZER_INTEGRATION_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"
#include "Shards_CellTopology.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_PointRule.hpp"

#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {

  class CellData;

  class IntegrationRule : public PointRule {
  public:
    
    //! if side = -1 then we use the cell volume integration rule.
    IntegrationRule(int cubature_degree, const panzer::CellData& cell_data);

    void setup(int cubature_degree, const panzer::CellData& cell_data);
  
    int cubature_degree;

    //! print information about the integration rule
    virtual void print(std::ostream & os);
  
  private:

  };

}

#endif
