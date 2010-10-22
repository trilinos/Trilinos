
#ifndef PANZER_INTEGRATION_RULE_HPP
#define PANZER_INTEGRATION_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"
#include "Shards_CellTopology.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {

  class CellData;

  struct IntegrationRule {
    
    //! if side = -1 then we use the cell volume integration rule.
    IntegrationRule(int cubature_degree, const panzer::CellData& cell_data);
  
    // Returns true if this Integration rule is for a sideset
    bool isSide();

    Teuchos::RCP<shards::CellTopology> topology;
    
    Teuchos::RCP<shards::CellTopology> side_topology;
    
    //! Data layout for scalar fields
    Teuchos::RCP<PHX::DataLayout> dl_scalar;
    //! Data layout for vector fields
    Teuchos::RCP<PHX::DataLayout> dl_vector;
    //! Data layout for rank-2 tensor fields
    Teuchos::RCP<PHX::DataLayout> dl_tensor;
    
    int cubature_degree;
    int spatial_dimension;
    int workset_size;
    int num_points;
    //! Defaults to -1 if this is volume and not sideset
    int side;

  };

}

#endif
