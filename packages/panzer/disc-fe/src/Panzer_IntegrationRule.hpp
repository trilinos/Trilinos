// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_INTEGRATION_RULE_HPP
#define PANZER_INTEGRATION_RULE_HPP

#include "Teuchos_RCP.hpp"
//#include "Teuchos_ArrayRCP.hpp"
#include "Shards_CellTopology.hpp"
//#include "Phalanx_DataLayout.hpp"

#include "Panzer_PointRule.hpp"
//
//#include "Intrepid2_DefaultCubatureFactory.hpp"
//#include "Intrepid2_CubatureControlVolume.hpp"
//#include "Intrepid2_CubatureControlVolumeSide.hpp"
//#include "Intrepid2_CubatureControlVolumeBoundary.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Panzer_IntegrationDescriptor.hpp"

#include <ostream>
#include <string>

namespace panzer {

  class CellData;

  /** Derived class for building a point rule based on a finite element integration rule.

      \param[in] cubature_degree Order of the cubature integration.
      \param[in] cell_data Description of the cell.
   */
  class IntegrationRule : public PointRule, public IntegrationDescriptor {
  public:
    
    //! if side = -1 then we use the cell volume integration rule.
    //TEUCHOS_DEPRECATED
    IntegrationRule(int cubature_degree, const panzer::CellData& cell_data);

    //TEUCHOS_DEPRECATED
    IntegrationRule(const panzer::CellData& cell_data, const std::string & cv_type);

    IntegrationRule(const panzer::IntegrationDescriptor& description,
                    const Teuchos::RCP<const shards::CellTopology> & cell_topology,
                    const int num_cells,
                    const int num_faces=-1);

    // TODO: Move to protected
    void setup(int cubature_degree, const panzer::CellData& cell_data);

    // TODO: Move to protected
    void setup_cv(const panzer::CellData& cell_data, std::string cv_type);
  
    //! Returns the order of integration (cubature degree in intrepid lingo)
    // Use getOrder() from base class
    //TEUCHOS_DEPRECATED
    int order() const;

    // Use _cubature_order if inside class, use getOrder() if outside
    //TEUCHOS_DEPRECATED
    int cubature_degree;

    // Use _type if inside class, use getType() if outside class
    //TEUCHOS_DEPRECATED
    std::string cv_type;

    //! print information about the integration rule
    virtual void print(std::ostream & os);

    //! Construct an array containing the reference coordinates 
    void referenceCoordinates(Kokkos::DynRankView<double,PHX::Device> & container);

    //! Returns the integration point offset for a given subcell_index (i.e. local face index)
    int getPointOffset(const int subcell_index) const;
  
  protected:

    /// Setup a surface integration
    void setup_surface(const Teuchos::RCP<const shards::CellTopology> & cell_topology,
                       const int num_cells,
                       const int num_faces);

    std::vector<int> _point_offsets;

  private:

  };

}

#endif
