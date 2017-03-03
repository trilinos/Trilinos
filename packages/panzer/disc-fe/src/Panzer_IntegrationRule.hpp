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
