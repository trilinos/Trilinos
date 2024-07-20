// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BasisInformation_HPP
#define PANZER_BasisInformation_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

  //! Description and data layouts associated with a particular basis
  class BasisInformation { 

  public:
    typedef enum { HGRAD=0, HCURL=1, HDIV=2, CONST=3 } EElementSpace;
    
    /** Build a basis information given a type and order
      \param[in] basis_type String name that describes the type of basis
      \param[in] basis_order Order of the basis
      \param[in] cell_topo A shards topology description
    */
    BasisInformation(const std::string & basis_type,const int basis_order,
              const shards::CellTopology & cell_topo);

    //! Returns the basis type
    std::string type() const
    { return basis_type_; }
    
    //! Returns the polynomial order of the basis
    int order() const
    { return basis_order_; }

    const shards::CellTopology & getCellTopology() const
    { return topology_; }

    //! Returns the dimension of the basis from the topology
    int dimension() const
    { return topology_.getDimension(); }
    
    EElementSpace getElementSpace() const
    { return element_space_; }

    bool requiresOrientations() const
    { return getElementSpace()==HCURL || getElementSpace()==HDIV || (getElementSpace()==HGRAD && basis_order_ > 2) ; }

    bool supportsGrad() const
    { return getElementSpace()==HGRAD; }

    bool supportsCurl() const
    { return getElementSpace()==HCURL; }

    bool supportsDiv() const
    { return getElementSpace()==HDIV; }

    bool isVectorBasis() const
    { return getElementSpace()==HCURL || getElementSpace()==HDIV; }

    bool isScalarBasis() const
    { return getElementSpace()==HGRAD || getElementSpace()==CONST; }

  private:

    shards::CellTopology topology_;

    std::string basis_type_;
    int basis_order_;

    EElementSpace element_space_;
  };

}

#endif
