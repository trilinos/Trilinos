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
