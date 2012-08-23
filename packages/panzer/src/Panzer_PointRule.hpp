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

#ifndef PANZER_POINT_RULE_HPP
#define PANZER_POINT_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

  class CellData;

  /** Base class useful for constructing data layouts
    * for points on a reference cell.
    */
  class PointRule {
  public:
    
    /** if side = -1 then we use the cell as an reference frame
      *
      * \param[in] ptName Name of the point rule.
      * \param[in] np Number of points per cell
      * \param[in] cell_data Description of the cell
      */
    PointRule(const std::string & ptName,int np, const panzer::CellData& cell_data);

    //! Destructor (Satisfying the compiler)
    virtual ~PointRule() {}

    void setup(const std::string & ptName,int np, const panzer::CellData& cell_data);
  
    // Returns true if this point rule is for a sideset
    bool isSide() const;

    /** Get the name of this point rule.
      */
    const std::string & getName() const;

    Teuchos::RCP<const shards::CellTopology> topology;
    
    Teuchos::RCP<shards::CellTopology> side_topology;
    
    //! Data layout for scalar fields
    Teuchos::RCP<PHX::DataLayout> dl_scalar;
    //! Data layout for vector fields
    Teuchos::RCP<PHX::DataLayout> dl_vector;
    //! Data layout for rank-2 tensor fields
    Teuchos::RCP<PHX::DataLayout> dl_tensor;
    
    int spatial_dimension;
    int workset_size;
    int num_points;

    //! Defaults to -1 if this is volume and not sideset
    int side;

    //! print information about the integration rule
    virtual void print(std::ostream & os);
  
  protected:
    PointRule() : side(-1) {}

    /** Look up side topology for a cell_data object. Returns null if
      * cell data does not correspond to a side object.
      */
    static Teuchos::RCP<shards::CellTopology> getSideTopology(const CellData & cell_data);

  private:
    std::string point_name;
  };

}

#endif
