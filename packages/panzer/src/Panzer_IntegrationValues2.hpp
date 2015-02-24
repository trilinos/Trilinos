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


#ifndef PANZER_INTEGRATION_VALUES2_HPP
#define PANZER_INTEGRATION_VALUES2_HPP

#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_Dimension.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

  template <typename Scalar>
  class IntegrationValues2 {
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;

    typedef PHX::MDField<Scalar> ArrayDynamic;
    typedef PHX::MDField<double> DblArrayDynamic;

    typedef PHX::MDField<Scalar,IP> Array_IP;
    typedef PHX::MDField<Scalar,IP,Dim> Array_IPDim;

    typedef PHX::MDField<Scalar,Cell,IP> Array_CellIP;
    typedef PHX::MDField<Scalar,Cell,IP,Dim> Array_CellIPDim;
    typedef PHX::MDField<Scalar,Cell,IP,Dim,Dim> Array_CellIPDimDim;

    typedef PHX::MDField<Scalar,Cell,BASIS,Dim> Array_CellBASISDim;

    IntegrationValues2(const std::string & pre="",bool allocArrays=false) 
        : alloc_arrays(allocArrays), prefix(pre), ddims_(1,0) {}
    
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    void setupArraysForNodeRule(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    //! Cell vertex coordinates, not basis coordinates
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates);

    Array_IPDim cub_points;              // <IP,Dim>
    Array_IPDim side_cub_points;         // <IP,Dim> points on face topology (dim-1)
    Array_IP cub_weights;                // <IP>
    Array_CellBASISDim node_coordinates; // <Cell,BASIS,Dim>
    Array_CellIPDimDim jac;              // <Cell,IP,Dim,Dim>
    Array_CellIPDimDim jac_inv;          // <Cell,IP,Dim,Dim>
    Array_CellIP jac_det;                // <Cell,IP>
    Array_CellIP weighted_measure;       // <Cell,IP>

    Teuchos::RCP<const panzer::IntegrationRule> int_rule;

    Teuchos::RCP< Intrepid::Cubature<double,DblArrayDynamic> > intrepid_cubature;

    // for Shakib stabilization <Cell,IP,Dim,Dim>
    Array_CellIPDimDim covarient; 
    Array_CellIPDimDim contravarient; 
    Array_CellIP norm_contravarient; 

    // integration points
    Array_CellIPDim ip_coordinates;      // <Cell,IP,Dim>

    DblArrayDynamic dyn_cub_points, dyn_side_cub_points, dyn_cub_weights;

  private:
    bool alloc_arrays;
    std::string prefix;
    std::vector<PHX::index_size_type> ddims_;
  };

} // namespace panzer

#endif
