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

#ifndef PANZER_POINT_VALUES_HPP
#define PANZER_POINT_VALUES_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_ArrayTraits.hpp"

namespace panzer {

  template <typename Scalar,typename Array>
  struct PointValues {
    typedef typename ArrayTraits<Scalar,Array>::size_type size_type;
    
    //! Sizes/allocates memory for arrays
    template <typename ArrayFactory>
    void setupArrays(const Teuchos::RCP<const panzer::PointRule>& pr,const ArrayFactory & af);

    template <typename NodeCoordinateArray,typename PointCoordinateArray>
    inline void evaluateValues(const NodeCoordinateArray & node_coordinates,const PointCoordinateArray & point_coordinates);

    template <typename CoordinateArray>
    void copyNodeCoords(const CoordinateArray& in_node_coords);

    template <typename CoordinateArray>
    void copyPointCoords(const CoordinateArray& in_point_coords);

    Array coords_ref;          // <Point,Dim>
    Array node_coordinates;    // <Cell,NODE,Dim>
    Array jac;                 // <Cell,Point,Dim,Dim>
    Array jac_inv;             // <Cell,Point,Dim,Dim>
    Array jac_det;             // <Cell,Point>

    // cell points
    Array point_coords;      // <Cell,Point,Dim>

    Teuchos::RCP<const panzer::PointRule> point_rule;
  };

} // namespace panzer

#include "Panzer_PointValues_impl.hpp"

#endif
