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

#ifndef PANZER_POINT_VALUES2_HPP
#define PANZER_POINT_VALUES2_HPP

#include "PanzerDiscFE_config.hpp"

#include "Panzer_PointRule.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_Dimension.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

  template <typename Scalar>
  class PointValues2 {
  public:
    typedef typename ArrayTraits<Scalar, PHX::MDField<Scalar> >::size_type size_type;

    template<typename SourceScalar>
    PointValues2<Scalar>&
    operator=(const PointValues2<SourceScalar>& source);

    PointValues2(const std::string & pre="",
                 bool allocArrays=false)
       : alloc_arrays_(allocArrays), prefix_(pre) {}

    PointValues2(const std::string & pre,
                 const std::vector<PHX::index_size_type> & ddims,
                 bool allocArrays=false)
       : alloc_arrays_(allocArrays), prefix_(pre), ddims_(ddims) {}
    
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::PointRule>& pr);

    /** Evaluate teh jacobian and derivative information at the requested reference
      * points.
      *
      * \param[in] node_coords Cell vertices
      * \param[in] point_coords Reference cell coordinates
      */
    template <typename CoordinateArray,typename PointArray>
    void evaluateValues(const CoordinateArray & node_coords,
                        const PointArray & in_point_coords)
    { copyNodeCoords(node_coords);
      copyPointCoords(in_point_coords);
      evaluateValues(); }

    /** Evaluate teh jacobian and derivative information at the requested reference
      * points. This version allows a shallow copy of the vertex coordinates. 
      *
      * \param[in] node_coords Cell vertices
      * \param[in] point_coords Reference cell coordinates
      * \param[in] shallow_copy_nodes Enable or disable a shallow copy of the vertices
      */ 
    template <typename PointArray>
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & node_coords,
                        const PointArray & in_point_coords, 
                        bool shallow_copy_nodes)
    { if(shallow_copy_nodes)
        node_coordinates = node_coords;
      else
        copyNodeCoords(node_coords);
      copyPointCoords(in_point_coords);
      evaluateValues(); }

    //! Return reference cell coordinates this class uses (IP,Dim) sized
    PHX::MDField<Scalar,IP,Dim> & getRefCoordinates() const 
    { return coords_ref; }

    //! Return the vertex coordinates this class uses (Cell,NODE,Dim) sized
    PHX::MDField<Scalar,Cell,NODE,Dim> & getVertexCoordinates() const
    { return node_coordinates; }

    // input fields: both mutable because of getRefCoordinates/getVertexCoordinates
    //               Not sure this is the best design, but works for this iteration
    mutable PHX::MDField<Scalar,IP,Dim>        coords_ref;       // <IP,Dim>
    mutable PHX::MDField<Scalar,Cell,NODE,Dim> node_coordinates; // <Cell,NODE,Dim>

    // output fields
    PHX::MDField<Scalar,Cell,IP,Dim,Dim>       jac;              // <Cell,IP,Dim,Dim>
    PHX::MDField<Scalar,Cell,IP,Dim,Dim>       jac_inv;          // <Cell,IP,Dim,Dim>
    PHX::MDField<Scalar,Cell,IP>               jac_det;          // <Cell,IP>
    PHX::MDField<Scalar,Cell,IP,Dim>           point_coords;     // <Cell,IP,Dim> // cell points

    Teuchos::RCP<const panzer::PointRule> point_rule;

  private:
    void evaluateValues();

    template <typename CoordinateArray>
    void copyNodeCoords(const CoordinateArray& in_node_coords);

    template <typename CoordinateArray>
    void copyPointCoords(const CoordinateArray& in_point_coords);

    bool alloc_arrays_;
    std::string prefix_;
    std::vector<PHX::index_size_type> ddims_;

  };

  template <typename Scalar>
  template<typename SourceScalar>
  PointValues2<Scalar>&
  PointValues2<Scalar>::operator=(const PointValues2<SourceScalar>& source)
  {
    // The separate template parameter for SourceScalar allows for
    // assignment to a "const Scalar" from a non-const "Scalar", but
    // we still need to enforce that the underlying scalar type is the
    // same.
    static_assert(std::is_same<typename std::decay<Scalar>::type,typename std::decay<SourceScalar>::type>::value,
                  "ERROR: PointValues assignment requires consistent scalar types!");

    coords_ref = source.coords_ref;
    node_coordinates = source.node_coordinates;
    jac = source.jac;
    jac_inv = source.jac_inv;
    jac_det = source.jac_det;
    point_coords = source.point_coords;
    point_rule = source.point_rule;
    return *this;
  }

} // namespace panzer

#endif
