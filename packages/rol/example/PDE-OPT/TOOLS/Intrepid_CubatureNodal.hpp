// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureNodal.hpp
    \brief  Header file for the Intrepid::CubatureNodal class.
    \author Created by D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_NODAL_HPP
#define INTREPID_CUBATURE_NODAL_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_Cubature.hpp"
#include "Teuchos_Assert.hpp"


namespace Intrepid {

/** \class Intrepid::CubatureNodal
    \brief Defines nodal cubature on select cell topologies, often used for
           lower-order integration or underintegration.  The use of this
           cubature with P1 elements yields a diagonal mass matrix.
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubatureNodal : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {
  private:
  /** \brief  Base topology of the cells for which the cubature is defined. See  the Shards package
              http://trilinos.sandia.gov/packages/shards for definition of base cell topology.
  */
  shards::CellTopology cellTopo_;

  unsigned numPoints_;
  unsigned cellDim_;

  Scalar weightVal_;

  public:

  ~CubatureNodal() {}

  /** \brief Constructor.
       \param cellTopo             [in]        - Cell topology.
  */
  CubatureNodal(const shards::CellTopology & cellTopo);

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Vector containing the cubature points.
      \param cubWeights      [out]     - Vector of corresponding cubature weights.
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const override;

  /** \brief Returns cubature points and weights.
              Method for physical space cubature, throws an exception.

       \param cubPoints             [out]        - Array containing the cubature points.
       \param cubWeights            [out]        - Array of corresponding cubature weights.
       \param cellCoords             [in]        - Array of cell coordinates
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights,
                           ArrayPoint  & cellCoords) const override;

  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const override;

  /** \brief Returns dimension of integration domain.
  */
  virtual int getDimension() const override;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has the size of the degree_ vector.
  */
  virtual void getAccuracy(std::vector<int> & degree) const override;

}; // end class CubatureNodal


} // end namespace Intrepid


// include templated definitions
#include "Intrepid_CubatureNodalDef.hpp"

#endif
