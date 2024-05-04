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

/** \file   Intrepid_CubatureDirect.hpp
    \brief  Header file for the Intrepid::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_DIRECT_HPP
#define INTREPID_CUBATURE_DIRECT_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Teuchos_Assert.hpp"


namespace Intrepid {

/** \class Intrepid::CubatureDirect
    \brief Defines direct cubature (integration) rules in Intrepid.

    Cubature template (rule) consists of cubature points and cubature weights.
    Intrepid provides a small collection of frequently used cubature rule templates
    for FEM reconstructions on simplices (edge, tri, tet) and pyramid cells.

    For quad, hex, and triprism rules, see tensor-product rules
    defined in the class CubatureTensor, and its derived classes.

    Cubature rules for simplices and the pyramid are stored in the
    <var>cubature_data_</var> array.

    All templates are defined on a reference cell and can be mapped to physical space
    cells by the methods available in the MultiCell class.
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubatureDirect : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {
  private:

  protected:

  /** \brief The degree of polynomials that are integrated
             exactly by this cubature rule.
  */
  int degree_;

  /** \brief Dimension of integration domain.
  */
  int dimension_;


  public:

  virtual ~CubatureDirect() {}

  CubatureDirect() {}  


  /** \brief Extracts points and weights from cubData.

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
      \param cubData         [in]      - Pointer to raw cubature data.
  */
  virtual void getCubatureData(ArrayPoint  &              cubPoints,
                               ArrayWeight &              cubWeights,
                               const CubatureTemplate * cubData) const;

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const;

  /** \brief Returns cubature points and weights.
              Method for physical space cubature, throws an exception.

       \param cubPoints             [out]        - Array containing the cubature points.
       \param cubWeights            [out]        - Array of corresponding cubature weights.
       \param cellCoords             [in]        - Array of cell coordinates
  */
  virtual void getCubature(ArrayPoint& cubPoints,
                           ArrayWeight& cubWeights,
                           ArrayPoint& cellCoords) const;

  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const;

  /** \brief Returns dimension of integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has size 1.
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const;

  /** \brief Returns cubature name.
  */
  virtual const char* getName() const = 0;

  /** \brief Exposes cubature data.
  */
  virtual const CubatureTemplate * exposeCubatureData() const = 0;

  /** \brief Returns maximum cubature accuracy.
  */
  virtual int getMaxAccuracy() const = 0;

}; // end class CubatureDirect 

} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureDirectDef.hpp>

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

