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

/** \file   Intrepid_CubatureTensor.hpp
    \brief  Header file for the Intrepid::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_TENSOR_HPP
#define INTREPID_CUBATURE_TENSOR_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Teuchos_TestForException.hpp"


namespace Intrepid {

/** \class Intrepid::CubatureTensor
    \brief Defines tensor-product cubature (integration) rules in Intrepid.
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubatureTensor : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {
  private:

  /** \brief Degree of polynomials that are integrated exactly by
             each cubature rule within the tensor product.
  */
  std::vector<int> degree_;

  /** \brief Dimension of integration domain.
  */
  int dimension_;

  /** \brief Array of cubature rules, stored as FieldContainers.
  */
  std::vector< Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > > cubatures_;
  
  public:

  ~CubatureTensor() {}

  /** \brief Constructor.

      \param cubatures        [in]     - Array of cubatures that represent the building blocks
                                         of the tensor product.
  */
  CubatureTensor( std::vector< Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > > cubatures);

  /** \brief Constructor.

      \param cubature1        [in]     - First direct cubature rule.
      \param cubature2        [in]     - Second direct cubature rule.
  */
  CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                 Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2);

  /** \brief Constructor.

      \param cubature1        [in]     - First direct cubature rule.
      \param cubature2        [in]     - Second direct cubature rule.
      \param cubature3        [in]     - Third direct cubature rule.
  */
  CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                 Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2,
                 Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature3);

  /** \brief Constructor.

      \param cubature         [in]     - Direct cubature rule.
      \param n                [in]     - Number of copies of the cubature rule in the tensor product.
  */
  CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature, int n);

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Vector containing the cubature points.
      \param cubWeights      [out]     - Vector of corresponding cubature weights.
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const;

  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const;

  /** \brief Returns dimension of integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has the size of the degree_ vector.
  */
  virtual void getAccuracy(std::vector<int> & degree) const;

}; // end class CubatureTensor 


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureTensorDef.hpp>

#endif
