// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
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
