// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataFunctors.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/31/23.
//

#ifndef Intrepid2_DataFunctors_hpp
#define Intrepid2_DataFunctors_hpp

/** \file  Intrepid2_DataFunctors.hpp
   \brief  Defines functors for use with Data objects: so far, we include simple arithmetical functors for sum, difference, product, and quotient.
   \author Created by N.V. Roberts.
*/

namespace Intrepid2 {

template<class Scalar>
struct ScalarSumFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a + b;
  }
};

template<class Scalar>
struct ScalarDifferenceFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a - b;
  }
};

template<class Scalar>
struct ScalarProductFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a * b;
  }
};

template<class Scalar>
struct ScalarQuotientFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a / b;
  }
};

} // end namespace Intrepid2

#endif /* Intrepid2_DataFunctors_hpp */
