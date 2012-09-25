// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_PCE_SPECIALIZATIONS_HPP
#define KOKKOS_PCE_SPECIALIZATIONS_HPP

#include "Stokhos_Sacado.hpp"
#include "Kokkos_MultiVectorKernelOps.hpp"

// Specializations of several Kokkos kernels for PCE scalar types

namespace Kokkos {

  template <typename ScalarType, typename StorageType>
  struct WeightNormOp< Sacado::PCE::OrthogPoly<ScalarType,StorageType> > {
    typedef Sacado::PCE::OrthogPoly<ScalarType,StorageType> Scalar;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename SCT::innerProductType   ipType;
    const Scalar *x, *w;
    typedef  ipType ReductionType;
    inline static ipType KERNEL_PREFIX identity() {
      return Teuchos::ScalarTraits<ipType>::zero();
    }
    inline static ipType KERNEL_PREFIX reduce(ipType x, ipType y) {
      return x+y;
    }
    inline        ipType KERNEL_PREFIX generate(int i) {
      Scalar tmp = x[i] / w[i];
      return SCT::innerProduct(tmp,tmp);
    }
  };

  template <typename ScalarType, typename StorageType>
  struct DotOp1< Sacado::PCE::OrthogPoly<ScalarType,StorageType> > {
    typedef Sacado::PCE::OrthogPoly<ScalarType,StorageType> Scalar;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename SCT::innerProductType   ipType;
    const Scalar *x;
    typedef  ipType ReductionType;
    inline static ipType KERNEL_PREFIX identity() {
      return Teuchos::ScalarTraits<ipType>::zero();
    }
    inline static ipType KERNEL_PREFIX reduce(ipType x, ipType y) {
      return x+y;
    }
    inline        ipType KERNEL_PREFIX generate(int i) {
      Scalar xi = x[i]; 
      return SCT::innerProduct(xi,xi);
    }
  };

  template <typename ScalarType, typename StorageType>
  struct DotOp2< Sacado::PCE::OrthogPoly<ScalarType,StorageType> > {
    typedef Sacado::PCE::OrthogPoly<ScalarType,StorageType> Scalar;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename SCT::innerProductType ipType;
    const Scalar *x, *y;
    typedef ipType ReductionType;
     inline static ipType KERNEL_PREFIX identity() {
      return Teuchos::ScalarTraits<ipType>::zero();
    }
    inline static ipType KERNEL_PREFIX reduce(ipType x, ipType y) {
      return x+y;
    }
    inline        ipType KERNEL_PREFIX generate(int i) {
      Scalar xi = x[i]; Scalar yi = y[i];
      return SCT::innerProduct(xi, yi);
    }
  };

}


#endif
