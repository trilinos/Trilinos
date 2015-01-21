// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
#define STOKHOS_KL_PRODUCT_EIGENPAIR_HPP

#include "Teuchos_Array.hpp"

namespace Stokhos {

  namespace KL {

    //! Container for multi-dimensional product of 1-D eigenfunctions/values
    template <typename eigen_function_type, typename Device>
    struct ProductEigenPair {

      typedef typename eigen_function_type::value_type value_type;
      typedef Device device_type;
      typedef OneDEigenPair<eigen_function_type> one_d_eigen_pair_type;
      typedef Teuchos::Array<one_d_eigen_pair_type> eig_pair_type;

      //! Product eigenvalue
      value_type eig_val;

      //! Eigenpairs for each dimension
      eig_pair_type eig_pairs;

      //! Default constructor
      ProductEigenPair() : eig_val(0.0), eig_pairs() {}

      //! Copy constructor
      ProductEigenPair(const ProductEigenPair& ep) :
        eig_val(ep.eig_val), eig_pairs(ep.eig_pairs) {}

      //! Assignment
      ProductEigenPair& operator=(const ProductEigenPair& ep) {
        if (this != &ep) {
          eig_val = ep.eig_val;
          eig_pairs = ep.eig_pairs;
        }
        return *this;
      }

      //! Set eigen pairs
      void set(const Teuchos::Array<one_d_eigen_pair_type>& ep) {
        eig_val = 1.0;
        eig_pairs = ep;
        std::size_t sz = eig_pairs.size();
        for (std::size_t i=0; i<sz; i++)
          eig_val *= eig_pairs[i].eig_val;
      }

      //! Evaluate eigenfunction at a given point
      template <typename point_type>
      KOKKOS_INLINE_FUNCTION
      value_type evalEigenfunction(const point_type& x) const {
        value_type result = 1.0;
        std::size_t sz = eig_pairs.size();
        for (std::size_t i=0; i<sz; i++)
          result *= eig_pairs[i].eig_func.evaluate(x[i]);
        return result;
      }

      //! Print eigenpair
      void print(std::ostream& os) const {
        os << eig_val << ", ";
        std::size_t sz = eig_pairs.size();
        for (std::size_t i=0; i<sz-1; i++) {
          os << "(";
          eig_pairs[i].eig_func.print(os);
          os << ") * ";
        }
        os << "(";
        eig_pairs[eig_pairs.size()-1].eig_func.print(os);
        os << ")";
      }
    };

    template <typename E, typename D>
    std::ostream&
    operator << (std::ostream& os, const ProductEigenPair<E,D>& f) {
      f.print(os);
      return os;
    }

    //! Predicate class for sorting product eigenfunctions based on eigenvalue
    template <typename E, typename D>
    struct ProductEigenPairGreater :
      public std::binary_function<ProductEigenPair<E,D>,
                                  ProductEigenPair<E,D>,
                                  bool> {
      bool operator() (const ProductEigenPair<E,D>& a,
                       const ProductEigenPair<E,D>& b) {
        return a.eig_val > b.eig_val;
      }
    }; // struct ProductEigenPairGreater

  } // namespace KL

} // namespace Stokhos

#endif // STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
