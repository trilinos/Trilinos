// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
#define STOKHOS_KL_PRODUCT_EIGENPAIR_HPP

#include "Teuchos_Array.hpp"

namespace Stokhos {

  namespace KL {

    //! Container for multi-dimensional product of 1-D eigenfunctions/values
    template <typename eigen_function_type, typename ExecutionSpace>
    struct ProductEigenPair {

      typedef typename eigen_function_type::value_type value_type;
      typedef ExecutionSpace execution_space;
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
    struct ProductEigenPairGreater
    {
      bool operator() (const ProductEigenPair<E,D>& a,
                       const ProductEigenPair<E,D>& b) {
        return a.eig_val > b.eig_val;
      }
    }; // struct ProductEigenPairGreater

  } // namespace KL

} // namespace Stokhos

#endif // STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
