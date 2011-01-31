// $Id: Stokhos_Quadrature.hpp,v 1.4 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_Quadrature.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
#define STOKHOS_KL_PRODUCT_EIGENPAIR_HPP

#include "Teuchos_Array.hpp"
#include "Stokhos_KL_OneDEigenPair.hpp"

namespace Stokhos {

  namespace KL {

    //! Container for multi-dimensional product of 1-D eigenfunctions/values
    template <typename value_type>
    struct ProductEigenPair {

      //! Product eigenvalue
      value_type eig_val;

      //! Eigenpairs for each dimension
      Teuchos::Array< OneDEigenPair<value_type> > eig_pairs;

      //! Default constructor
      ProductEigenPair() : eig_val(0.0), eig_pairs() {}

      //! Constructor
      ProductEigenPair(
	const Teuchos::Array< OneDEigenPair<value_type> >& eig_pairs_) 
	: eig_val(1.0), eig_pairs(eig_pairs_) {
	std::size_t sz = eig_pairs.size();
	for (std::size_t i=0; i<sz; i++)
	  eig_val *= eig_pairs[i].eig_val;
      }

      //! Evaluate eigenfunction at a given point
      value_type evalEigenfunction(const Teuchos::Array<value_type>& x) const {
	value_type result = 1.0;
	std::size_t sz = eig_pairs.size();
	for (std::size_t i=0; i<sz; i++)
	  result *= eig_pairs[i].eig_func->evaluate(x[i]);
	return result;
      }
      
      //! Print eigenpair
      void print(std::ostream& os) const {
	os << eig_val << ", ";
	std::size_t sz = eig_pairs.size();
	for (std::size_t i=0; i<sz-1; i++) {
	  os << "(";
	  eig_pairs[i].eig_func->print(os);
	  os << ") * ";
	}
	os << "(";
	eig_pairs[eig_pairs.size()-1].eig_func->print(os);
	os << ")";
      }
    };
    
    template <typename value_type>
    std::ostream& 
    operator << (std::ostream& os, const ProductEigenPair<value_type>& f) {
      f.print(os);
      return os;
    }

    //! Predicate class for sorting product eigenfunctions based on eigenvalue
    template <typename value_type>
    struct ProductEigenPairGreater : 
      public std::binary_function<ProductEigenPair<value_type>,
				  ProductEigenPair<value_type>,
				  bool> {
      bool operator() (const ProductEigenPair<value_type>& a, 
		       const ProductEigenPair<value_type>& b) {
	return a.eig_val > b.eig_val;
      }
    }; // struct ProductEigenPairGreater

  } // namespace KL

} // namespace Stokhos

#endif // STOKHOS_KL_PRODUCT_EIGENPAIR_HPP
