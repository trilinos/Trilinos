// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_ORTHOGPOLYAPPROX_HPP
#define STOKHOS_ORTHOGPOLYAPPROX_HPP

#include "Teuchos_RCP.hpp"             // class data member
#include "Stokhos_OrthogPolyBasis.hpp" // class data member
#include <ostream>	               // for std::ostream

#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*! 
   * \brief Class to store coefficients of a projection onto an orthogonal
   * polynomial basis.
   */
  template <typename ordinal_type, typename value_type, 
	    typename storage_type = Stokhos::StandardStorage<ordinal_type, 
							     value_type> >
  class OrthogPolyApprox {
  public:

    typedef typename storage_type::reference reference;
    typedef typename storage_type::const_reference const_reference;
    typedef typename storage_type::pointer pointer;
    typedef typename storage_type::const_pointer const_pointer;

    //! Constructor with supplied size \c sz
    /*!
     * Normally \c sz should equal the basis size, however this is not
     * enforced since other situations can arise, e.g., using a size of 1
     * for a constant expansion.  If \c sz = 0, it is computed from the basis
     * size (and is 1 if \c basis is null).
     */
    OrthogPolyApprox(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis = Teuchos::null,
      ordinal_type sz = 0,
      const value_type* vals = NULL);
    
    //! Copy constructor
    OrthogPolyApprox(const OrthogPolyApprox& x);
    
    //! Destructor
    ~OrthogPolyApprox();
    
    //! Assignment operator (deep copy)
    OrthogPolyApprox& operator=(const OrthogPolyApprox& x);

    //! Assignment operator with scalar
    OrthogPolyApprox& operator=(const value_type& v);

    //! Initialize coefficients to value
    void init(const value_type& v);

    //! Initialize coefficients to an array of values
    void init(const value_type* v);

    //! Initialize coefficients from an OrthogPolyApprox with different storage
    template <typename S>
    void init(const OrthogPolyApprox<ordinal_type, value_type, S>& v) {
      coeff_.init(v.coeff());
    }

    //! Load coefficients to an array of values
    void load(value_type* v);

    //! Load coefficients into an OrthogPolyApprox with different storage
    template <typename S>
    void load(OrthogPolyApprox<ordinal_type, value_type, S>& v) {
      coeff_.load(v.coeff());
    }

    //! Return basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > 
    basis() const;

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.  Coefficients are preserved.
     */
    void reset(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis, ordinal_type sz = 0);

    //! Resize coefficient array (coefficients are preserved)
    void resize(ordinal_type sz);

    //! Return size
    ordinal_type size() const;

    //! Return coefficient array
    pointer coeff();

    //! Return coefficient array
    const_pointer coeff() const;

    //! Array access
    reference operator[](ordinal_type i);

    //! Array access
    const_reference operator[](ordinal_type i) const;

    //! Get coefficient term for given dimension and order
    reference term(ordinal_type dimension, ordinal_type order);

    //! Get coefficient term for given dimension and order
    const_reference term(ordinal_type dimension, ordinal_type order) const;

    //! Get orders for a given term
    const MultiIndex<ordinal_type>& order(ordinal_type term) const;

    //! Evaluate polynomial approximation at a point
    value_type evaluate(const Teuchos::Array<value_type>& point) const;

    //! Evaluate polynomial approximation at a point with supplied basis values
    value_type evaluate(const Teuchos::Array<value_type>& point,
			const Teuchos::Array<value_type>& basis_vals) const;

    //! Compute mean of expansion
    value_type mean() const;

    //! Compute standard deviation of expansion
    value_type standard_deviation() const;

    //! Compute the two-norm of expansion
    value_type two_norm() const;

    //! Compute the squared two-norm of expansion
    value_type two_norm_squared() const;

    //! Compute the L2 inner product of 2 PCEs
    value_type inner_product(const OrthogPolyApprox& b) const;

    //! Print approximation in basis
    std::ostream& print(std::ostream& os) const;

  protected:

    //! Basis expansion is relative to
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis_;

    //! OrthogPolyApprox coefficients
    storage_type coeff_;
    
  }; // class OrthogPolyApprox

  //! Prints the array of coefficients (more compact than print())
  template <typename ordinal_type, typename value_type, typename node_type> 
  std::ostream& 
  operator << (std::ostream& os, 
	       const OrthogPolyApprox<ordinal_type,value_type,node_type>& a);

} // namespace Stokhos

#include "Stokhos_OrthogPolyApproxImp.hpp"

#endif //  STOKHOS_ORTHOGPOLYAPPROX_HPP
