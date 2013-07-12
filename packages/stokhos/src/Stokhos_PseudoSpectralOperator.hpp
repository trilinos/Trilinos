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

#ifndef STOKHOS_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_PSEUDO_SPECTRAL_OPERATOR_HPP

#include <map>
#include "Teuchos_Array.hpp"

#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  //! Struct defining default point compare type
  template <typename ordinal_type, typename value_type>
  struct DefaultPointCompare {
    typedef LexographicLess< TensorProductElement<ordinal_type,value_type>, 
			     FloatingPointLess<value_type> > type;
  };

  /*!
   * \brief An operator interface for building pseudo-spectral approximations
   */
  template <typename ordinal_t, 
	    typename value_t, 
	    typename point_compare_type = 
	    typename DefaultPointCompare<ordinal_t,value_t>::type>
  class PseudoSpectralOperator {
  public:

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef TensorProductElement<ordinal_type,value_type> point_type;
    typedef std::map<point_type, std::pair<value_type,ordinal_type>, 
		     point_compare_type> point_set_type;
    typedef Teuchos::Array<point_type> point_map_type;
    
    typedef typename point_map_type::iterator iterator;
    typedef typename point_map_type::const_iterator const_iterator;
    typedef typename point_set_type::iterator set_iterator;
    typedef typename point_set_type::const_iterator const_set_iterator;

    //! Constructor
    PseudoSpectralOperator() {}

    //! Destructor
    virtual ~PseudoSpectralOperator() {}

    //! Number of points
    virtual ordinal_type point_size() const = 0;

    //! Number of coefficients
    virtual ordinal_type coeff_size() const = 0;

    //! Iterator to begining of point set
    virtual iterator begin() = 0;

    //! Iterator to end of point set
    virtual iterator end() = 0;

    //! Iterator to begining of point set
    virtual const_iterator begin() const = 0;

    //! Iterator to end of point set
    virtual const_iterator end() const = 0;

    //! Iterator to begining of point set
    virtual set_iterator set_begin() = 0;

    //! Iterator to end of point set
    virtual set_iterator set_end() = 0;

    //! Iterator to begining of point set
    virtual const_set_iterator set_begin() const = 0;

    //! Iterator to end of point set
    virtual const_set_iterator set_end() const = 0;

    //! Get point index for given point
    virtual ordinal_type index(const point_type& point) const = 0;

    //! Get point for given index
    virtual const point_type& point(ordinal_type n) const  = 0;

    //! Transform values at quadrature points to PCE coefficients
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    virtual void transformQP2PCE(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans = false) const = 0;

    //! Transform PCE coefficients to quadrature values
    /*!
     * \c input is a vector storing polynomial chaos coefficients
     * and \c result will contain the resulting values at the quadrature points.
     * \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    virtual void transformPCE2QP(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans = false) const = 0;

  };

}

#endif
