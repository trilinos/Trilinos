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

#ifndef STOKHOS_SMOLYAK_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_SMOLYAK_PSEUDO_SPECTRAL_OPERATOR_HPP

#include "Stokhos_PseudoSpectralOperator.hpp"
#include "Stokhos_SmolyakBasis.hpp"
#include "Stokhos_TensorProductPseudoSpectralOperator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  /*!
   * \brief An operator for building pseudo-spectral coefficients using
   * a sparse Smolyak construction.
   */
  template <typename ordinal_t, 
	    typename value_t, 
	    typename point_compare_type = 
	    typename DefaultPointCompare<ordinal_t,value_t>::type>
  class SmolyakPseudoSpectralOperator : 
    public PseudoSpectralOperator<ordinal_t,value_t,point_compare_type> {
  public:

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef PseudoSpectralOperator<ordinal_type,value_type,point_compare_type> base_type;
    typedef typename base_type::point_type point_type;    
    typedef typename base_type::point_set_type point_set_type;
    typedef typename base_type::point_map_type point_map_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::set_iterator set_iterator;
    typedef typename base_type::const_set_iterator const_set_iterator;

    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef TensorProductPseudoSpectralOperator<ordinal_type, value_type> operator_type;
    typedef Teuchos::Array< Teuchos::RCP<operator_type> > operator_set_type;
    
    //! Constructor
    template <typename coeff_compare_type,
	      typename point_growth_rule_type>
    SmolyakPseudoSpectralOperator(
      const SmolyakBasis<ordinal_type,value_type,coeff_compare_type>& smolyak_basis, 
      const point_growth_rule_type& point_growth_rule,
      bool use_smolyak_apply = true,
      bool use_pst = true,
      const point_compare_type& point_compare = point_compare_type());

    //! Destructor
    virtual ~SmolyakPseudoSpectralOperator() {}

    //! Number of points
    ordinal_type point_size() const { return points.size(); }

    //! Number of coefficients
    ordinal_type coeff_size() const { return coeff_sz; }

    //! Iterator to begining of point set
    iterator begin() { return point_map.begin(); }
    
    //! Iterator to end of point set
    iterator end() { return point_map.end(); }

    //! Iterator to begining of point set
    const_iterator begin() const { return point_map.begin(); }

    //! Iterator to end of point set
    const_iterator end() const { return point_map.end(); }

    //! Iterator to begining of point set
    set_iterator set_begin() { return points.begin(); }

    //! Iterator to end of point set
    set_iterator set_end() { return points.end(); }

    //! Iterator to begining of point set
    const_set_iterator set_begin() const { return points.begin(); }

    //! Iterator to end of point set
    const_set_iterator set_end() const { return points.end(); }

    //! Get point index for given point
    ordinal_type index(const point_type& point) const { 
      const_set_iterator it = points.find(point);
      TEUCHOS_TEST_FOR_EXCEPTION(
	it == points.end(), std::logic_error, "Invalid term " << point);
      return it->second.second;
    }

    //! Get point for given index
    const point_type& point(ordinal_type n) const { return point_map[n]; }

    //! Apply pseudo-spectral quadrature operator
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void apply(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans = false) const;

  protected:

    //! Apply Smolyak pseudo-spectral operator using direct quadrature
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void apply_direct(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans) const;

    //! Apply Smolyak pseudo-spectral operator using Smolyak formula
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void apply_smolyak(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans) const;

    void gather(
      const Teuchos::Array<ordinal_type>& map, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
      bool trans, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const;

    void scatter(
      const Teuchos::Array<ordinal_type>& map, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
      bool trans, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const;

  protected:

    //! Use Smolyak apply method
    bool use_smolyak;

    //! Number of coefficients
    ordinal_type coeff_sz;

    //! Smolyak sparse grid
    point_set_type points;

    //! Map index to sparse grid term
    point_map_type point_map;

    //! Smolyak coefficients
    Teuchos::Array<value_type> smolyak_coeffs;

    //! Tensor pseudospectral operators
    operator_set_type operators;

    //! Gather maps for each operator for Smolyak apply
    Teuchos::Array< Teuchos::Array<ordinal_type> > gather_maps;

    //! Scatter maps for each operator for Smolyak apply
    Teuchos::Array< Teuchos::Array<ordinal_type> > scatter_maps;

    //! Matrix mapping points to coefficients for direct apply
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> A;

  };

}

#include "Stokhos_SmolyakPseudoSpectralOperatorImp.hpp"

#endif
