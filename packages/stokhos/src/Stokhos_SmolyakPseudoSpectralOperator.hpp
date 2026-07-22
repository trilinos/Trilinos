// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SMOLYAK_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_SMOLYAK_PSEUDO_SPECTRAL_OPERATOR_HPP

#include "Stokhos_PseudoSpectralOperator.hpp"
#include "Stokhos_SmolyakBasis.hpp"
#include "Stokhos_TensorProductPseudoSpectralOperator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

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
    template <typename coeff_compare_type>
    SmolyakPseudoSpectralOperator(
      const SmolyakBasis<ordinal_type,value_type,coeff_compare_type>& smolyak_basis, 
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

    //! Transform values at quadrature points to PCE coefficients
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void transformQP2PCE(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans = false) const;

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
      bool trans = false) const;

  protected:

    //! Apply transformation operator using direct method
    void apply_direct(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A,
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans) const;

    /*! 
     * \brief Transform values at quadrature points to PCE coefficients 
     * using Smolyak formula
     */
    void transformQP2PCE_smolyak(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans) const;

    /*! 
     * \brief Transform PCE coefficients to values at quadrature points
     * using Smolyak formula
     */
    void transformPCE2QP_smolyak(
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

    //! Matrix mapping points to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> qp2pce;

    //! Matrix mapping coefficients to points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> pce2qp;

    //! BLAS wrappers
    Teuchos::BLAS<ordinal_type,value_type> blas;

  };

}

#include "Stokhos_SmolyakPseudoSpectralOperatorImp.hpp"

#endif
