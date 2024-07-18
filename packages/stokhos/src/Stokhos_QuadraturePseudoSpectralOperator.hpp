// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_QUADRATURE_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_QUADRATURE_PSEUDO_SPECTRAL_OPERATOR_HPP

#include "Stokhos_PseudoSpectralOperator.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

namespace Stokhos {

  /*!
   * \brief An operator for building pseudo-spectral coefficients using
   * an arbitrary quadrature rule.
   */
  template <typename ordinal_t,
            typename value_t,
            typename point_compare_type =
            typename DefaultPointCompare<ordinal_t,value_t>::type>
  class QuadraturePseudoSpectralOperator :
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

    //! Constructor
    QuadraturePseudoSpectralOperator(
      const OrthogPolyBasis<ordinal_type,value_type>& basis,
      const Quadrature<ordinal_type,value_type>& quad,
      const point_compare_type& point_compare = point_compare_type()) :
      coeff_sz(basis.size()),
      points(point_compare) {

      const Teuchos::Array<value_type>& quad_weights = quad.getQuadWeights();
      const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
        quad.getQuadPoints();
      const Teuchos::Array< Teuchos::Array<value_type> > & quad_vals =
        quad.getBasisAtQuadPoints();

      ordinal_type nqp = quad.size();
      ordinal_type npc = basis.size();
      ordinal_type dim = basis.dimension();

      // Generate point set
      point_type thePoint(dim);
      for (ordinal_type i=0; i<nqp; ++i) {
        for (ordinal_type k=0; k<dim; k++) {
          thePoint[k] = quad_points[i][k];
        }
        points[thePoint] = std::make_pair(quad_weights[i],i);
      }

      // Generate quadrature operator
      // This has to happen before we change the global index since
      // the point set likely will reorder the points
      qp2pce.reshape(npc,nqp);
      pce2qp.reshape(nqp,npc);
      typename point_set_type::iterator di = points.begin();
      typename point_set_type::iterator di_end = points.end();
      ordinal_type jdx = 0;
      for (; di != di_end; ++di) {
        ordinal_type j = di->second.second;
        for (ordinal_type i=0; i<npc; i++) {
          qp2pce(i,jdx) = quad_weights[j]*quad_vals[j][i] /
            basis.norm_squared(i);
          pce2qp(jdx,i) = quad_vals[j][i];
        }
        ++jdx;
      }

      // Generate linear ordering of points
      point_map.resize(nqp);
      ordinal_type idx=0;
      di = points.begin();
      while (di != di_end) {
        di->second.second = idx;
        point_map[idx] = di->first;
        ++idx;
        ++di;
      }


    }

    //! Destructor
    virtual ~QuadraturePseudoSpectralOperator() {}

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
    virtual void transformQP2PCE(
      const value_type& alpha,
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result,
      const value_type& beta,
      bool trans = false) const {
      ordinal_type ret;
      if (trans)
        ret = result.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha,
                              input, qp2pce, beta);
      else
        ret = result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha,
                              qp2pce, input, beta);
      TEUCHOS_ASSERT(ret == 0);
    }

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
      bool trans = false) const {
      if (trans) {
        TEUCHOS_ASSERT(input.numCols() <= pce2qp.numCols());
        TEUCHOS_ASSERT(result.numCols() == pce2qp.numRows());
        TEUCHOS_ASSERT(result.numRows() == input.numRows());
        blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, input.numRows(),
                  pce2qp.numRows(), input.numCols(), alpha, input.values(),
                  input.stride(), pce2qp.values(), pce2qp.stride(), beta,
                  result.values(), result.stride());
      }
      else {
        TEUCHOS_ASSERT(input.numRows() <= pce2qp.numCols());
        TEUCHOS_ASSERT(result.numRows() == pce2qp.numRows());
        TEUCHOS_ASSERT(result.numCols() == input.numCols());
        blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, pce2qp.numRows(),
                  input.numCols(), input.numRows(), alpha, pce2qp.values(),
                  pce2qp.stride(), input.values(), input.stride(), beta,
                  result.values(), result.stride());
      }
    }

  protected:

    //! Number of coefficients
    ordinal_type coeff_sz;

    //! Quadrature points
    point_set_type points;

    //! Map index to point term
    point_map_type point_map;

    //! Matrix mapping points to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> qp2pce;

    //! Matrix mapping coefficients to points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> pce2qp;

    //! BLAS wrappers
    Teuchos::BLAS<ordinal_type,value_type> blas;

  };

}

#endif
