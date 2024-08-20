// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TENSOR_PRODUCT_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_TENSOR_PRODUCT_PSEUDO_SPECTRAL_OPERATOR_HPP

#include "Stokhos_PseudoSpectralOperator.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_TensorProductBasis.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

namespace Stokhos {

  /*!
   * \brief An operator for building pseudo-spectral coefficients using 
   * tensor-product quadrature.
   */
  template <typename ordinal_t, 
	    typename value_t, 
	    typename point_compare_type = 
	    typename DefaultPointCompare<ordinal_t,value_t>::type>
  class TensorProductPseudoSpectralOperator : 
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
    TensorProductPseudoSpectralOperator(
      const ProductBasis<ordinal_type,value_type>& product_basis, 
      bool use_pst_ = false, 
      multiindex_type multiindex = multiindex_type(),
      const point_compare_type& point_compare = point_compare_type()) : 
      use_pst(use_pst_),
      coeff_sz(product_basis.size()),
      dim(product_basis.dimension()),
      points(point_compare),
      reorder(false) {

      // Check if we need to reorder for PST
      // This isn't a perfect test since it doesn't catch not tensor-product
      // bases
      typedef LexographicLess< MultiIndex<ordinal_type> > lexo_type;
      if (use_pst) {
	typedef TensorProductBasis<ordinal_type,value_type,lexo_type> tb_type;
	const tb_type *tb = dynamic_cast<const tb_type*>(&product_basis);
	if (tb == NULL) 
	  reorder = true;
      }

      Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > bases_ = product_basis.getCoordinateBases();

      // Make sure order of each 1-D basis is large enough for 
      // supplied set of coefficients
      Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(bases_);
      multiindex_type max_orders = product_basis.getMaxOrders();
      for (ordinal_type i=0; i<dim; ++i) {
	if (bases[i]->order() < max_orders[i])
	  bases[i] = bases[i]->cloneWithOrder(max_orders[i]);
      }

      // Check for default quadrature order
      if (multiindex.dimension() == 0)
	multiindex = max_orders;

      // Compute quad points, weights, values
      Teuchos::Array< Teuchos::Array<value_type> > gp(dim);
      Teuchos::Array< Teuchos::Array<value_type> > gw(dim);
      Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > gv(dim);
      multiindex_type n(dim);
      if (use_pst) {
	qp2pce_k.resize(dim);
	pce2qp_k.resize(dim);
      }
      for (ordinal_type k=0; k<dim; k++) {
	bases[k]->getQuadPoints(2*multiindex[k], gp[k], gw[k], gv[k]);
	n[k] = gp[k].size()-1;

	// Generate quadrature operators for PST
	if (use_pst) {
	  ordinal_type npc = bases[k]->size();
	  ordinal_type nqp = gp[k].size();
	  qp2pce_k[k].reshape(npc,nqp);
	  pce2qp_k[k].reshape(nqp,npc);
	  qp2pce_k[k].putScalar(1.0);
	  pce2qp_k[k].putScalar(1.0);
	  for (ordinal_type j=0; j<nqp; j++) {
	    for (ordinal_type i=0; i<npc; i++) {
	      qp2pce_k[k](i,j) *= gw[k][j]*gv[k][j][i] / 
		bases[k]->norm_squared(i);
	      pce2qp_k[k](j,i) *= gv[k][j][i];
	    }
	  }
	}
      }

      // Generate points and weights
      Stokhos::TensorProductIndexSet<ordinal_type> pointIndexSet(n);
      typedef typename TensorProductIndexSet<ordinal_type>::iterator index_iterator;
      index_iterator point_iterator = pointIndexSet.begin();
      index_iterator point_end = pointIndexSet.end();
      point_type point(dim);
      while (point_iterator != point_end) {
	value_type w = 1.0;
	for (ordinal_type k=0; k<dim; k++) {
	  point[k] = gp[k][(*point_iterator)[k]];
	  w *= gw[k][(*point_iterator)[k]];
	}
	points[point] = std::make_pair(w,ordinal_type(0));
	++point_iterator;
      }

      // Generate linear ordering of points
      ordinal_type nqp = points.size();
      point_map.resize(nqp);      
      ordinal_type idx=0;
      typename point_set_type::iterator di = points.begin();
      typename point_set_type::iterator di_end = points.end();
      while (di != di_end) {
	di->second.second = idx;
	point_map[idx] = di->first;
	++idx;
	++di;
      }

      // Build permutation array for reordering if coefficients aren't
      // lexographically ordered
      if (use_pst && reorder) {
	typedef std::map<multiindex_type,ordinal_type,lexo_type> reorder_type;
	reorder_type reorder_map;
	for (ordinal_type i=0; i<coeff_sz; ++i)
	  reorder_map[product_basis.term(i)] = i;
	perm.resize(coeff_sz);
	ordinal_type idx = 0;
	for (typename reorder_type::iterator it = reorder_map.begin();
	     it != reorder_map.end(); ++it)
	  perm[idx++] = it->second;
	TEUCHOS_ASSERT(idx == coeff_sz);
      }

      // Generate quadrature operator for direct apply
      // Always do this because we may not use PST in some situations
      ordinal_type npc = product_basis.size();
      qp2pce.reshape(npc,nqp);
      pce2qp.reshape(nqp,npc);
      qp2pce.putScalar(1.0);
      pce2qp.putScalar(1.0);
      for (point_iterator = pointIndexSet.begin(); 
	   point_iterator != point_end; 
	   ++point_iterator) {
	for (ordinal_type k=0; k<dim; k++)
	  point[k] = gp[k][(*point_iterator)[k]];
	ordinal_type j = points[point].second;
	for (ordinal_type i=0; i<npc; i++) {
	  for (ordinal_type k=0; k<dim; k++) {
	    ordinal_type l = (*point_iterator)[k];
	    ordinal_type m = product_basis.term(i)[k];
	    qp2pce(i,j) *= gw[k][l]*gv[k][l][m] / bases[k]->norm_squared(m);
	    pce2qp(j,i) *= gv[k][l][m];
	  }
	}
      }
      
    }

    //! Destructor
    virtual ~TensorProductPseudoSpectralOperator() {}

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
      if (use_pst)
	apply_pst(qp2pce_k, alpha, input, result, beta, trans, false, reorder);
      else
	apply_direct(qp2pce, alpha, input, result, beta, trans);
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

      // Only use PST if it was requested and number of rows/columns match
      // (since we may allow fewer rows in input than columns of pce2qp)
      bool can_use_pst = use_pst;
      if (trans) 
	can_use_pst = can_use_pst && (input.numCols() == pce2qp.numRows());
      else
	can_use_pst = can_use_pst && (input.numRows() == pce2qp.numRows());

      if (can_use_pst)	
	apply_pst(pce2qp_k, alpha, input, result, beta, trans, reorder, false);
      else 
	apply_direct(pce2qp, alpha, input, result, beta, trans);
    }

  protected:

    //! Apply transformation operator using direct method
    void apply_direct(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A,
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans) const {
      if (trans) {
	TEUCHOS_ASSERT(input.numCols() <= A.numCols());
	TEUCHOS_ASSERT(result.numCols() == A.numRows());
	TEUCHOS_ASSERT(result.numRows() == input.numRows());
	blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, input.numRows(), 
		  A.numRows(), input.numCols(), alpha, input.values(), 
		  input.stride(), A.values(), A.stride(), beta, 
		  result.values(), result.stride());
      }
      else {
	TEUCHOS_ASSERT(input.numRows() <= A.numCols());
	TEUCHOS_ASSERT(result.numRows() == A.numRows());
	TEUCHOS_ASSERT(result.numCols() == input.numCols());
	blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, A.numRows(), 
		  input.numCols(), input.numRows(), alpha, A.values(), 
		  A.stride(), input.values(), input.stride(), beta, 
		  result.values(), result.stride());
      }
    }

    //! Apply tranformation operator using PST method
    /*!
     * For k=1,...,d, let A_k be the m_k-by-n_k quadrature operator for
     * dimension k and A = A_1 \otimes ... \otimes A_d be m-\by-n where 
     * m = m_1...m_d and n = n_1...n_d.  For any two matrices
     * B and C and vector x, (B \otimes C)x = vec( C*X*B^T ) = 
     * vec( C*(B*X^T)^T ) where x = vec(X), and the vec() operator makes a 
     * vector out of a matrix by stacking columns.  Applying this formula
     * recursively to A yields the simple algorithm for computing y = A*x
     * (x is n-by-1 and y is m-by-1):
     *    X = x;
     *    for k=1:d
     *      n = n / n_k;
     *      X = reshape(X, n, n_k);
     *      X = A_k*X';
     *      n = n * m_k;
     *    end
     *    y = reshape(X, m, 1);
     *
     * When x has p columns, it is somehwat more complicated because the
     * standard transpose above isn't correct.  Instead the transpose
     * needs to be applied to each p block in X.
     *
     * The strategy here for dealing with transposed input and result
     * is to transpose them when copying to/from the temporary buffers
     * used in the algorithm.
     */
    void apply_pst(
      const Teuchos::Array< Teuchos::SerialDenseMatrix<ordinal_type,value_type> >& Ak,
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans,
      bool reorder_input,
      bool reorder_result) const {

      ordinal_type n, m, p;
      if (trans) {
	TEUCHOS_ASSERT(input.numRows() == result.numRows());
	n = input.numCols();
	p = input.numRows();
	m = result.numCols();
      }
      else {
	TEUCHOS_ASSERT(input.numCols() == result.numCols());
	n = input.numRows();
	p = input.numCols();
	m = result.numRows();
      }    
      ordinal_type M = 1;
      ordinal_type N = 1;
      for (ordinal_type k=0; k<dim; ++k) {
	M *= Ak[k].numRows();
	N *= Ak[k].numCols();
      }
      TEUCHOS_ASSERT(n == N);
      TEUCHOS_ASSERT(m == M);
      ordinal_type sz = std::max(n,m);  

      // Temporary buffers used in algorithm
      Teuchos::Array<value_type> tmp1(sz*p),  tmp2(sz*p);

      // Copy input to tmp1 (transpose if necessary)
      if (trans) {
	if (reorder_input) {
	  ordinal_type idx;
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<n; i++) {
	      idx = perm[i];
	      tmp1[i+j*n] = input(j,idx);
	    }
	}
	else {
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<n; i++)
	      tmp1[i+j*n] = input(j,i);
	}
      }
      else {
	if (reorder_input) {
	  ordinal_type idx;
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<n; i++) {
	      idx = perm[i];
	      tmp1[i+j*n] = input(idx,j);
	    }
	}
	else {
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<n; i++)
	      tmp1[i+j*n] = input(i,j);
	}
      }
      
      // Loop over each term in Kronecker product
      for (ordinal_type k=0; k<dim; k++) {
	ordinal_type mk = Ak[k].numRows();
	ordinal_type nk = Ak[k].numCols();
	n = n / nk;
	
	// x = reshape(x, n, n_k)
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(
	  Teuchos::View, &tmp1[0], n, n, nk*p);

	// y = x'
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> y(
	  Teuchos::View, &tmp2[0], nk, nk, n*p);
	for (ordinal_type l=0; l<p; l++)
	  for (ordinal_type j=0; j<n; j++)
	    for (ordinal_type i=0; i<nk; i++)
	      y(i,j+l*n) = x(j,i+l*nk);

	// x = A_k*x
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> z(
	  Teuchos::View, &tmp1[0], mk, mk, n*p);
	ordinal_type ret = 
	  z.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ak[k], y, 0.0);
	TEUCHOS_ASSERT(ret == 0);

	n = n * mk;
      }

      // Sum answer into result (transposing and/or reordering if necessary)
      if (trans) {
	if (reorder_result) {
	  ordinal_type idx;
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<m; i++) {
	      idx = perm[i];
	      result(j,idx) = beta*result(j,idx) + alpha*tmp1[i+j*m];
	    }
	}
	else {
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<m; i++)
	      result(j,i) = beta*result(j,i) + alpha*tmp1[i+j*m];
	}
      }
      else {
	if (reorder_result) {
	  ordinal_type idx;
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<m; i++) {
	      idx = perm[i];
	      result(idx,j) = beta*result(idx,j) + alpha*tmp1[i+j*m];
	    }
	}
	else {
	  for (ordinal_type j=0; j<p; j++)
	    for (ordinal_type i=0; i<m; i++)
	      result(i,j) = beta*result(i,j) + alpha*tmp1[i+j*m];
	}
      }
    }

  protected:

    //! Use partial-summation-transformation
    bool use_pst;

    //! Number of coefficients
    ordinal_type coeff_sz;

    //! Dimension
    ordinal_type dim;

    //! Quadrature points
    point_set_type points;

    //! Map index to point term
    point_map_type point_map;

    //! Do we need to reorder coefficients for PST
    bool reorder;

    //! Permutation array when reordering for PST
    Teuchos::Array<ordinal_type> perm;

    //! Matrix mapping points to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> qp2pce;

    //! Matrix mapping coefficients to points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> pce2qp;

    //! Matrix mapping points to coefficients for each dimension for PST
    Teuchos::Array< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > qp2pce_k;

    //! Matrix mapping coefficients to points for each dimension for PST
    Teuchos::Array< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > pce2qp_k;

    //! BLAS wrappers
    Teuchos::BLAS<ordinal_type,value_type> blas;

  };

}

#endif
