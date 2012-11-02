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

#ifndef STOKHOS_QUADRATURE_PSEUDO_SPECTRAL_OPERATOR_HPP
#define STOKHOS_QUADRATURE_PSEUDO_SPECTRAL_OPERATOR_HPP

#include "Stokhos_PseudoSpectralOperator.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

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
      point_type point(dim);
      for (ordinal_type i=0; i<nqp; ++i) {
	for (ordinal_type k=0; k<dim; k++) {
	  point[k] = quad_points[i][k];
	}
	points[point] = std::make_pair(quad_weights[i],i);
      }

      // Generate quadrature operator
      // This has to happen before we change the global index since
      // the point set likely will reorder the points
      A.reshape(npc,nqp);
      typename point_set_type::iterator di = points.begin();
      typename point_set_type::iterator di_end = points.end();
      ordinal_type jdx = 0;
      for (; di != di_end; ++di) {
	ordinal_type j = di->second.second;
	for (ordinal_type i=0; i<npc; i++) 
	  A(i,jdx) = quad_weights[j]*quad_vals[j][i] / basis.norm_squared(i);
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

    //! Apply pseudo-spectral quadrature operator
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void apply(const value_type& alpha, 
	       const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
	       Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
	       const value_type& beta,
	       bool trans = false) const {
      ordinal_type ret;
      if (trans)
	ret = result.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha, input,
			      A, beta);
      else
	ret = result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A,
			      input, beta);
      TEUCHOS_ASSERT(ret == 0);
    }

  protected:

    //! Number of coefficients
    ordinal_type coeff_sz;

    //! Quadrature points
    point_set_type points;

    //! Map index to point term
    point_map_type point_map;

    //! Matrix mapping point to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> A;

  };

}

#endif
