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

#ifndef STOKHOS_ONEDORTHOGPOLYBASIS_HPP
#define STOKHOS_ONEDORTHOGPOLYBASIS_HPP

#include <ostream>
#include <string>
#include "Stokhos_Dense3Tensor.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_DAKOTA
#include "pecos_global_defs.hpp"
#endif

//! Top-level namespace for Stokhos classes and functions.
namespace Stokhos {

  //! Abstract base class for 1-D orthogonal polynomials.
  /*!
   * This class provides an abstract interface for univariate orthogonal
   * polynomials.  Orthogonality is defined by the inner product
   * \f[
   *      (f,g) = \langle fg \rangle = 
   *              \int_{-\infty}^{\infty} f(x)g(x) \rho(x) dx
   * \f]
   * where \f$\rho\f$ is the density function of the measure associated with
   * the orthogonal polynomials.
   * See Stokhos::RecurrenceBasis for a general implementation
   * of this interface based on the three-term recurrence satisfied by
   * these polynomials.  Multivariate polynomials can be formed from
   * a collection of univariate polynomials through tensor products (see 
   * Stokhos::CompletePolynomialBasis).
   *
   * Like most classes in Stokhos, the class is templated on the ordinal
   * and value types.  Typically \c ordinal_type = \c int and \c value_type
   * = \c double.
   */
  template <typename ordinal_type, typename value_type>
  class OneDOrthogPolyBasis {
  public:

    //! Default constructor
    OneDOrthogPolyBasis() {};

    //! Destructor
    virtual ~OneDOrthogPolyBasis() {};

    //! Return order of basis (largest monomial degree \f$P\f$).
    virtual ordinal_type order() const = 0;

    //! Return total size of basis (given by order() + 1).
    virtual ordinal_type size() const = 0;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is given by order().
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const = 0;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const = 0;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual 
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const = 0;

    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeSparseTripleProductTensor(ordinal_type order) const = 0;

    //! Compute derivative double product tensor
    /*!
     * The \f$(i,j)\f$ entry of the tensor \f$B_{ij}\f$ is given by
     * \f$B_{ij} = \langle\psi_i'\psi_j\rangle\f$ where \f$\psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is the order of the basis.
     */
    virtual 
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > 
    computeDerivDoubleProductTensor() const = 0;

    //! Evaluate each basis polynomial at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to order order().
     */
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const = 0;

    /*! 
     * \brief Evaluate basis polynomial given by order \c order at given 
     * point \c point.
     */
    virtual value_type evaluate(const value_type& point, 
				ordinal_type order) const = 0;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const {};

    //! Return string name of basis
    virtual const std::string& getName() const = 0;

    /*! 
     * \brief Compute quadrature points, weights, and values of 
     * basis polynomials at given set of points \c points.
     */
    /*!
     * \c quad_order specifies the order to which the quadrature should be
     * accurate, not the number of quadrature points.  The number of points
     * is given by (\c quad_order + 1) / 2.   Note however the passed arrays
     * do NOT need to be sized correctly on input as they will be resized
     * appropriately.
     */
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const = 0;

    /*! 
     * \brief Clone this object with the option of building a higher order
     * basis.
     */
    /*!
     * This method is following the Prototype pattern (see Design Pattern's textbook).
     * The slight variation is that it allows the order of the polynomial to be modified,
     * otherwise an exact copy is formed. The use case for this is creating basis functions
     * for column indices in a spatially varying adaptive refinement context.
     */
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const = 0;

    //! Get sparse grid rule number as defined by Dakota's \c webbur package
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.  The current rule definitions are:
     * 1  Clenshaw-Curtis
     * 2  Fejer Type 2
     * 3  Gauss-Patterson
     * 4  Gauss-Legendre
     * 5  Gauss-Hermite
     * 6  Generalized Gauss-Hermite
     * 7  Gauss-Laguerre
     * 8  Generalized Gauss-Laguerre
     * 9  Gauss-Jacobi
     * 10 Golub-Welsch (Gauss points for arbitrary weight function)
     * 11 Genz-Keister (Gauss-Patterson-type Hermite)
     * 12 Newton-Cotes
     */
    virtual int getSparseGridRule() const = 0;

    //! Set sparse grid rule
    virtual void setSparseGridRule(int rule) = 0;

    /*! 
     * \brief Get sparse grid rule growth rule as defined by 
     * Dakota's \c webbur package
     */
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.  The current rule definitions are:
     * 1  Default growth
     * 2  Slow linear
     * 3  Slow linear odd
     * 4  Moderate linear
     * 5  Slow exponential
     * 6  Moderate exponential
     * 7  Full exponential
     */
    virtual int getSparseGridGrowthRule() const = 0;

    //! Set sparse grid growth rule
    virtual void setSparseGridGrowthRule(int rule) = 0;

  private:

    // Prohibit copying
    OneDOrthogPolyBasis(const OneDOrthogPolyBasis&);

    // Prohibit Assignment
    OneDOrthogPolyBasis& operator=(const OneDOrthogPolyBasis& b);
  

  }; // class OrthogPolyBasis

  //! Print basis to stream \c os.
  template <typename ordinal_type, typename value_type> 
  std::ostream& 
  operator << (std::ostream& os, 
	       const OneDOrthogPolyBasis<ordinal_type, value_type>& b) {
    b.print(os);
    return os;
  }

} // Namespace Stokhos

#endif
