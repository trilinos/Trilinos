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

#ifndef STOKHOS_RECURRENCEBASIS_HPP
#define STOKHOS_RECURRENCEBASIS_HPP

#include "Stokhos_OneDOrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Implementation of OneDOrthogPolyBasis based on the general
   * three-term recurrence relationship:
   * \f[
   *    \psi_{k+1}(x) = 
   *       \gamma_{k+1}\big( (\alpha_k - \delta_k x)\psi_k(x) - 
   *                         \beta_k\psi_{k-1}(x) \big)
   * \f]
   * for \f$k=0,\dots,P\f$ where \f$\psi_{-1}(x) = 0\f$, \f$\psi_{0}(x) = 1\f$,
   * and \f$\beta_{0} = 1\f$.  
   */
  /*!Derived classes implement the recurrence
   * relationship by implementing computeRecurrenceCoefficients().  If
   * \c normalize = \c true in the constructor, then 
   * \f$\gamma_k = 1/\sqrt{\langle\psi_k^2\rangle}\f$, otherwise 
   * \f$\gamma_k = 1\f$.
   * Note that a three term recurrence can always be defined with 
   * \f$\delta_k = 1\f$ in which case the polynomials are monic.  However
   * typical normalizations of some polynomial families (see 
   * Stokhos::LegendreBasis) require the extra terms.
   */
  template <typename ordinal_type, typename value_type>
  class RecurrenceBasis : 
    public OneDOrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Destructor
    virtual ~RecurrenceBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{

    //! Return order of basis (largest monomial degree \f$P\f$).
    virtual ordinal_type order() const;

    //! Return total size of basis (given by order() + 1).
    virtual ordinal_type size() const;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is given by order().
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     *
     * This method is implemented by computing \f$C_{ijk}\f$ using Gaussian
     * quadrature.
     */
    virtual Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const;

    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeSparseTripleProductTensor(ordinal_type order) const;

    //! Compute derivative double product tensor
    /*!
     * The \f$(i,j)\f$ entry of the tensor \f$B_{ij}\f$ is given by
     * \f$B_{ij} = \langle\psi_i'\psi_j\rangle\f$ where \f$\psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is the order of the basis.
     *
     * This method is implemented by computing \f$B_{ij}\f$ using Gaussian
     * quadrature.
     */
    virtual Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > computeDerivDoubleProductTensor() const;

    //! Evaluate each basis polynomial at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to order order().
     */
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const;

    /*! 
     * \brief Evaluate basis polynomial given by order \c order at given 
     * point \c point.
     */
    virtual value_type evaluate(const value_type& point, 
				ordinal_type order) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

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
     *
     * The quadrature points and weights are computed from the three-term
     * recurrence by solving a tri-diagional symmetric eigenvalue problem
     * (see Gene H. Golub and John H. Welsch, "Calculation of Gauss Quadrature 
     * Rules", Mathematics of Computation, Vol. 23, No. 106 (Apr., 1969), 
     * pp. 221-230).
     */
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    //! Get sparse grid rule number as defined by Dakota's \c webbur package
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.  A rule number of 10 is not defined by the webbur package, and
     * this rule number is used internally by Stokhos::SparseGridQuadrature
     * to pass an arbitrary one-dimensional basis to that package.
     */
    virtual int getSparseGridRule() const { return sparse_grid_rule; }

    //! Set sparse grid rule
    virtual void setSparseGridRule(int rule) { sparse_grid_rule = rule; }

    /*! 
     * \brief Get sparse grid rule growth rule as defined by 
     * Dakota's \c webbur package
     */
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.  Returns growth rule appropriate for Gaussian quadrature points.
     */
    virtual int getSparseGridGrowthRule() const { 
      return sparse_grid_growth_rule; };

    //! Set sparse grid growth rule
    virtual void setSparseGridGrowthRule(int rule) { 
      sparse_grid_growth_rule = rule; }

    //@}

    //! Return recurrence coefficients defined by above formula
    virtual void getRecurrenceCoefficients(Teuchos::Array<value_type>& alpha,
					   Teuchos::Array<value_type>& beta,
					   Teuchos::Array<value_type>& delta,
					   Teuchos::Array<value_type>& gamma) const;

    //! Evaluate basis polynomials and their derivatives at given point \c point
    virtual void evaluateBasesAndDerivatives(const value_type& point,
					     Teuchos::Array<value_type>& vals,
					     Teuchos::Array<value_type>& derivs) const;

    //! Set tolerance for zero in quad point generation
    virtual void setQuadZeroTol(value_type tol) {
      quad_zero_tol = tol; }

  protected:

    //! Constructor to be called by derived classes
    /*!
     * \c name is the name for the basis that will be displayed when
     * printing the basis, \c p is the order of the basis, \c normalize
     * indicates whether the basis polynomials should have unit-norm, and
     * \c quad_zero_tol is used to replace any quadrature point within this
     * tolerance with zero (which can help with duplicate removal in sparse
     * grid calculations).
     */
     RecurrenceBasis(const std::string& name, ordinal_type p, bool normalize);

     //! Copy constructor with specified order
     RecurrenceBasis(ordinal_type p, const RecurrenceBasis& basis);

    //! Compute recurrence coefficients
    /*!
     * Derived classes should implement this method to compute their
     * recurrence coefficients.  \c n is the number of coefficients to compute.
     * Derived classes should call this method in their constructor, followed
     * be setup() to fully setup the basis.  The method is also called by
     * getQuadPoints() if a quadrature order greater than twice the polynomial
     * order is requested.
     */
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const = 0;

    //! Setup basis after computing recurrence coefficients
    /*!
     * Derived classes should call this method after computing their recurrence
     * coefficients in their constructor to finish setting up the basis.
     */
    void setup();

  private:

    // Prohibit copying
    RecurrenceBasis(const RecurrenceBasis&);

    // Prohibit Assignment
    RecurrenceBasis& operator=(const RecurrenceBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Order of basis
    ordinal_type p;

    //! Normalize basis
    bool normalize;

    //! Tolerance for quadrature points near zero
    value_type quad_zero_tol;

    //! Sparse grid rule (as determined by Pecos)
    int sparse_grid_rule;

    //! Sparse grid growth rule (as determined by Pecos)
    int sparse_grid_growth_rule;

    //! Recurrence \f$\alpha\f$ coefficients
    Teuchos::Array<value_type> alpha;

    //! Recurrence \f$\beta\f$ coefficients
    Teuchos::Array<value_type> beta;

    //! Recurrence \f$\delta\f$ coefficients
    Teuchos::Array<value_type> delta;

    //! Recurrence \f$\gamma\f$ coefficients
    Teuchos::Array<value_type> gamma;

    //! Norms
    Teuchos::Array<value_type> norms;

  }; // class RecurrenceBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_RecurrenceBasisImp.hpp"

#endif
