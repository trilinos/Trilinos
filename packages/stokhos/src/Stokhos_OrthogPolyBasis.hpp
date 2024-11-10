// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_ORTHOGPOLYBASIS_HPP
#define STOKHOS_ORTHOGPOLYBASIS_HPP

#include <ostream>
#include <string>
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_Dense3Tensor.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  //! Abstract base class for multivariate orthogonal polynomials.
  /*!
   * This class provides an abstract interface for multivariate orthogonal
   * polynomials.  Orthogonality is defined by the inner product
   * \f[
   *      (f,g) = \langle fg \rangle = 
   *              \int_{R^d} f(x)g(x) \rho(x) dx
   * \f]
   * where \f$\rho\f$ is the density function of the measure associated with
   * the orthogonal polynomials and \f$d\f$ is the dimension of the domain.
   *
   * Like most classes in Stokhos, the class is templated on the ordinal
   * and value types.  Typically \c ordinal_type = \c int and \c value_type
   * = \c double.
   */
  template <typename ordinal_type, typename value_type>
  class OrthogPolyBasis {
  public:

    //! Constructor
    OrthogPolyBasis() {};

    //! Destructor
    virtual ~OrthogPolyBasis() {};

    //! Return order of basis
    virtual ordinal_type order() const = 0;

    //! Return dimension of basis
    virtual ordinal_type dimension() const = 0;

    //! Return total size of basis
    virtual ordinal_type size() const = 0;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\Psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is size()-1.
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const = 0;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const = 0;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j,k=0,\dots,P\f$ where
     * \f$P\f$ is size()-1.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const = 0;

    //! Compute linear triple product tensor where k = 0,1
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeLinearTripleProductTensor() const = 0;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const = 0;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(
      const Teuchos::ArrayView<const value_type>& point,
      Teuchos::Array<value_type>& basis_vals) const = 0;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const = 0;

    //! Return string name of basis
    virtual const std::string& getName() const = 0;

  private:

    // Prohibit copying
    OrthogPolyBasis(const OrthogPolyBasis&);

    // Prohibit Assignment
    OrthogPolyBasis& operator=(const OrthogPolyBasis& b);

  }; // class OrthogPolyBasis

  //! Print basis to stream \c os.
  template <typename ordinal_type, typename value_type> 
  std::ostream& operator << (std::ostream& os, 
			     const OrthogPolyBasis<ordinal_type, value_type>& b)
  {
    b.print(os);
    return os;
  }

} // Namespace Stokhos

#endif // STOKHOS_ORTHOGPOLYBASIS
