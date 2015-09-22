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

#ifndef STOKHOS_DERIVBASIS_HPP
#define STOKHOS_DERIVBASIS_HPP

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {
  
  /*! 
   * \brief Abstract base class for multivariate orthogonal polynomials
   * that support computing double and triple products involving derivatives
   * of the basis polynomials.
   */
  template <typename ordinal_type, typename value_type>
  class DerivBasis : 
    public virtual OrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    DerivBasis() {};

    //! Destructor
    virtual ~DerivBasis() {};

    /*! 
     * \brief Compute triple product tensor 
     * \f$D_{ijk} = \langle\Psi_i\Psi_j D_v\Psi_k\rangle\f$ where 
     * \f$D_v\Psi_k\f$ represents the derivative of \f$\Psi_k\f$ in the 
     * direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeDerivTripleProductTensor(
      const Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
      const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const = 0;

    /*! 
     * \brief Compute double product tensor 
     * \f$B_{ij} = \langle \Psi_i D_v\Psi_j\rangle\f$ where \f$D_v\Psi_j\f$
     * represents the derivative of \f$\Psi_j\f$ in the direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > 
    computeDerivDoubleProductTensor() const = 0;

  private:

    // Prohibit copying
    DerivBasis(const DerivBasis&);

    // Prohibit Assignment
    DerivBasis& operator=(const DerivBasis& b);

  }; // class DerivBasis

} // Namespace Stokhos

#endif // STOKHOS_DERIVBASIS
