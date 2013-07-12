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

#ifndef STOKHOS_MONOPROJPCEBASIS_HPP
#define STOKHOS_MONOPROJPCEBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "Stokhos_RecurrenceBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {

  /*! 
   * \brief Generates three-term recurrence using the Lanczos 
   * procedure applied to a polynomial chaos expansion in another basis.
   */
  template <typename ordinal_type, typename value_type>
  class MonoProjPCEBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansion defining new density function
     * \param quad quadrature data for basis of PC expansion
     */
    MonoProjPCEBasis(
      ordinal_type p,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
      const Stokhos::Quadrature<ordinal_type, value_type>& quad,
      const Stokhos::Sparse3Tensor<ordinal_type, value_type>& Cijk,
      bool limit_integration_order = false);

    //! Destructor
    ~MonoProjPCEBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

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
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const;

    //@}

    //! Get new coefficients in this new basis
    value_type getNewCoeffs(ordinal_type i) const;

    //! Map expansion coefficients from this basis to original
    void transformCoeffs(const value_type *in, value_type *out) const;

  protected:

    //! \name Implementation of Stokhos::RecurrenceBasis methods
    //@{ 

    //! Compute recurrence coefficients
    virtual bool
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta,
				  Teuchos::Array<value_type>& gamma) const;

    //@}

  private:

    // Prohibit copying
    MonoProjPCEBasis(const MonoProjPCEBasis&);

    // Prohibit Assignment
    MonoProjPCEBasis& operator=(const MonoProjPCEBasis& b);

    //! Copy constructor with specified order
    MonoProjPCEBasis(ordinal_type p, const MonoProjPCEBasis& basis);

  protected:

    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> matrix_type;
    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> vector_type;

    //! Flag indicating whether to limit the integration order
    bool limit_integration_order;

    //! Size of PC expansion
    ordinal_type pce_sz;

    //! Basis norms
    Teuchos::Array<value_type> pce_norms;

    //! Stores full set of alpha coefficients
    Teuchos::Array<value_type> a;

    //! Stores full set of beta coefficients
    Teuchos::Array<value_type> b;

    //! Basis vectors
    matrix_type basis_vecs;

    //! Projection of pce in new basis
    vector_type new_pce;

  }; // class MonoProjPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_MonoProjPCEBasisImp.hpp"

#endif
