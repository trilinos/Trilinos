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

#ifndef STOKHOS_CLENSHAWCURTISLEGENDREBASIS_HPP
#define STOKHOS_CLENSHAWCURTISLEGENDREBASIS_HPP

#include "Stokhos_LegendreBasis.hpp"

namespace Stokhos {

  //! Legendre polynomial basis using Clenshaw-Curtis quadrature points
  /*!
   * This is the same as Stokhos::LegendreBasis, but uses Clenshaw-Curtis
   * quadrature points (instead of Gauss-Legendre).
   */
  template <typename ordinal_type, typename value_type>
  class ClenshawCurtisLegendreBasis : 
    public LegendreBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param normalize whether polynomials should be given unit norm
     */
    ClenshawCurtisLegendreBasis(ordinal_type p, bool normalize = false,
				bool isotropic = false);

    //! Destructor
    ~ClenshawCurtisLegendreBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{ 

    /*! 
     * \brief Compute quadrature points, weights, and values of 
     * basis polynomials at given set of points \c points.
     */
    /*!
     * \c quad_order specifies the order to which the quadrature should be
     * accurate, not the number of quadrature points which for Clenshaw-Curtis
     * points is equal to the number of points
     */
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    /*!
     * Return polynomial degree of exactness for a given number of quadrature
     * points.
     */
    virtual ordinal_type quadDegreeOfExactness(ordinal_type n) const;

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

    //! Evaluate coefficient growth rule for Smolyak-type bases
    virtual ordinal_type coefficientGrowth(ordinal_type n) const;

    //! Evaluate point growth rule for Smolyak-type bases
    virtual ordinal_type pointGrowth(ordinal_type n) const;

    //@}

  protected:

    //! Copy constructor with specified order
    ClenshawCurtisLegendreBasis(ordinal_type p, 
				const ClenshawCurtisLegendreBasis& basis);

  private:

    // Prohibit copying
    ClenshawCurtisLegendreBasis(const ClenshawCurtisLegendreBasis&);

    // Prohibit Assignment
    ClenshawCurtisLegendreBasis& operator=(const ClenshawCurtisLegendreBasis& b);

  protected:

    //! Flag determining if expansion is iostropic (same basis in every dim)
    bool isotropic;

  }; // class ClenshawCurtisLegendreBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ClenshawCurtisLegendreBasisImp.hpp"

#endif
