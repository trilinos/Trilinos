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

#ifndef STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP
#define STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_SmolyakPseudoSpectralOperator.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  /*!
   * \brief Defines quadrature for a tensor product basis by Smolyak sparse
   * grids.
   */
  /*!
   * This class generates the sparse grids using the
   * SmolyakPseudoSpectralOperator and doesn't rely on Dakota.
   */
  template < typename ordinal_type, typename value_type,
             typename point_compare_type =
             typename DefaultPointCompare<ordinal_type,value_type>::type >
  class SmolyakSparseGridQuadrature :
    public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param product_basis product basis
     * \param index_set     index set defining growth levels
     * \param point_compare comparison functor used in ordering points
     */
    template <typename index_set_type>
    SmolyakSparseGridQuadrature(
      const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
      const index_set_type& index_set,
      const value_type duplicate_tol = 1.0e-12,
      const point_compare_type& point_compare = point_compare_type());

    //! Destructor
    virtual ~SmolyakSparseGridQuadrature() {}

    //! Get number of quadrature points
    virtual ordinal_type size() const { return quad_weights.size(); }

    //! Get quadrature points
    /*!
     * Array is dimensioned Q-by-d where Q is the number of quadrature
     * points and d is the dimension of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> >&
    getQuadPoints() const;

    //! Get quadrature weights
    /*!
     * Array is of size Q where Q is the number of quadrature points.
     */
    virtual const Teuchos::Array<value_type>&
    getQuadWeights() const;

    //! Get values of basis at quadrature points
    /*!
     * Array is dimensioned Q-by-P where Q is the number of quadrature
     * points and P is the size of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> > &
    getBasisAtQuadPoints() const;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const;

  private:

    // Prohibit copying
    SmolyakSparseGridQuadrature(const SmolyakSparseGridQuadrature&);

    // Prohibit Assignment
    SmolyakSparseGridQuadrature& operator=(const SmolyakSparseGridQuadrature& b);

  protected:

    //! Quadrature points
    Teuchos::Array< Teuchos::Array<value_type> > quad_points;

    //! Quadrature weights
    Teuchos::Array<value_type> quad_weights;

    //! Quadrature values
    Teuchos::Array< Teuchos::Array<value_type> >  quad_values;

  }; // class SmolyakSparseGridQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_SmolyakSparseGridQuadratureImp.hpp"

#endif //STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP
