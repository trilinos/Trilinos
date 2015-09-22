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

#ifndef STOKHOS_ANISOSPARSEGRIDQUADRATURE
#define STOKHOS_ANISOSPARSEGRIDQUADRATURE

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_DAKOTA

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_RCP.hpp"
#include "pecos_global_defs.hpp"

namespace Stokhos {

  /*!
   * \brief Defines quadrature for a tensor product basis by anisotropic
   * Smolyak sparse grids.
   */
  /*!
   * Requires Dakota webbur quadrature package, which is currently provided
   * through TriKota.  To enable, configure Stokhos with TriKota enabled
   * and see the TriKota instructions for building TriKota with Dakota.
   */
  template <typename ordinal_type, typename value_type>
  class AnisoSparseGridQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param product_basis product basis
     * \param sparse_grid_level sparse grid level defining the order of the 
     *                          quadrature.  If equal to 0, the level is
     *                          calculated using a heuristic formula.
     */
    AnisoSparseGridQuadrature(
      const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
      ordinal_type sparse_grid_level,
      value_type dim_weights[],
      value_type duplicate_tol = 1.0e-12,
      ordinal_type growth_rate = Pecos::MODERATE_RESTRICTED_GROWTH);

    //! Destructor
    virtual ~AnisoSparseGridQuadrature() {}

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

  private:

    // Prohibit copying
    AnisoSparseGridQuadrature(const AnisoSparseGridQuadrature&);

    // Prohibit Assignment
    AnisoSparseGridQuadrature& operator=(const AnisoSparseGridQuadrature& b);

    static void getMyPoints( int order, int dim, double x[] );
  
    static void getMyWeights( int order, int dim, double w[] );

  protected:

    //! Coordinate bases
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > coordinate_bases;

    //! Quadrature points
    Teuchos::Array< Teuchos::Array<value_type> > quad_points;

    //! Quadrature weights
    Teuchos::Array<value_type> quad_weights;

    //! Quadrature values
    Teuchos::Array< Teuchos::Array<value_type> >  quad_values;

    //! Static pointer for VPISparseGrid interface
    static AnisoSparseGridQuadrature *sgq;

  }; // class AnisoSparseGridQuadrature

} // namespace Stokhos
#include "Stokhos_AnisoSparseGridQuadratureImp.hpp"

#endif // HAVE_STOKHOS_DAKOTA

#endif // STOKHOS_ANISOSPARSEGRIDQUADRATURE
