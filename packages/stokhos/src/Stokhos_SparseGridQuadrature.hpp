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

#ifndef STOKHOS_SPARSEGRIDQUADRATURE
#define STOKHOS_SPARSEGRIDQUADRATURE

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_DAKOTA

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_RCP.hpp"
#include "pecos_global_defs.hpp"

namespace Stokhos {

  /*!
   * \brief Defines quadrature for a tensor product basis by Smolyak sparse
   * grids.
   */
  /*!
   * Requires Dakota webbur quadrature package, which is currently provided
   * through TriKota.  To enable, configure Stokhos with TriKota enabled
   * and see the TriKota instructions for building TriKota with Dakota.
   */
  template <typename ordinal_type, typename value_type>
  class SparseGridQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param product_basis product basis
     * \param sparse_grid_level sparse grid level defining the order of the 
     *                          quadrature.  If equal to 0, the level is
     *                          calculated using a heuristic formula.
     */
    SparseGridQuadrature(
      const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis, 
      ordinal_type sparse_grid_level = 0,
      value_type duplicate_tol = 1.0e-12,
      ordinal_type growth_rate = Pecos::SLOW_RESTRICTED_GROWTH);

    //! Destructor
    virtual ~SparseGridQuadrature() {}

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
    SparseGridQuadrature(const SparseGridQuadrature&);

    // Prohibit Assignment
    SparseGridQuadrature& operator=(const SparseGridQuadrature& b);

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
    static SparseGridQuadrature *sgq;

  }; // class SparseGridQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_SparseGridQuadratureImp.hpp"

#endif // HAVE_STOKHOS_DAKOTA

#endif // STOKHOS_SPARSEGRIDQUADRATURE
