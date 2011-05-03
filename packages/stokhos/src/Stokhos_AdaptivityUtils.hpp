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

#ifndef STOKHOS_AdaptivityUtils_HPP
#define STOKHOS_AdaptivityUtils_HPP

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_ProductBasis.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include <vector>

namespace Stokhos {
namespace adapt_utils {

  /** Builds and returns an adapted graph given a set of row basis functions.
    * This does all required global communication to construct this graph.
    */
  Teuchos::RCP<Epetra_CrsGraph> buildAdaptedGraph(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        bool onlyUseLinear = false,
        int kExpOrder = -1);

  /** Builds and returns an adapted graph given a set of row basis functions.
    * This does all required global communication to construct this graph.
    */
  Teuchos::RCP<Epetra_CrsGraph> buildAdaptedGraph(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<int> & myRowGidOffsets,std::vector<int> & myColGidOffsets,
        bool onlyUseLinear = false,
        int kExpOrder=-1);

   /** Construct a row map that is appropriate for the set of adapted
     * basis functions.
     *
     * \param[in] determGraph Graph for the deterministic problem
     * \param[in] per_dof_row_basis Stochastic basis functions for each deterministic degree of freedom
     * \param[out] myRowGidOffsets Will be of length <code>per_dof_row_basis.size()</code> on exit. All data
     *                          will be overwritten. This contains the starting GID 
     *                          of each deterministic degree of freedom.
     *
     * \returns The adapted row map.
     */
   Teuchos::RCP<Epetra_Map> buildAdaptedRowMapAndOffsets(
        const Epetra_Comm & Comm,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<int> & myRowGidOffsets);

   /** Construct a row map that is appropriate for the set of adapted
     * basis functions.
     */
   Teuchos::RCP<Epetra_Map> buildAdaptedRowMap(
        const Epetra_Comm & Comm,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis);

   /** Build offsets mapping a local column id to a GID. Note that this function requires parallel
     * communication (for parallel code).
     *
     * \param[in] determGraph Graph for the deterministic problem
     * \param[in] myRowGidOffsets Computed by <code>buildAdaptedRowMapAndColOffsets</code>
     * \param[out] myRowGidOffsets Will be of length <code>determGraph.ColMap().NumMyElements()</code> on exit. All data
     *                             will be overwritten. This contains the starting GID of each deterministic degree of
     *                             freedom in the column map.
     */
   void buildAdaptedColOffsets(
           const Epetra_CrsGraph & determGraph,
           const std::vector<int> & myRowGidOffsets,
           std::vector<int> & myColGidOffsets);

   /** Use the deterministic graph and the basis function on this processor to determine the
     * basis functions for all the column entries on this processor. Note that this function requires
     * parallel communication (for parallel code).
     *
     * \param[in] determGraph Graph for the deterministic problem
     * \param[in] per_dof_row_basis Stochastic basis functions for each deterministic degree of freedom
     * \param[out] per_dof_col_basis Stochastic basis functions for each degree of freedom
     *                               in the deterministic column map.
     */
   void buildColBasisFunctions(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_col_basis);
}
    
} // namespace Stokhos

#endif // STOKHOS_AdaptivityUtils_HPP
