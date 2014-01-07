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
