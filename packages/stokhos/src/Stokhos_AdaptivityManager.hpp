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

#ifndef STOKHOS_AdaptivityManager_HPP
#define STOKHOS_AdaptivityManager_HPP

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_SGOperator.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include <vector>

namespace Stokhos {

   /** Describes and constructs all things needed for adaptivity.
     */
   class AdaptivityManager {
   public:
      AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_CrsGraph & determ_graph, 
         bool onlyUseLinear,int kExpOrder,
         bool scaleOp=true);

      AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_Comm & comm,
         bool scaleOp=true);

      /** Given a deterministic local column ID and a basis index
        * determine the global column ID in the fully assembled system.
        *
        * \param[in] determLid Deterministic local column ID
        * \param[in] basisIndex Index into the stochastic basis
        *                       associated with this column ID.
        *
        * \returns Associated index in the global system.
        *
        * \note basisIndex must be less than getRowStochasticBasisSize(determLid)
        */
      inline int getGlobalColId(int determLid,int basisIndex) const
      { return myColGidOffsets_[determLid]+basisIndex; }

      /** Given a deterministic local row ID and a basis index
        * determine the global row ID in the fully assembled system.
        *
        * \param[in] determLid Deterministic local row ID
        * \param[in] basisIndex Index into the stochastic basis
        *                       associated with this row ID.
        *
        * \returns Associated index in the global system.
        *
        * \note basisIndex must be less than getRowStochasticBasisSize(determLid)
        */
      inline int getGlobalRowId(int determLid,int basisIndex) const
      { return myRowGidOffsets_[determLid]+basisIndex; }

      /** Build a CRS matrix from the internally constructed graph
        */
      Teuchos::RCP<Epetra_CrsMatrix> buildMatrixFromGraph() const;

      /** Build a CRS graph from a determinstic graph
        */
      void setupWithGraph(const Epetra_CrsGraph & graph,bool onlyUseLinear,int kExpOrder);

      /** Get map associated with this set of adaptive indices
        */
      Teuchos::RCP<const Epetra_Map> getAdaptedMap() const
      { return rowMap_; }

      /** Setup operator
        */
      void setupOperator(Epetra_CrsMatrix & A,const Sparse3Tensor<int,double> & Cijk,Stokhos::VectorOrthogPoly<Epetra_Operator> & poly,
                         bool onlyUseLinear=false,bool includeMean=true) const;

      /** Sum into a matrix constructed from <code>buildMatrixFromGraph</code>
        * using the Cjik tensor a matrix J_k
        */ 
      void sumInOperator(Epetra_CrsMatrix & A,const Stokhos::Sparse3Tensor<int,double> & Cijk,int k,const Epetra_CrsMatrix & J_k) const;

      /** Copy to an adaptive vector from a set of blocked vectors
        */
      void copyToAdaptiveVector(const Stokhos::EpetraVectorOrthogPoly & x_sg,Epetra_Vector & x) const;

      /** Copy from an adaptive vector to a set of blocked vectors
        */
      void copyFromAdaptiveVector(const Epetra_Vector & x,Stokhos::EpetraVectorOrthogPoly & x_sg) const;

      /** How many stochastic degrees of freedom are associated
        * with a particular deterministic row degree of freedom.
        */
      int getRowStochasticBasisSize(int determLid) const
      { return sg_basis_row_dof_[determLid]->size(); }
  
      /** How many stochastic degrees of freedom are associated
        * with a particular deterministic row degree of freedom.
        */
      int getColStochasticBasisSize(int determLid) const
      { return sg_basis_col_dof_[determLid]->size(); }
  
      /** Get master stochastic basis
        */
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > getMasterStochasticBasis() const
      { return sg_master_basis_; } 
  
      /** Get stochastic basis associated with a particular deterministic row local id.
        */
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > getRowStochasticBasis(int determLid) const
      { return sg_basis_row_dof_[determLid]; } 
  
      /** Get stochastic basis associated with a particular deterministic column local id.
        */
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > getColStochasticBasis(int determLid) const
      { return sg_basis_col_dof_[determLid]; } 

   private:

      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > sg_master_basis_;
      std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sg_basis_row_dof_;
      std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sg_basis_col_dof_;

      std::vector<int> myRowGidOffsets_;
      std::vector<int> myColGidOffsets_;

      Teuchos::RCP<Epetra_CrsGraph> graph_;
      Teuchos::RCP<Epetra_Map> rowMap_;

      bool scaleOp_;
   };
    
} // namespace Stokhos

#endif // STOKHOS_AdaptivityUtils_HPP
