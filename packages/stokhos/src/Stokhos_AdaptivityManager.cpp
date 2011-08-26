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

#include "Stokhos_AdaptivityManager.hpp"
#include "Stokhos_AdaptivityUtils.hpp"
#include "Stokhos_BasisInteractionGraph.hpp"

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_RowMatrixOut.h"

Stokhos::AdaptivityManager::AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_CrsGraph & determ_graph,
         bool onlyUseLinear,int kExpOrder,
         bool scaleOp)
   : sg_master_basis_(sg_master_basis), sg_basis_row_dof_(sg_basis_row_dof), scaleOp_(scaleOp)
{
    rowMap_ = adapt_utils::buildAdaptedRowMapAndOffsets(determ_graph.Comm(),sg_basis_row_dof_,myRowGidOffsets_);

    setupWithGraph(determ_graph,onlyUseLinear,kExpOrder);
}

Stokhos::AdaptivityManager::AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_Comm & comm,
         bool scaleOp)
   : sg_master_basis_(sg_master_basis), sg_basis_row_dof_(sg_basis_row_dof), scaleOp_(scaleOp)
{
   rowMap_ = adapt_utils::buildAdaptedRowMapAndOffsets(comm,sg_basis_row_dof_,myRowGidOffsets_);
}

Teuchos::RCP<Epetra_CrsMatrix> 
Stokhos::AdaptivityManager::buildMatrixFromGraph() const
{
   return Teuchos::rcp(new Epetra_CrsMatrix(Copy,*graph_));
}

void Stokhos::AdaptivityManager::setupWithGraph(const Epetra_CrsGraph & determGraph,bool onlyUseLinear,int kExpOrder) 
{
   graph_ = adapt_utils::buildAdaptedGraph(determGraph, sg_master_basis_, sg_basis_row_dof_, onlyUseLinear, kExpOrder);

   adapt_utils::buildAdaptedColOffsets(determGraph,myRowGidOffsets_,myColGidOffsets_);
   adapt_utils::buildColBasisFunctions(determGraph,sg_master_basis_,sg_basis_row_dof_,sg_basis_col_dof_);
}

/** Setup operator
  */
void 
Stokhos::AdaptivityManager::
setupOperator(Epetra_CrsMatrix & A,const Sparse3Tensor<int,double> & Cijk,Stokhos::EpetraOperatorOrthogPoly & poly,
              bool onlyUseLinear,bool includeMean) const
{
   typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

   // Zero out matrix
   A.PutScalar(0.0);

   // Compute loop bounds
   Cijk_type::k_iterator k_begin = Cijk.k_begin();
   Cijk_type::k_iterator k_end = Cijk.k_end();
   if (!includeMean && index(k_begin) == 0)
     ++k_begin;
   if (onlyUseLinear) {
     int dim = sg_master_basis_->dimension();
     k_end = Cijk.find_k(dim+1);
   }
 
   // Assemble matrix
   // const Teuchos::Array<double>& norms = sg_basis->norm_squared();
   for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = index(k_it);
      Teuchos::RCP<Epetra_CrsMatrix> block = 
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(poly.getCoeffPtr(k), 
 						    true);

      // add in matrix k
      sumInOperator(A,Cijk,k,*block);
   }
}

void
Stokhos::AdaptivityManager::
sumInOperator(Epetra_CrsMatrix & A,const Stokhos::Sparse3Tensor<int,double> & Cijk,int k,const Epetra_CrsMatrix & J_k) const
{
   TEUCHOS_ASSERT(J_k.NumMyRows() == int(sg_basis_row_dof_.size()));
   TEUCHOS_ASSERT(J_k.NumMyCols() == int(sg_basis_col_dof_.size()));

   const Teuchos::Array<double> & normValues = sg_master_basis_->norm_squared();

   // loop over deterministic rows 
   for(int localM=0;localM<J_k.NumMyRows();localM++) {
      int m = J_k.GRID(localM);

      // grab row basis
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > rowStochBasis 
            = sg_basis_row_dof_[localM]; 
 
      // grab row from deterministic system
      int d_numEntries;
      int * d_Indices;
      double * d_Values;
     
      J_k.ExtractMyRowView(m,d_numEntries,d_Values,d_Indices);
      
      // loop over stochastic degrees of freedom of this row
      for(int rb_i=0;rb_i<rowStochBasis->size();rb_i++) {
         int i = sg_master_basis_->getIndex(rowStochBasis->getTerm(rb_i));

         double normValue = normValues[i]; // sg_master_basis->norm_squared(i);
         
         int sg_m = getGlobalRowId(localM,rb_i);

         // we wipe out old values, capacity should gurantee
         // we don't allocate more often than neccessary!
         std::vector<int> sg_indices;
         std::vector<double> sg_values;

         // sg_indices.resize(0); 
         // sg_values.resize(0);

         // loop over each column
         for(int colInd=0;colInd<d_numEntries;colInd++) {
            int n = d_Indices[colInd]; // grab global deterministic column id
            int localN = J_k.LCID(n);  // grab local deterministic column id

            // grab row basis
            Teuchos::RCP<const Stokhos::ProductBasis<int,double> > colStochBasis 
                  = sg_basis_col_dof_[localN]; 

            // build values array
            for(int cb_j=0;cb_j<colStochBasis->size();cb_j++) {
               int j = sg_master_basis_->getIndex(colStochBasis->getTerm(cb_j));
               int sg_n = getGlobalColId(localN,cb_j);
               double cijk = Cijk.getValue(i,j,k); 

               // no reason to work it in!
               if(cijk==0) continue;

               if(scaleOp_)
                  cijk = cijk/normValue;

               sg_indices.push_back(sg_n);
               sg_values.push_back(cijk*d_Values[colInd]);
            }
         }

         // add in matrix values
         A.SumIntoGlobalValues(sg_m,sg_indices.size(),&sg_values[0],&sg_indices[0]);
      }
   }
}

/** Copy to an adaptive vector from a set of blocked vectors
  */
void Stokhos::AdaptivityManager::copyToAdaptiveVector(const Stokhos::EpetraVectorOrthogPoly & x_sg,Epetra_Vector & x) const
{
   Teuchos::RCP<const EpetraExt::BlockVector> x_sg_bv = x_sg.getBlockVector();

   // copy from adapted vector to deterministic
   for(std::size_t i=0;i<sg_basis_row_dof_.size();i++) {
      int P_i = getRowStochasticBasisSize(i); 
      int localId = rowMap_->LID(getGlobalRowId(i,0));

      for(int j=0;j<P_i;j++,localId++) {
         int blk = sg_master_basis_->getIndex(sg_basis_row_dof_[i]->getTerm(j));
         x[localId] = x_sg_bv->GetBlock(blk)->operator[](i);
      }
   }
}

/** Copy from an adaptive vector to a set of blocked vectors
  */
void Stokhos::AdaptivityManager::copyFromAdaptiveVector(const Epetra_Vector & x,Stokhos::EpetraVectorOrthogPoly & x_sg) const
{
   int numBlocks = x_sg.size();
   Teuchos::RCP<EpetraExt::BlockVector> x_sg_bv = x_sg.getBlockVector();

   // zero out determinstic vectors
   for(int blk=0;blk<numBlocks;blk++)
      x_sg_bv->GetBlock(blk)->PutScalar(0.0);

   // copy from adapted vector to deterministic
   for(std::size_t i=0;i<sg_basis_row_dof_.size();i++) {
      int P_i = getRowStochasticBasisSize(i); 
      int localId = rowMap_->LID(getGlobalRowId(i,0));

      for(int j=0;j<P_i;j++,localId++) {
         int blk = sg_master_basis_->getIndex(sg_basis_row_dof_[i]->getTerm(j));
         x_sg_bv->GetBlock(blk)->operator[](i) = x[localId];
      }
   }
}
