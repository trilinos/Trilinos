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

#include "Stokhos_BasisInteractionGraph.hpp"

Stokhos::BasisInteractionGraph::BasisInteractionGraph() 
{}

Stokhos::BasisInteractionGraph::BasisInteractionGraph(const BasisInteractionGraph & flm)
   : vecLookup_(flm.vecLookup_), onlyUseLinear_(flm.onlyUseLinear_)
{ }
   
Stokhos::BasisInteractionGraph::BasisInteractionGraph(const Stokhos::OrthogPolyBasis<int,double> & max_basis,bool onlyUseLinear,int porder)
   : onlyUseLinear_(onlyUseLinear)
{ 
   using Teuchos::RCP;

   if(porder<0)
      porder = max_basis.size();
   RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = max_basis.computeTripleProductTensor(porder);
   initialize(max_basis,*Cijk,porder); 
}

Stokhos::BasisInteractionGraph::BasisInteractionGraph(const Stokhos::ProductBasis<int,double> & masterBasis,
                                                      const Stokhos::ProductBasis<int,double> & rowBasis,
                                                      const Stokhos::ProductBasis<int,double> & colBasis,bool onlyUseLinear,int porder)
   : onlyUseLinear_(onlyUseLinear)
{
   using Teuchos::RCP;

   if(porder<0)
      porder = masterBasis.size();
   RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = masterBasis.computeTripleProductTensor(porder);
   initialize(masterBasis,*Cijk,rowBasis,colBasis,porder);
}

Stokhos::BasisInteractionGraph::BasisInteractionGraph(const Stokhos::OrthogPolyBasis<int,double> & max_basis,
                                                      const Stokhos::Sparse3Tensor<int,double> & Cijk,
                                                      bool onlyUseLinear,int porder)
   : onlyUseLinear_(onlyUseLinear)
{ 
   initialize(max_basis,Cijk,porder); 
}

Stokhos::BasisInteractionGraph::BasisInteractionGraph(const Stokhos::ProductBasis<int,double> & masterBasis,
                                                      const Stokhos::Sparse3Tensor<int,double> & Cijk,
                                                      const Stokhos::ProductBasis<int,double> & rowBasis,
                                                      const Stokhos::ProductBasis<int,double> & colBasis,bool onlyUseLinear,int porder)
   : onlyUseLinear_(onlyUseLinear)
{
   using Teuchos::RCP;

   initialize(masterBasis,Cijk,rowBasis,colBasis,porder);
}

void Stokhos::BasisInteractionGraph::initialize(const Stokhos::OrthogPolyBasis<int,double> & max_basis,
                                                const Stokhos::Sparse3Tensor<int,double> & Cijk,
                                                int porder)
{
   using Teuchos::RCP;
   typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

   // // max it out if defaulted
   // if(porder<0)
   //    porder = max_basis.size();

   // RCP<Stokhos::Sparse3Tensor<int,double> > Cijk = max_basis.computeTripleProductTensor(porder);

   Cijk_type::k_iterator k_end = Cijk.k_end();
   if (onlyUseLinear_) {
      int dim = max_basis.dimension();
      k_end = Cijk.find_k(dim+1);
   }

   vecLookup_.resize(max_basis.size()); // defines number of rows
   numCols_ = vecLookup_.size(); // set number of columns

   // Loop over Cijk entries including a non-zero in the graph at
   // indices (i,j) if there is any k for which Cijk is non-zero
   for(Cijk_type::k_iterator k_it=Cijk.k_begin(); k_it!=k_end; ++k_it) {
      for(Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it); ++j_it) {
         int j = index(j_it);
         for(Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it); i_it != Cijk.i_end(j_it); ++i_it) {
            int i = index(i_it);
            vecLookup_[i].push_back(j);
         }
      }
   }
}

void Stokhos::BasisInteractionGraph::initialize(const Stokhos::ProductBasis<int,double> & masterBasis,
                                                const Stokhos::Sparse3Tensor<int,double> & Cijk,
                                                const Stokhos::ProductBasis<int,double> & rowBasis,
                                                const Stokhos::ProductBasis<int,double> & colBasis,int porder)
{
   // for determining if their is an interaction or not
   Stokhos::BasisInteractionGraph masterBig(masterBasis,Cijk,onlyUseLinear_,porder);

   vecLookup_.resize(rowBasis.size()); // defines number of rows

   // set number of columns
   numCols_ = colBasis.size();

   // build row basis terms
   std::vector<int> rowIndexToMasterIndex(rowBasis.size());
   for(int i=0;i<rowBasis.size();i++) 
      rowIndexToMasterIndex[i] = masterBasis.index(rowBasis.term(i));

   // build column basis terms
   std::vector<int> colIndexToMasterIndex(colBasis.size());
   for(int i=0;i<colBasis.size();i++) 
      colIndexToMasterIndex[i] = masterBasis.index(colBasis.term(i));

   // build graph by looking up sparsity in master basis
   for(int r=0;r<rowBasis.size();r++) {
      int masterRow = rowIndexToMasterIndex[r];
      for(int c=0;c<colBasis.size();c++) {
         int masterCol = colIndexToMasterIndex[c];

         // is row and column active in master element
         bool activeRC = masterBig(masterRow,masterCol); 

         // if active add to local graph
         if(activeRC)
            vecLookup_[r].push_back(c);
      } 
   }
}

bool Stokhos::BasisInteractionGraph::operator()(std::size_t i,std::size_t j) const
{
   const std::vector<std::size_t> & indices = activeIndices(i);
   return indices.end() != std::find(indices.begin(),indices.end(),j);
}

void Stokhos::BasisInteractionGraph::printGraph(std::ostream & os) const
{
   for(std::size_t r=0;r<rowCount();r++) {
      for(std::size_t c=0;c<colCount();c++)
         if(operator()(r,c)) os << " * ";
         else os << "   ";
      os << std::endl;
   }
}

std::size_t Stokhos::BasisInteractionGraph::numNonZeros() const
{
   std::size_t nnz = 0;
   for(std::size_t r=0;r<rowCount();r++)
      nnz += activeIndices(r).size();
   return nnz;
}
