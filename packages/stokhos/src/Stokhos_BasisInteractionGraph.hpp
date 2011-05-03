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

#ifndef STOKHOS_BASIS_INTERACTION_GRAPH_HPP
#define STOKHOS_BASIS_INTERACTION_GRAPH_HPP

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_RCP.hpp"

#include <vector>

namespace Stokhos {

  class BasisInteractionGraph {
  public:
     BasisInteractionGraph();
     BasisInteractionGraph(const BasisInteractionGraph & flm);
     BasisInteractionGraph(const Stokhos::OrthogPolyBasis<int,double> & max_basis,bool onlyUseLinear=false,int porder=-1);
     BasisInteractionGraph(const Stokhos::ProductBasis<int,double> & masterBasis,
                           const Stokhos::ProductBasis<int,double> & rowBasis,
                           const Stokhos::ProductBasis<int,double> & colBasis,bool onlyUseLinear=false,int porder=-1);

     //! Setup the lookup graph
     void initialize(const Stokhos::OrthogPolyBasis<int,double> & max_basis,int porder=-1);

     //! Setup the lookup graph
     void initialize(const Stokhos::ProductBasis<int,double> & max_basis,
                     const Stokhos::ProductBasis<int,double> & rowBasis,
                     const Stokhos::ProductBasis<int,double> & colBasis,int porder=-1);

     //! Grab active indicies in graph for row i
     const std::vector<std::size_t> & activeIndices(std::size_t i) const
     { return vecLookup_[i]; }

     //! Is there an entry for (i,j) in the graph
     bool operator()(std::size_t i,std::size_t j) const;

     //! How many non zeros are in this graph
     std::size_t numNonZeros() const;

     //! What is the number of rows
     std::size_t rowCount() const
     { return vecLookup_.size(); }

     //! What is the number of columns
     std::size_t colCount() const
     { return numCols_; }

     void printGraph(std::ostream & os) const;

  protected:
     std::size_t numCols_;

     std::vector<std::vector<std::size_t> > vecLookup_;
 
     bool onlyUseLinear_;
  };
  
} // namespace Stokhos

#endif // STOKHOS_BASIS_INTERACTION_GRAPH_HPP
