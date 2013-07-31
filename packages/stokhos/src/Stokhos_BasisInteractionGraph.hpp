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

     BasisInteractionGraph(const Stokhos::OrthogPolyBasis<int,double> & max_basis,
                           const Stokhos::Sparse3Tensor<int,double> & Cijk,
                           bool onlyUseLinear=false,int porder=-1);

     BasisInteractionGraph(const Stokhos::ProductBasis<int,double> & masterBasis,
                           const Stokhos::Sparse3Tensor<int,double> & Cijk,
                           const Stokhos::ProductBasis<int,double> & rowBasis,
                           const Stokhos::ProductBasis<int,double> & colBasis,bool onlyUseLinear=false,int porder=-1);

     //! Setup the lookup graph
     void initialize(const Stokhos::OrthogPolyBasis<int,double> & max_basis,
                     const Stokhos::Sparse3Tensor<int,double> & Cijk,
                     int porder=-1);

     //! Setup the lookup graph
     void initialize(const Stokhos::ProductBasis<int,double> & max_basis,
                     const Stokhos::Sparse3Tensor<int,double> & Cijk,
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
