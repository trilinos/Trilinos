// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
