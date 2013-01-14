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

#include "Stokhos_Epetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// sparsity_example
//
//  usage: 
//     sparsity_example [options]
//
//  output:  
//     prints the sparsity of the sparse 3 tensor specified by the basis,
//     dimension, order, given by summing over the third index, to a matrix
//     market file.  This sparsity pattern yields the sparsity of the block
//     stochastic Galerkin matrix which can be visualized, e.g., by matlab.
//     The full/linear flag determines whether the third index ranges over
//     the full polynomial dimension, or only over the zeroth and first order
//     terms.

// Basis types
enum BasisType { HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS, JACOBI };
const int num_basis_types = 6;
const BasisType basis_type_values[] = { 
  HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS, JACOBI };
const char *basis_type_names[] = { 
  "hermite", "legendre", "clenshaw-curtis", "gauss-patterson", "rys", "jacobi" };

// Growth policies
const int num_growth_types = 2;
const Stokhos::GrowthPolicy growth_type_values[] = {
  Stokhos::SLOW_GROWTH, Stokhos::MODERATE_GROWTH };
const char *growth_type_names[] = { "slow", "moderate" };

// Product Basis types
enum ProductBasisType { COMPLETE, TENSOR, TOTAL, SMOLYAK };
const int num_prod_basis_types = 4;
const ProductBasisType prod_basis_type_values[] = { 
  COMPLETE, TENSOR, TOTAL, SMOLYAK };
const char *prod_basis_type_names[] = { 
  "complete", "tensor", "total", "smolyak" };

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;

struct CijkNonzeros {
  int i, total_nz;
  Array< Array<int> > nz_tiles;
};

struct NZCompare {
  bool operator() (const CijkNonzeros& a, 
		   const CijkNonzeros& b) const {
    return (a.total_nz == b.total_nz) ? (a.i < b.i) : (a.total_nz > b.total_nz) ;
  }
};

struct NZPairCompare {
  bool operator() (const std::pair<int,int>& a, 
		   const std::pair<int,int>& b) const {
    return (a.second == b.second) ? (a.first < b.first) : (a.second > b.second);
  }
};

int main(int argc, char **argv)
{
  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example generates the sparsity pattern for the block stochastic Galerkin matrix.\n");
    int d = 3;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
    std::string file = "A.mm";
    CLP.setOption("filename", &file, "Matrix Market filename");
    BasisType basis_type = LEGENDRE;
    CLP.setOption("basis", &basis_type, 
		  num_basis_types, basis_type_values, basis_type_names, 
		  "Basis type");
    Stokhos::GrowthPolicy growth_type = Stokhos::SLOW_GROWTH;
    CLP.setOption("growth", &growth_type, 
		  num_growth_types, growth_type_values, growth_type_names, 
		  "Growth type");
    ProductBasisType prod_basis_type = COMPLETE;
    CLP.setOption("product_basis", &prod_basis_type, 
		  num_prod_basis_types, prod_basis_type_values, 
		  prod_basis_type_names, 
		  "Product basis type");
    double alpha = 1.0;
    CLP.setOption("alpha", &alpha, "Jacobi alpha index");
    double beta = 1.0;
    CLP.setOption("beta", &beta, "Jacobi beta index");
    bool full = true;
    CLP.setOption("full", "linear", &full, "Use full or linear expansion");
    bool use_old = false;
    CLP.setOption("old", "new", &use_old, "Use old or new Cijk algorithm");
    int tile_size = 100;
    CLP.setOption("tile_size", &tile_size, "Tile size");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      if (basis_type == HERMITE)
	bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(
				  p, true, growth_type));
      else if (basis_type == LEGENDRE)
	bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(
				  p, true, growth_type));
      else if (basis_type == CC_LEGENDRE)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<int,double>(
			 p, true));
      else if (basis_type == GP_LEGENDRE)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<int,double>(
			 p, true));
      else if (basis_type == RYS)
	bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(
				  p, 1.0, true, growth_type));
      else if (basis_type == JACOBI)
	bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<int,double>(
				  p, alpha, beta, true, growth_type));
    }
    RCP<const Stokhos::ProductBasis<int,double> > basis;
    if (prod_basis_type == COMPLETE)
      basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(
		       bases, drop, use_old));
    else if (prod_basis_type == TENSOR)
      basis = 
	Teuchos::rcp(new Stokhos::TensorProductBasis<int,double>(
		       bases, drop));
    else if (prod_basis_type == TOTAL)
      basis = 
	Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double>(
		       bases, drop));
    else if (prod_basis_type == SMOLYAK) {
      Stokhos::TotalOrderIndexSet<int> index_set(d, p);
      basis = 
	Teuchos::rcp(new Stokhos::SmolyakBasis<int,double>(
		       bases, index_set, drop));
    }

    // Triple product tensor
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    RCP<Cijk_type> Cijk;
    if (prod_basis_type == COMPLETE) {
      if (full)
	Cijk = basis->computeTripleProductTensor(basis->size());
      else
	Cijk = basis->computeTripleProductTensor(basis->dimension()+1);
    }
    else {
      if (full)
	Cijk = basis->computeTripleProductTensor(p);
      else
	Cijk = basis->computeTripleProductTensor(1);
    }

    int sz = basis->size();
    std::cout << "basis size = " << sz
	      << " num nonzero Cijk entries = " << Cijk->num_entries() 
	      << std::endl;

    // Setup tiles
    if (tile_size > sz)
      tile_size = sz;
    int j_sz = sz;
    int k_sz = sz;
    if (!full)
      k_sz = basis->dimension()+1;
    int nj_tiles = j_sz / tile_size;
    int nk_tiles = k_sz / tile_size;
    if (j_sz - nj_tiles*tile_size > 0)
      ++nj_tiles;
    if (k_sz - nk_tiles*tile_size > 0)
      ++nk_tiles;
    Array<CijkNonzeros> nz(sz);
    for (int i=0; i<sz; ++i) {
      nz[i].i = i;
      nz[i].nz_tiles.resize(nj_tiles);
      for (int j=0; j<nj_tiles; ++j)
	nz[i].nz_tiles[j].resize(nk_tiles);
    }

    // Get number of nonzeros in Cijk for each i
    Cijk_type::k_iterator k_begin = Cijk->k_begin();
    Cijk_type::k_iterator k_end = Cijk->k_end();
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = index(k_it);
      int k_tile = k / tile_size;
      Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
      Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	int j_tile = j / tile_size;
	Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
	for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	  int i = index(i_it);
	  ++nz[i].total_nz;
	  ++nz[i].nz_tiles[j_tile][k_tile];
	}
      }
    }

    // Sort based on total number of nonzeros
    std::sort(nz.begin(), nz.end(), NZCompare());
    
    // Print nonzeros
    int w_index = 3;
    int w_nz = 5;
    int w_tile = 4;
    for (int i=0; i<nz.size(); ++i) {
      int idx = nz[i].i;
      std::cout << setw(w_index) << idx << " " 
		<< basis->term(idx) << ": " 
		<< setw(w_nz) << nz[i].total_nz
		<< ", ";
      for (int j=0; j<nj_tiles; ++j)
	for (int k=0; k<nk_tiles; ++k)
	  std::cout << setw(w_tile) << nz[i].nz_tiles[j][k] << " ";
      std::cout << std::endl;
    }

    // Add up the nonzeros for each (j,k) tile
    Array< Array<int> > total_nz_tiles(nj_tiles);
    int total_nz = 0;
    for (int j=0; j<nj_tiles; ++j)
      total_nz_tiles[j].resize(nk_tiles);
    for (int i=0; i<nz.size(); ++i) {
      total_nz += nz[i].total_nz;
      for (int j=0; j<nj_tiles; ++j)
	for (int k=0; k<nk_tiles; ++k)
	  total_nz_tiles[j][k] += nz[i].nz_tiles[j][k];
    }
    int w_total = (w_index+1) + (2*basis->dimension()+5) + w_nz;
    std::cout << std::endl << setw(w_total) << total_nz << ", ";
    for (int j=0; j<nj_tiles; ++j)
      for (int k=0; k<nk_tiles; ++k)
	std::cout << setw(w_tile) << total_nz_tiles[j][k] << " ";
    std::cout << std::endl;

    // Now partition Cijk for each tile
    Array< Array< RCP<Cijk_type> > > Cijk_tile(nj_tiles);
    for (int j=0; j<nj_tiles; ++j) {
      Cijk_tile[j].resize(nk_tiles);
      for (int k=0; k<nk_tiles; ++k)
	Cijk_tile[j][k] = rcp(new Cijk_type);
    }
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = index(k_it);
      int k_tile = k / tile_size;
      Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
      Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	int j_tile = j / tile_size;
	Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
	for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	  int i = index(i_it);
	  double c = value(i_it);
	  Cijk_tile[j_tile][k_tile]->add_term(i,j,k,c);
	}
      }
    }
    for (int j=0; j<nj_tiles; ++j)
      for (int k=0; k<nk_tiles; ++k)
	Cijk_tile[j][k]->fillComplete();

    
    Array< Array< std::map<int,int> > > nz_tile(nj_tiles);
    Array< Array< Array< std::pair<int,int> > > > sorted_nz_tile(nj_tiles);
    for (int j_tile=0; j_tile<nj_tiles; ++j_tile) {
      nz_tile[j_tile].resize(nk_tiles); 
      sorted_nz_tile[j_tile].resize(nk_tiles); 
      for (int k_tile=0; k_tile<nk_tiles; ++k_tile) {

	// Count nonzeros for each i, for each tile
	Cijk_type::k_iterator k_begin = Cijk_tile[j_tile][k_tile]->k_begin();
	Cijk_type::k_iterator k_end = Cijk_tile[j_tile][k_tile]->k_end();
	for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
	  //int k = index(k_it);
	  Cijk_type::kj_iterator j_begin = 
	    Cijk_tile[j_tile][k_tile]->j_begin(k_it);
	  Cijk_type::kj_iterator j_end = 
	    Cijk_tile[j_tile][k_tile]->j_end(k_it);
	  for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	    //int j = index(j_it);
	    Cijk_type::kji_iterator i_begin = 
	      Cijk_tile[j_tile][k_tile]->i_begin(j_it);
	    Cijk_type::kji_iterator i_end = 
	      Cijk_tile[j_tile][k_tile]->i_end(j_it);
	    for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it){
	      int i = index(i_it);
	      if (nz_tile[j_tile][k_tile].count(i) == 0)
		nz_tile[j_tile][k_tile][i] = 1;
	      else
		++(nz_tile[j_tile][k_tile][i]);
	    }
	  }
	}

	// Sort based on non-zeros for each i, for each tile
	sorted_nz_tile[j_tile][k_tile].resize(nz_tile[j_tile][k_tile].size());
	int idx=0;
	for (std::map<int,int>::iterator it = nz_tile[j_tile][k_tile].begin();
	     it != nz_tile[j_tile][k_tile].end(); ++it) {
	  sorted_nz_tile[j_tile][k_tile][idx] = 
	    std::make_pair(it->first, it->second);
	  ++idx;
	}
	std::sort( sorted_nz_tile[j_tile][k_tile].begin(),
		   sorted_nz_tile[j_tile][k_tile].end(),
		   NZPairCompare() );

	// Print number of non-zeros for each i, for each tile
	std::cout << std::endl 
		  << "Tile (" << j_tile << ", " << k_tile << "):" << std::endl;
	for (int i=0; i<sorted_nz_tile[j_tile][k_tile].size(); ++i) {
	  int idx = sorted_nz_tile[j_tile][k_tile][i].first;
	  std::cout << setw(w_index) << idx << " " 
		    << basis->term(idx) << ": " 
		    << setw(w_nz) << sorted_nz_tile[j_tile][k_tile][i].second
		    << std::endl;
	  if (i % 32 == 31)
	    std::cout << std::endl;
	}
      }
    }
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
