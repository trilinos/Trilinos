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
#include "Teuchos_ParameterList.hpp"
#include "Isorropia_EpetraPartitioner.hpp"

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

// Ordering types
enum OrderingType { TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING };
const int num_ordering_types = 2;
const OrderingType ordering_type_values[] = {
  TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING };
const char *ordering_type_names[] = {
  "total", "lexicographic" };

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::Array;

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
    OrderingType ordering_type = TOTAL_ORDERING;
    CLP.setOption("ordering", &ordering_type,
                  num_ordering_types, ordering_type_values,
                  ordering_type_names,
                  "Product basis ordering");
    double alpha = 1.0;
    CLP.setOption("alpha", &alpha, "Jacobi alpha index");
    double beta = 1.0;
    CLP.setOption("beta", &beta, "Jacobi beta index");
    bool full = true;
    CLP.setOption("full", "linear", &full, "Use full or linear expansion");
    int tile_size = 128;
    CLP.setOption("tile_size", &tile_size, "Tile size");
    bool save_3tensor = false;
    CLP.setOption("save_3tensor", "no-save_3tensor", &save_3tensor,
                  "Save full 3tensor to file");
    std::string file_3tensor = "Cijk.dat";
    CLP.setOption("filename_3tensor", &file_3tensor,
                  "Filename to store full 3-tensor");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
      if (basis_type == HERMITE)
        bases[i] = rcp(new Stokhos::HermiteBasis<int,double>(
                                  p, true, growth_type));
      else if (basis_type == LEGENDRE)
        bases[i] = rcp(new Stokhos::LegendreBasis<int,double>(
                                  p, true, growth_type));
      else if (basis_type == CC_LEGENDRE)
        bases[i] =
          rcp(new Stokhos::ClenshawCurtisLegendreBasis<int,double>(
                         p, true));
      else if (basis_type == GP_LEGENDRE)
        bases[i] =
          rcp(new Stokhos::GaussPattersonLegendreBasis<int,double>(
                         p, true));
      else if (basis_type == RYS)
        bases[i] = rcp(new Stokhos::RysBasis<int,double>(
                                  p, 1.0, true, growth_type));
      else if (basis_type == JACOBI)
        bases[i] = rcp(new Stokhos::JacobiBasis<int,double>(
                                  p, alpha, beta, true, growth_type));
    }
    RCP<const Stokhos::ProductBasis<int,double> > basis;
    typedef Stokhos::TotalOrderLess< Stokhos::MultiIndex<int> > total_less;
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > lexo_less;
    if (prod_basis_type == COMPLETE)
      basis =
        rcp(new Stokhos::CompletePolynomialBasis<int,double>(
                       bases, drop));
    else if (prod_basis_type == TENSOR) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          rcp(new Stokhos::TensorProductBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          rcp(new Stokhos::TensorProductBasis<int,double,lexo_less>(
                         bases, drop));
    }

    else if (prod_basis_type == TOTAL) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          rcp(new Stokhos::TotalOrderBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          rcp(new Stokhos::TotalOrderBasis<int,double,lexo_less>(
                         bases, drop));
    }
    else if (prod_basis_type == SMOLYAK) {
      Stokhos::TotalOrderIndexSet<int> index_set(d, p);
      if (ordering_type == TOTAL_ORDERING)
        basis =
          rcp(new Stokhos::SmolyakBasis<int,double,total_less>(
                         bases, index_set, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          rcp(new Stokhos::SmolyakBasis<int,double,lexo_less>(
                         bases, index_set, drop));
    }

    // Triple product tensor
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    RCP<Cijk_type> Cijk;
    if (full)
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();

    int basis_size = basis->size();
    std::cout << "basis size = " << basis_size
              << " num nonzero Cijk entries = " << Cijk->num_entries()
              << std::endl;

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    // Store Cijk (i,j,k) triples in Epetra_MultiVector
    int num_cijk_entries = Cijk->num_entries();
    Epetra_LocalMap cijk_map(num_cijk_entries, 0, comm);
    Epetra_MultiVector ijk_triples(cijk_map, 3);
    int idx = 0;
    Cijk_type::k_iterator k_begin = Cijk->k_begin();
    Cijk_type::k_iterator k_end = Cijk->k_end();
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = index(k_it);
      Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
      Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        int j = index(j_it);
        Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
        Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
        for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
          int i = index(i_it);
          ijk_triples[0][idx] = i;
          ijk_triples[1][idx] = j;
          ijk_triples[2][idx] = k;
          ++idx;
        }
      }
    }

    // Partition ijk_triples using isorropia
    ParameterList params;
    params.set("partitioning method", "rcb");
    int num_parts = num_cijk_entries / tile_size;
    if (num_cijk_entries % tile_size > 0)
      ++num_parts;
    std::stringstream ss;
    ss << num_parts;
    params.set<std::string>("num parts", ss.str());
    RCP<const Epetra_MultiVector> ijk_triples_rcp =
      rcp(&ijk_triples,false);
    Isorropia::Epetra::Partitioner partitioner(ijk_triples_rcp, params);
    partitioner.compute();

    std::cout << "num parts requested = " << num_parts 
              << " num parts computed = " << partitioner.numProperties() << std::endl;

    Teuchos::Array<int> parts(num_cijk_entries);
    int sz;
    partitioner.extractPartsCopy(num_cijk_entries, sz, &parts[0]);
    TEUCHOS_ASSERT(sz == num_cijk_entries);

    // Print full 3-tensor to file
    if (save_3tensor) {
      std::ofstream cijk_file(file_3tensor.c_str());
      cijk_file.precision(14);
      cijk_file.setf(std::ios::scientific);
      cijk_file << "i, j, k, part" << std::endl;
      for (int i=0; i<num_cijk_entries; ++i) {
        cijk_file << ijk_triples[0][i] << ", "
                  << ijk_triples[1][i] << ", "
                  << ijk_triples[2][i] << ", "
                  << parts[i] << std::endl;
      }
      cijk_file.close();
    }

    Teuchos::TimeMonitor::summarize(std::cout);

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
