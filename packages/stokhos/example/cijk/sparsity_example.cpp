// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_Epetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Ifpack_RCMReordering.h"

#include <fstream>
#include <iostream>

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
enum OrderingType { TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING, MORTON_Z_ORDERING };
const int num_ordering_types = 3;
const OrderingType ordering_type_values[] = {
  TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING, MORTON_Z_ORDERING };
const char *ordering_type_names[] = {
  "total", "lexicographic", "morton-z" };

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
    bool use_old = false;
    CLP.setOption("old", "new", &use_old, "Use old or new Cijk algorithm");
    bool print = false;
    CLP.setOption("print", "no-print", &print, "Print Cijk to screen");
    bool save_3tensor = false;
    CLP.setOption("save_3tensor", "no-save_3tensor", &save_3tensor,
                  "Save full 3tensor to file");
    std::string file_3tensor = "Cijk.dat";
    CLP.setOption("filename_3tensor", &file_3tensor,
                  "Filename to store full 3-tensor");
    bool unique = false;
    CLP.setOption("unique", "no-unique", &unique,
                  "Only save the unique non-zeros");
    bool rcm = false;
    CLP.setOption("rcm", "no-rcm", &rcm, "Reorder using RCM");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
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
    Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;
    typedef Stokhos::TotalOrderLess< Stokhos::MultiIndex<int> > total_less;
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > lexo_less;
    typedef Stokhos::MortonZLess< Stokhos::MultiIndex<int> > z_less;
    if (prod_basis_type == COMPLETE)
      basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(
                       bases, drop, use_old));
    else if (prod_basis_type == TENSOR) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,lexo_less>(
                         bases, drop));
      else if (ordering_type == MORTON_Z_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,z_less>(
                         bases, drop));
    }
    else if (prod_basis_type == TOTAL) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,lexo_less>(
                         bases, drop));
      else if (ordering_type == MORTON_Z_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,z_less>(
                         bases, drop));
    }
    else if (prod_basis_type == SMOLYAK) {
      Stokhos::TotalOrderIndexSet<int> index_set(d, p);
       if (ordering_type == TOTAL_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,total_less>(
                          bases, index_set, drop));
       else if (ordering_type == LEXICOGRAPHIC_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,lexo_less>(
                          bases, index_set, drop));
       else if (ordering_type == MORTON_Z_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,z_less>(
                          bases, index_set, drop));
    }

    // Triple product tensor
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    Teuchos::RCP<Cijk_type> Cijk;
    if (full)
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();

    std::cout << "basis size = " << basis->size()
              << " num nonzero Cijk entries = " << Cijk->num_entries()
              << std::endl;

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    if (rcm) {
      Teuchos::RCP<Cijk_type> Cijk3 = Teuchos::rcp(new Cijk_type);
      {
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
            double cijk = value(i_it);
            if (i != 0 || (i == 0 && j == 0 && k == 0))
              Cijk3->add_term(i, j, k, cijk);
          }
        }
      }
      }
      Cijk3->fillComplete();

      Teuchos::RCP<Epetra_CrsGraph> graph =
        Stokhos::sparse3Tensor2CrsGraph(*basis, *Cijk3, comm);
      Epetra_CrsMatrix mat(Copy, *graph);
      mat.FillComplete();
      mat.PutScalar(1.0);
      Ifpack_RCMReordering ifpack_rcm;
      ifpack_rcm.SetParameter("reorder: root node", basis->size()-1);
      ifpack_rcm.Compute(mat);

      Teuchos::RCP<Cijk_type> Cijk2 = Teuchos::rcp(new Cijk_type);
      Cijk_type::k_iterator k_begin = Cijk->k_begin();
      Cijk_type::k_iterator k_end = Cijk->k_end();
      for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
        int k = ifpack_rcm.Reorder(index(k_it));
        Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
        Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
        for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          int j = ifpack_rcm.Reorder(index(j_it));
          Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
          Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
          for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
            int i = ifpack_rcm.Reorder(index(i_it));
            double cijk = value(i_it);
            Cijk2->add_term(i, j, k, cijk);
          }
        }
      }
      Cijk2->fillComplete();
      Cijk = Cijk2;
    }

    if (print) {
      std::cout << *Cijk << std::endl;
    }

    // Print triple product sparsity to matrix market file
    Stokhos::sparse3Tensor2MatrixMarket(*basis, *Cijk, comm, file);

    // Print full 3-tensor to file
    if (save_3tensor) {
      std::ofstream cijk_file(file_3tensor.c_str());
      cijk_file.precision(14);
      cijk_file.setf(std::ios::scientific);
      cijk_file << "i, j, k, cijk" << std::endl;
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
            double cijk = value(i_it);
            if (!unique || ( i >= j && j >= k ))
              cijk_file << i << ", "
                        << j << ", "
                        << k << ", "
                        << cijk << std::endl;
          }
        }
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
