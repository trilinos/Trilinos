// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos.hpp"
#include "Stokhos_Sparse3TensorPartition.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

#include <fstream>
#include <iostream>

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

// Symmetry types
const int num_symmetry_types = 3;
const Stokhos::CijkSymmetryType symmetry_type_values[] = {
  Stokhos::CIJK_NO_SYMMETRY,
  Stokhos::CIJK_TWO_WAY_SYMMETRY,
  Stokhos::CIJK_SIX_WAY_SYMMETRY
};
const char *symmetry_type_names[] = {
  "none", "2-way", "6-way" };

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example generates the sparsity pattern for the block stochastic Galerkin matrix.\n");
    int d = 5;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 3;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
    bool symmetric = true;
    CLP.setOption("symmetric", "asymmetric", &symmetric, "Use basis polynomials with symmetric PDF");
    Stokhos::GrowthPolicy growth_type = Stokhos::SLOW_GROWTH;
    CLP.setOption("growth", &growth_type,
                  num_growth_types, growth_type_values, growth_type_names,
                  "Growth type");
    ProductBasisType prod_basis_type = TOTAL;
    CLP.setOption("product_basis", &prod_basis_type,
                  num_prod_basis_types, prod_basis_type_values,
                  prod_basis_type_names,
                  "Product basis type");
    OrderingType ordering_type = LEXICOGRAPHIC_ORDERING;
    CLP.setOption("ordering", &ordering_type,
                  num_ordering_types, ordering_type_values,
                  ordering_type_names,
                  "Product basis ordering");
    Stokhos::CijkSymmetryType symmetry_type = Stokhos::CIJK_NO_SYMMETRY;
    CLP.setOption("symmetry", &symmetry_type,
                  num_symmetry_types, symmetry_type_values,
                  symmetry_type_names,
                  "Cijk symmetry type");
    bool full = true;
    CLP.setOption("full", "linear", &full, "Use full or linear expansion");
    int tile_size = 32;
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
    const double alpha = 1.0;
    const double beta = symmetric ? 1.0 : 2.0 ;
    for (int i=0; i<d; i++) {
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

    // File for saving sparse Cijk tensor and parts
    std::ofstream cijk_file;
    if (save_3tensor) {
      cijk_file.open(file_3tensor.c_str());
      cijk_file.precision(14);
      cijk_file.setf(std::ios::scientific);
      cijk_file << "i, j, k, part" << std::endl;
    }

    // Build tensor data list
    Teuchos::ArrayRCP< Stokhos::CijkData<int,double> > coordinate_list =
      Stokhos::build_cijk_coordinate_list(*Cijk, symmetry_type);

    // Partition via RCB
    typedef Stokhos::RCB< Stokhos::CijkData<int,double> > rcb_type;
    rcb_type rcb(tile_size, 10000, coordinate_list());
    int num_parts = rcb.get_num_parts();
    std::cout << "num parts = " << num_parts << std::endl;

    // Print part sizes
    RCP< Array< RCP<rcb_type::Box> > > parts = rcb.get_parts();
    for (int i=0; i<num_parts; ++i) {
      RCP<rcb_type::Box> box = (*parts)[i];
      std::cout << "part " << i << " bounding box ="
                << " [ " << box->delta_x << ", " << box->delta_y << ", "
                  << box->delta_z << " ]" << " nnz = "
                  << box->coords.size() << std::endl;
    }

    if (save_3tensor) {
      ArrayRCP<int> part_ids = rcb.get_part_IDs();
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
            if ((symmetry_type == Stokhos::CIJK_NO_SYMMETRY) ||
                (symmetry_type == Stokhos::CIJK_TWO_WAY_SYMMETRY && j >= k) ||
                (symmetry_type == Stokhos::CIJK_SIX_WAY_SYMMETRY && i >= j && j >= k)) {
              cijk_file << i << ", " << j << ", " << k << ", "
                        << part_ids[idx++] << std::endl;
            }
          }
        }
      }
      cijk_file.close();
    }

    //Teuchos::TimeMonitor::summarize(std::cout);

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
