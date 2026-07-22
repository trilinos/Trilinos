// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

#include <algorithm>
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

struct Coord {
  int i, j, k;
  int gid;
};

template <typename coord_t>
struct Tile {
  int lower, upper;
  Array<coord_t> parts;
};

typedef Tile<Coord> KTile;
typedef Tile<KTile> JTile;
typedef Tile<JTile> ITile;

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
    bool symmetric = true;
    CLP.setOption("symmetric", "asymmetric", &symmetric,
                  "Use basis polynomials with symmetric PDF");
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
    int i_tile_size = 128;
    CLP.setOption("tile_size", &i_tile_size, "Tile size");
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
    RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

    int basis_size = basis->size();
    std::cout << "basis size = " << basis_size
              << " num nonzero Cijk entries = " << Cijk->num_entries()
              << std::endl;

    // Build 2-way symmetric Cijk tensor
    RCP<Cijk_type> Cijk_sym = rcp(new Cijk_type);
    Cijk_type::i_iterator i_begin = Cijk->i_begin();
    Cijk_type::i_iterator i_end = Cijk->i_end();
    for (Cijk_type::i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      int i = index(i_it);
      Cijk_type::ik_iterator k_begin = Cijk->k_begin(i_it);
      Cijk_type::ik_iterator k_end = Cijk->k_end(i_it);
      for (Cijk_type::ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        int k = index(k_it);
        Cijk_type::ikj_iterator j_begin = Cijk->j_begin(k_it);
        Cijk_type::ikj_iterator j_end = Cijk->j_end(k_it);
        for (Cijk_type::ikj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          int j = index(j_it);
          if (k <= j) {
            double c = value(j_it);
            Cijk_sym->add_term(i, j, k, c);
          }
        }
      }
    }
    Cijk_sym->fillComplete();

    // First partition based on i
    int j_tile_size = i_tile_size / 2;
    int num_i_parts = (basis_size + i_tile_size-1) / i_tile_size;
    int its = basis_size / num_i_parts;
    Array<ITile> i_tiles(num_i_parts);
    for (int i=0; i<num_i_parts; ++i) {
      i_tiles[i].lower = i*its;
      i_tiles[i].upper = (i+1)*its;
      i_tiles[i].parts.resize(1);
      i_tiles[i].parts[0].lower = basis_size;
      i_tiles[i].parts[0].upper = 0;
    }

    // Next partition j
    for (Cijk_type::i_iterator i_it=Cijk_sym->i_begin();
         i_it!=Cijk_sym->i_end(); ++i_it) {
      int i = index(i_it);

      // Find which part i belongs to
      int idx = 0;
      while (idx < num_i_parts && i >= i_tiles[idx].lower) ++idx;
      --idx;
      TEUCHOS_ASSERT(idx >= 0 && idx < num_i_parts);

      Cijk_type::ik_iterator k_begin = Cijk_sym->k_begin(i_it);
      Cijk_type::ik_iterator k_end = Cijk_sym->k_end(i_it);
      for (Cijk_type::ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        int j = index(k_it);  // using symmetry to interchange j and k

        if (j < i_tiles[idx].parts[0].lower)
          i_tiles[idx].parts[0].lower = j;
        if (j > i_tiles[idx].parts[0].upper)
          i_tiles[idx].parts[0].upper = j;
      }
    }
    for (int idx=0; idx<num_i_parts; ++idx) {
      int lower = i_tiles[idx].parts[0].lower;
      int upper = i_tiles[idx].parts[0].upper;
      int range = upper - lower + 1;
      int num_j_parts = (range + j_tile_size-1) / j_tile_size;
      int jts = range / num_j_parts;
      Array<JTile> j_tiles(num_j_parts);
      for (int j=0; j<num_j_parts; ++j) {
        j_tiles[j].lower = lower + j*jts;
        j_tiles[j].upper = lower + (j+1)*jts;
        j_tiles[j].parts.resize(1);
        j_tiles[j].parts[0].lower = basis_size;
        j_tiles[j].parts[0].upper = 0;
      }
      i_tiles[idx].parts.swap(j_tiles);
    }

    // Now partition k
    for (Cijk_type::i_iterator i_it=Cijk_sym->i_begin();
         i_it!=Cijk_sym->i_end(); ++i_it) {
      int i = index(i_it);

      // Find which part i belongs to
      int idx = 0;
      while (idx < num_i_parts && i >= i_tiles[idx].lower) ++idx;
      --idx;
      TEUCHOS_ASSERT(idx >= 0 && idx < num_i_parts);

      Cijk_type::ik_iterator k_begin = Cijk_sym->k_begin(i_it);
      Cijk_type::ik_iterator k_end = Cijk_sym->k_end(i_it);
      for (Cijk_type::ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        int j = index(k_it);  // using symmetry to interchange j and k

        // Find which part j belongs to
        int num_j_parts = i_tiles[idx].parts.size();
        int jdx = 0;
        while (jdx < num_j_parts && j >= i_tiles[idx].parts[jdx].lower) ++jdx;
        --jdx;
        TEUCHOS_ASSERT(jdx >= 0 && jdx < num_j_parts);

        Cijk_type::ikj_iterator j_begin = Cijk_sym->j_begin(k_it);
        Cijk_type::ikj_iterator j_end = Cijk_sym->j_end(k_it);
        for (Cijk_type::ikj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          int k = index(j_it);  // using symmetry to interchange j and k

          if (k >= j) {
            Coord coord;
            coord.i = i; coord.j = j; coord.k = k;
            i_tiles[idx].parts[jdx].parts[0].parts.push_back(coord);
            if (k < i_tiles[idx].parts[jdx].parts[0].lower)
              i_tiles[idx].parts[jdx].parts[0].lower = k;
            if (k > i_tiles[idx].parts[jdx].parts[0].upper)
              i_tiles[idx].parts[jdx].parts[0].upper = k;
          }
        }
      }
    }

    // Now need to divide up k-parts based on lower/upper bounds
    int num_parts = 0;
    int num_coord = 0;
    for (int idx=0; idx<num_i_parts; ++idx) {
      int num_j_parts = i_tiles[idx].parts.size();
      for (int jdx=0; jdx<num_j_parts; ++jdx) {
        int lower = i_tiles[idx].parts[jdx].parts[0].lower;
        int upper = i_tiles[idx].parts[jdx].parts[0].upper;
        int range = upper - lower + 1;
        int num_k_parts = (range + j_tile_size-1) / j_tile_size;
        int kts = range / num_k_parts;
        Array<KTile> k_tiles(num_k_parts);
        for (int k=0; k<num_k_parts; ++k) {
          k_tiles[k].lower = lower + k*kts;
          k_tiles[k].upper = lower + (k+1)*kts;
        }
        int num_k = i_tiles[idx].parts[jdx].parts[0].parts.size();
        for (int l=0; l<num_k; ++l) {
          int i = i_tiles[idx].parts[jdx].parts[0].parts[l].i;
          int j = i_tiles[idx].parts[jdx].parts[0].parts[l].j;
          int k = i_tiles[idx].parts[jdx].parts[0].parts[l].k;

          // Find which part k belongs to
          int kdx = 0;
          while (kdx < num_k_parts && k >= k_tiles[kdx].lower) ++kdx;
          --kdx;
          TEUCHOS_ASSERT(kdx >= 0 && kdx < num_k_parts);

          Coord coord;
          coord.i = i; coord.j = j; coord.k = k;
          k_tiles[kdx].parts.push_back(coord);
          ++num_coord;
          if (j != k) ++num_coord;
        }

        // Eliminate parts with zero size
        Array<KTile> k_tiles2;
        for (int k=0; k<num_k_parts; ++k) {
          if (k_tiles[k].parts.size() > 0)
            k_tiles2.push_back(k_tiles[k]);
        }
        num_parts += k_tiles2.size();
        i_tiles[idx].parts[jdx].parts.swap(k_tiles2);
      }
    }
    TEUCHOS_ASSERT(num_coord == Cijk->num_entries());

    std::cout << "num parts requested = " << num_parts << std::endl;

    // Form list of part IDs
    Teuchos::Array<int> part_IDs(num_parts);
    for (int i=0; i<num_parts; ++i)
      part_IDs[i] = i;
    std::random_shuffle(part_IDs.begin(), part_IDs.end());

    // Assign part IDs
    int pp = 0;
    for (int idx=0; idx<num_i_parts; ++idx) {
      int num_j_parts = i_tiles[idx].parts.size();
      for (int jdx=0; jdx<num_j_parts; ++jdx) {
        int num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (int kdx=0; kdx<num_k_parts; ++kdx) {
          int num_k = i_tiles[idx].parts[jdx].parts[kdx].parts.size();
          for (int l=0; l<num_k; ++l) {
            i_tiles[idx].parts[jdx].parts[kdx].parts[l].gid = part_IDs[pp];
          }
          ++pp;
        }
      }
    }

    int part = 0;
    for (int idx=0; idx<num_i_parts; ++idx) {
      int num_j_parts = i_tiles[idx].parts.size();
      for (int jdx=0; jdx<num_j_parts; ++jdx) {
        int num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (int kdx=0; kdx<num_k_parts; ++kdx) {
          std::cout << part++ << " : ["
                    << i_tiles[idx].lower << ","
                    << i_tiles[idx].upper << ") x ["
                    << i_tiles[idx].parts[jdx].lower << ","
                    << i_tiles[idx].parts[jdx].upper << ") x ["
                    << i_tiles[idx].parts[jdx].parts[kdx].lower << ","
                    << i_tiles[idx].parts[jdx].parts[kdx].upper << ") : "
                    << i_tiles[idx].parts[jdx].parts[kdx].parts.size()
                    << std::endl;
        }
      }
    }

    // Print full 3-tensor to file
    std::ofstream cijk_file;
    if (save_3tensor) {
      cijk_file.open(file_3tensor.c_str());
      cijk_file.precision(14);
      cijk_file.setf(std::ios::scientific);
      cijk_file << "i, j, k, part" << std::endl;
      for (int idx=0; idx<num_i_parts; ++idx) {
        int num_j_parts = i_tiles[idx].parts.size();
        for (int jdx=0; jdx<num_j_parts; ++jdx) {
          int num_k_parts = i_tiles[idx].parts[jdx].parts.size();
          for (int kdx=0; kdx<num_k_parts; ++kdx) {
            int num_k = i_tiles[idx].parts[jdx].parts[kdx].parts.size();
            for (int l=0; l<num_k; ++l) {
              Coord c = i_tiles[idx].parts[jdx].parts[kdx].parts[l];
              cijk_file << c.i << ", " << c.j << ", " << c.k << ", " << c.gid
                        << std::endl;
              if (c.j != c.k)
                cijk_file << c.i << ", " << c.k << ", " << c.j << ", " << c.gid
                          << std::endl;
            }
          }
        }
      }
      cijk_file.close();
    }

  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
