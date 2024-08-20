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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_toString.hpp"

#include <fstream>
#include <iostream>

extern "C" {
#include "zoltan.h"
}

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

// Partitioning types
enum PartitioningType { RCB, HG_FLAT_J };
const int num_partitioning_types = 2;
const PartitioningType partitioning_type_values[] = {
  RCB, HG_FLAT_J };
const char *partitioning_type_names[] = {
  "rcb", "hg_flat_j" };

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::Array;
using Teuchos::toString;

struct TensorData {
  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
  RCP<const Stokhos::ProductBasis<int,double> > basis;
  RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
};

// Functions implementing hypergraph for 3-D decomposition
// For this hypergraph model
//   * the nnz vertices are the Cijk non-zeros
//   * the 3*n hyperedges are the i, j, and k values
//       each Cijk non-zero belongs to 3 hyperedges:  i, j, and k
namespace HG_3D {

  // Return number of vertices
  int get_number_of_vertices(void *data, int *ierr) {
    TensorData *td = static_cast<TensorData*>(data);
    *ierr = ZOLTAN_OK;

    return td->Cijk->num_entries();
  }

  // Compute IDs and weights of each vertex
  void get_vertex_list(void *data, int sizeGID, int sizeLID,
                       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int wgt_dim, float *obj_wgts, int *ierr) {
    TensorData *td = static_cast<TensorData*>(data);
    *ierr = ZOLTAN_OK;

    int nnz = td->Cijk->num_entries();
    for (int i=0; i<nnz; ++i) {
      globalID[i] = i;
      localID[i] = i;
    }

    // Do not set weights so Zoltan assumes equally weighted vertices
  }

  // Compute number of hyperedges and pins
  void get_hypergraph_size(void *data, int *num_lists, int *num_pins,
                           int *format, int *ierr) {
    TensorData *td = static_cast<TensorData*>(data);
    *ierr = ZOLTAN_OK;

    //int n = td->basis->size();
    int nnz = td->Cijk->num_entries();

    // Number of vertices
    *num_lists = nnz;

    // Number of pins.  Each nonzero belongs creates 1 pin in 3 hyperedges
    // thus there are 3*nnz pins
    *num_pins = 3*nnz;

    // hypergraph will be stored in compressed-vertex format
    *format = ZOLTAN_COMPRESSED_VERTEX;
  }

  // Compute hypergraph
  void get_hypergraph(void *data, int sizeGID, int num_vtx, int num_pins,
                      int format, ZOLTAN_ID_PTR vtxGID, int *edgePtr,
                      ZOLTAN_ID_PTR edgeGID, int *ierr) {
    TensorData *td = static_cast<TensorData*>(data);
    *ierr = ZOLTAN_OK;

    int n = td->basis->size();

    TEUCHOS_ASSERT(sizeGID == 1);
    TEUCHOS_ASSERT(num_vtx == td->Cijk->num_entries());
    TEUCHOS_ASSERT(num_pins == 3*(td->Cijk->num_entries()));

    // Compute pins in each hyperedge stored in compressed-vertex format.
    // For each vertex we store the GIDs of the 3 edges that it connects to.
    // Edges are ordered as follows:
    //    [0,n)     -- i edges
    //    [n,2*n)   -- j edges
    //    [2*n,3*n) -- k edges
    int vtx_idx = 0;
    int pin_idx = 0;
    TensorData::Cijk_type::k_iterator k_begin = td->Cijk->k_begin();
    TensorData::Cijk_type::k_iterator k_end =   td->Cijk->k_end();
    for (TensorData::Cijk_type::k_iterator k_it=k_begin; k_it!=k_end;
         ++k_it) {
      int k = index(k_it);
       TensorData::Cijk_type::kj_iterator j_begin = td->Cijk->j_begin(k_it);
       TensorData::Cijk_type::kj_iterator j_end =   td->Cijk->j_end(k_it);
       for (TensorData::Cijk_type::kj_iterator j_it = j_begin; j_it != j_end;
            ++j_it) {
         int j = index(j_it);
         TensorData::Cijk_type::kji_iterator i_begin = td->Cijk->i_begin(j_it);
         TensorData::Cijk_type::kji_iterator i_end =   td->Cijk->i_end(j_it);
         for (TensorData::Cijk_type::kji_iterator i_it = i_begin; i_it != i_end;
              ++i_it) {
           int i = index(i_it);
           vtxGID[vtx_idx] = vtx_idx;
           edgePtr[vtx_idx++] = pin_idx;
           edgeGID[pin_idx++] = i;
           edgeGID[pin_idx++] = n + j;
           edgeGID[pin_idx++] = 2*n + k;
         }
       }
    }
  }
}


int main(int argc, char **argv)
{
  try {

    // Initialize Zoltan
    float version;
    int rc = Zoltan_Initialize(argc,argv,&version);
    TEUCHOS_ASSERT(rc == 0);

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
    PartitioningType partitioning_type = RCB;
    CLP.setOption("partitioning", &partitioning_type,
                  num_partitioning_types, partitioning_type_values,
                  partitioning_type_names,
                  "Partitioning Method");
    double imbalance_tol = 1.2;
    CLP.setOption("imbalance", &imbalance_tol, "Imbalance tolerance");
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

    // Create zoltan
    Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

    // Setup Zoltan parameters
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");

    // partitioning method
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); // version of method
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");// global IDs are integers
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");// local IDs are integers
    //Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); // export AND import lists
    Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); // use Zoltan default vertex weights
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");// use Zoltan default hyperedge weights
    int num_parts = basis_size / tile_size;
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", toString(num_parts).c_str());
    Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", toString(num_parts).c_str());
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", toString(imbalance_tol).c_str());

    // Set query functions
    TensorData td; td.basis = basis; td.Cijk = Cijk;
    Zoltan_Set_Num_Obj_Fn(zz, HG_3D::get_number_of_vertices, &td);
    Zoltan_Set_Obj_List_Fn(zz, HG_3D::get_vertex_list, &td);
    Zoltan_Set_HG_Size_CS_Fn(zz, HG_3D::get_hypergraph_size, &td);
    Zoltan_Set_HG_CS_Fn(zz, HG_3D::get_hypergraph, &td);

    // Partition
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    rc =
      Zoltan_LB_Partition(
        zz, // input (all remaining fields are output)
        &changes,        // 1 if partitioning was changed, 0 otherwise
        &numGidEntries,  // Number of integers used for a global ID
        &numLidEntries,  // Number of integers used for a local ID
        &numImport,      // Number of vertices to be sent to me
        &importGlobalGids,  // Global IDs of vertices to be sent to me
        &importLocalGids,   // Local IDs of vertices to be sent to me
        &importProcs,    // Process rank for source of each incoming vertex
        &importToPart,   // New partition for each incoming vertex
        &numExport,      // Number of vertices I must send to other processes*/
        &exportGlobalGids,  // Global IDs of the vertices I must send
        &exportLocalGids,   // Local IDs of the vertices I must send
        &exportProcs,    // Process to which I send each of the vertices
        &exportToPart);  // Partition to which each vertex will belong
    TEUCHOS_ASSERT(rc == 0);

    std::cout << "num parts requested = " << num_parts
              << " changes= " << changes
              << " num import = " << numImport
              << " num export = " << numExport << std::endl;

    // Build list of rows that belong to each part based on diagonal
    Array< Array<int> > part_map(num_parts);
    int idx = 0;
    int num_diag = 0;
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
          if (i == j && j == k) {
            part_map[ exportToPart[idx] ].push_back(i);
            ++num_diag;
          }
          idx++;
        }
      }
    }

    std::cout << "basis_size = " << basis_size << " num_diag = " << num_diag
              << std::endl;

    // Build permuation array mapping reoredered to original
    Array<int> perm_new_to_old;
    for (int part=0; part<num_parts; ++part) {
      int num_row = part_map[part].size();
      for (int i=0; i<num_row; ++i)
        perm_new_to_old.push_back(part_map[part][i]);
    }
    TEUCHOS_ASSERT(perm_new_to_old.size() == basis_size);

    // Build permuation array mapping original to reordered
    Array<int> perm_old_to_new(basis_size);
    for (int i=0; i<basis_size; ++i)
      perm_old_to_new[ perm_new_to_old[i] ] = i;

    if (save_3tensor) {
      idx = 0;
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
            cijk_file << perm_old_to_new[i] << ", "
                      << perm_old_to_new[j] << ", "
                      << perm_old_to_new[k] << ", "
                      << exportToPart[idx++] << std::endl;
            // cijk_file << i << ", "
            //           << j << ", "
            //           << k << ", "
            //           << exportToPart[idx++] << std::endl;
          }
        }
      }
      cijk_file.close();
    }

    // Clean-up
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                        &exportProcs, &exportToPart);
    Zoltan_Destroy(&zz);

    //Teuchos::TimeMonitor::summarize(std::cout);

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
