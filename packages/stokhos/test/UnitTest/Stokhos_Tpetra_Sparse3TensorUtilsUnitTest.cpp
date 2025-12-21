// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_Tpetra_Sparse3TensorUtilities.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace TpetraSparse3TensorUnitTest {

  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol, sparse_tol;
    OrdinalType d, p, sz;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases;
    Teuchos::RCP<Stokhos::OrthogPolyBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<Cijk_type> Cijk;
    
    UnitTestSetup(): d(3), p(4), bases(d) {
      rtol = 1e-14;
      atol = 1e-14;
      sparse_tol = 1e-14;
      
      // Create product basis
      for (OrdinalType i=0; i<d; i++)
        bases[i] = 
          Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
      basis =
	      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases, sparse_tol));
      sz = basis->size();

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
    }
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_TpetraSparse3Tensor, Graph ) {
    using TpetraLocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type;
    using TpetraGlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type;
    using TpetraNode = Tpetra::Details::DefaultTypes::node_type;
    using TpetraGraph = Tpetra::CrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>;
    using TpetraComm = Teuchos::Comm<int>;
    using host_view_type = TpetraGraph::nonconst_local_inds_host_view_type;

    // Create MPI comm
    Teuchos::RCP<const TpetraComm> comm = Tpetra::getDefaultComm();

    // Create CrsGraph from Cijk tensor
    Teuchos::RCP<TpetraGraph> graph = 
      Stokhos::sparse3Tensor2TpetraCrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>(
        *setup.basis, *setup.Cijk, comm);

    // Check all entries in Cijk are in the graph
    success = true;
    for (Cijk_type::k_iterator k_it=setup.Cijk->k_begin(); k_it!=setup.Cijk->k_end(); ++k_it) {
      for (Cijk_type::kj_iterator j_it = setup.Cijk->j_begin(k_it); j_it != setup.Cijk->j_end(k_it); ++j_it) {
        TpetraLocalOrdinal j = Stokhos::index(j_it);
        for (Cijk_type::kji_iterator i_it = setup.Cijk->i_begin(j_it); i_it != setup.Cijk->i_end(j_it); ++i_it) {
          TpetraLocalOrdinal i = Stokhos::index(i_it);
          size_t num_col_inds = graph->getNumAllocatedEntriesInLocalRow(i);
          host_view_type col_inds("col_inds", num_col_inds);
          graph->getLocalRowCopy(i, col_inds, num_col_inds);
          bool my_success = false;
          for (size_t l=0; l<num_col_inds; ++l) {
            if (col_inds(l) == j) {
              my_success = true;
              break;
            }
          }
          if (!my_success) {
            out << std::endl << "For row index " << i << ", could not find column index " << j << std::endl;
            out << "column indices: ";
            for (size_t l=0; l<num_col_inds; ++l) {
              out << col_inds(l) << " ";
            }
            out << std::endl;
          }
          success = success && my_success;
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_TpetraSparse3Tensor, MatrixMarketWriter ) {
    using TpetraLocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type;
    using TpetraGlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type;
    using TpetraNode = Tpetra::Details::DefaultTypes::node_type;
    using TpetraGraph = Tpetra::CrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>;
    using TpetraMatrix = Tpetra::CrsMatrix<double,TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>;
    using TpetraReader = Tpetra::MatrixMarket::Reader<TpetraMatrix>;
    using TpetraComm = Teuchos::Comm<int>;
    using host_view_type = TpetraGraph::nonconst_local_inds_host_view_type;

    // Create MPI comm
    Teuchos::RCP<const TpetraComm> comm = Tpetra::getDefaultComm();

    // Write Cijk to matrix market file
    const std::string mm_file = "Cijk.mm";
    Stokhos::sparse3Tensor2TpetraMatrixMarket<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>(
      *setup.basis, *setup.Cijk, comm, mm_file);

    // Read graph from file
    Teuchos::RCP<TpetraGraph> graph = TpetraReader::readSparseGraphFile(mm_file, comm);

    // Check all entries in Cijk are in the graph
    success = true;
    for (Cijk_type::k_iterator k_it=setup.Cijk->k_begin(); k_it!=setup.Cijk->k_end(); ++k_it) {
      for (Cijk_type::kj_iterator j_it = setup.Cijk->j_begin(k_it); j_it != setup.Cijk->j_end(k_it); ++j_it) {
        TpetraLocalOrdinal j = Stokhos::index(j_it);
        for (Cijk_type::kji_iterator i_it = setup.Cijk->i_begin(j_it); i_it != setup.Cijk->i_end(j_it); ++i_it) {
          TpetraLocalOrdinal i = Stokhos::index(i_it);
          size_t num_col_inds = graph->getNumAllocatedEntriesInLocalRow(i);
          host_view_type col_inds("col_inds", num_col_inds);
          graph->getLocalRowCopy(i, col_inds, num_col_inds);
          bool my_success = false;
          for (size_t l=0; l<num_col_inds; ++l) {
            if (col_inds(l) == j) {
              my_success = true;
              break;
            }
          }
          if (!my_success) {
            out << std::endl << "For row index " << i << ", could not find column index " << j << std::endl;
            out << "column indices: ";
            for (size_t l=0; l<num_col_inds; ++l) {
              out << col_inds(l) << " ";
            }
            out << std::endl;
          }
          success = success && my_success;
        }
      }
    }
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
