// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestTemplate.cpp

\brief Ifpack2 Unit testing template.

This file demonstrates how you create a unit test for template code.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_IlukGraph.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2IlukGraph, IlukGraphTest0, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef typename crs_graph_type::local_graph_device_type local_graph_type;
  typedef typename local_graph_type::row_map_type lno_row_view_t;
  typedef typename local_graph_type::entries_type lno_nonzero_view_t;
  typedef typename local_graph_type::device_type::memory_space TemporaryMemorySpace;
  typedef typename local_graph_type::device_type::memory_space PersistentMemorySpace;
  typedef typename local_graph_type::device_type::execution_space HandleExecSpace;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle
      <typename lno_row_view_t::const_value_type, typename lno_nonzero_view_t::const_value_type, double,
       HandleExecSpace, TemporaryMemorySpace,PersistentMemorySpace > kk_handle_type;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    tif_utest::create_tridiag_graph<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  int num_procs = crsgraph->getRowMap()->getComm()->getSize();
  TEST_EQUALITY( crsgraph->getRowMap()->getLocalNumElements(), num_rows_per_proc)
  TEST_EQUALITY( crsgraph->getRowMap()->getGlobalNumElements(), num_rows_per_proc*num_procs)

  TEST_EQUALITY( crsgraph->getGlobalNumRows(),crsgraph->getRowMap()->getGlobalNumElements())
  
  LocalOrdinal overlap_levels = 2;

  //Original Serial implementation
  {
    LocalOrdinal fill_levels = 0;

    Ifpack2::IlukGraph<crs_graph_type, kk_handle_type> iluk0_graph(crsgraph, fill_levels, overlap_levels);
    iluk0_graph.initialize();

    //The number of nonzeros in an ILU(0) graph should be the same as the
    //number of nonzeros in the input graph:

    size_t nnz0 = iluk0_graph.getL_Graph()->getGlobalNumEntries() +
                  iluk0_graph.getU_Graph()->getGlobalNumEntries() +
                  iluk0_graph.getNumGlobalDiagonals();

    fill_levels = 2;

    Ifpack2::IlukGraph<crs_graph_type, kk_handle_type> iluk2_graph(crsgraph, fill_levels, overlap_levels);
    iluk2_graph.initialize();

    //The number of nonzeros in an ILU(2) graph should be greater than the
    //number of nonzeros in the ILU(0) graph:

    size_t nnz2 = iluk2_graph.getL_Graph()->getGlobalNumEntries() +
                  iluk2_graph.getU_Graph()->getGlobalNumEntries() +
                  iluk2_graph.getNumGlobalDiagonals();

    TEST_EQUALITY( nnz2 >= nnz0, true)
  }

  //Kokkos Kernels KSPILUK implementation
  {
    Teuchos::RCP<kk_handle_type> KernelHandle0, KernelHandle2;

    LocalOrdinal fill_levels = 0;

    KernelHandle0 = Teuchos::rcp (new kk_handle_type ());
    KernelHandle0->create_spiluk_handle( KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1, 
                                         crsgraph->getLocalNumRows(),
                                         2*crsgraph->getLocalNumEntries()*(fill_levels+1), 
                                         2*crsgraph->getLocalNumEntries()*(fill_levels+1) );

    Ifpack2::IlukGraph<crs_graph_type, kk_handle_type> iluk0_graph(crsgraph, fill_levels, overlap_levels);
    iluk0_graph.initialize(KernelHandle0);

    //The number of nonzeros in an ILU(0) graph should be the same as the
    //number of nonzeros in the input graph:

    size_t nnz0 = iluk0_graph.getL_Graph()->getGlobalNumEntries() +
                  iluk0_graph.getU_Graph()->getGlobalNumEntries() -
                  iluk0_graph.getL_Graph()->getGlobalNumRows();

    fill_levels = 2;

    KernelHandle2 = Teuchos::rcp (new kk_handle_type ());
    KernelHandle2->create_spiluk_handle( KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1, 
                                         crsgraph->getLocalNumRows(),
                                         2*crsgraph->getLocalNumEntries()*(fill_levels+1), 
                                         2*crsgraph->getLocalNumEntries()*(fill_levels+1) );


    Ifpack2::IlukGraph<crs_graph_type, kk_handle_type> iluk2_graph(crsgraph, fill_levels, overlap_levels);
    iluk2_graph.initialize(KernelHandle2);

    //The number of nonzeros in an ILU(2) graph should be greater than the
    //number of nonzeros in the ILU(0) graph:

    size_t nnz2 = iluk2_graph.getL_Graph()->getGlobalNumEntries() +
                  iluk2_graph.getU_Graph()->getGlobalNumEntries() -
                  iluk2_graph.getL_Graph()->getGlobalNumRows();

    TEST_EQUALITY( nnz2 >= nnz0, true)
  }
}

#define UNIT_TEST_GROUP_LO_GO(LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2IlukGraph, IlukGraphTest0, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

IFPACK2_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

} // namespace (anonymous)
