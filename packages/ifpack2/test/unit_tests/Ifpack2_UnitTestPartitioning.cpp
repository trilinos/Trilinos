// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Ifpack2_OverlappingPartitioner.hpp>
#include <Ifpack2_LinearPartitioner.hpp>





using Tpetra::global_size_t;
typedef tif_utest::Node Node;
using namespace std;
using Teuchos::rcp;
using Teuchos::RCP;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Partitioning, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> CRSG;

  // Useful matrices and such (tridiagonal test)
  global_size_t num_rows_per_proc = 5;
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  // ====================================== //
  //            point blocking              //
  // ====================================== //
  {
    Teuchos::ParameterList List;
    List.set("partitioner: local parts",(int) num_rows_per_proc);
    Ifpack2::LinearPartitioner<CRSG > MyPart(Matrix->getGraph());
    MyPart.setParameters(List);
    MyPart.compute();
    const Teuchos::ArrayView<const LocalOrdinal>  & myview = MyPart.nonOverlappingPartition();

    Teuchos::Array<LocalOrdinal> correct_solution((int)num_rows_per_proc);
    for(int i=0;i<(int)num_rows_per_proc;i++)
      correct_solution[i]=i;

    TEST_COMPARE_ARRAYS(myview,correct_solution);
  }

  // ====================================== //
  //            full blocking               //
  // ====================================== //
  // Point blocking
  {
    Teuchos::ParameterList List;
    List.set("partitioner: local parts",(int) 1);
    Ifpack2::LinearPartitioner<CRSG > MyPart(Matrix->getGraph());
    MyPart.setParameters(List);
    MyPart.compute();
    const Teuchos::ArrayView<const LocalOrdinal>  & myview = MyPart.nonOverlappingPartition();

    Teuchos::Array<LocalOrdinal> correct_solution((int)num_rows_per_proc);
    for(int i=0;i<(int)num_rows_per_proc;i++)
      correct_solution[i]=0;

    TEST_COMPARE_ARRAYS(myview,correct_solution);
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Partitioning, Test0, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types.

IFPACK2_INSTANTIATE_SLG( UNIT_TEST_GROUP_SC_LO_GO )
