/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
