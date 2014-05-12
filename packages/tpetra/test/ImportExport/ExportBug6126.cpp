/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_USE_CUDA_BUILD
  #define DO_COMPILATION
#else
  #ifndef KOKKOS_HAVE_CUDA
    #define DO_COMPILATION
  #endif
#endif
#else
  #define DO_COMPILATION
#endif

#ifdef DO_COMPILATION

#include "Teuchos_UnitTestHarness.hpp"

#include <Tpetra_TestingUtilities.hpp>

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <iterator>

namespace {
  using Tpetra::TestingUtilities::getNode;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::DefaultPlatform;
  using Tpetra::global_size_t;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::tuple;
  using Teuchos::Range1D;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::Export;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;

  using Tpetra::createContigMap;
  using Tpetra::createContigMapWithNode;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

 TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Export, Bug6126, Node )
  {
    typedef Tpetra::Map<int,int, Node> map_type;
    typedef Tpetra::Export<int,int, Node> export_type;
     RCP<const Comm<int> > comm = getDefaultComm();
    int NumProc = comm->getSize();
    int MyPID   = comm->getRank();
    
    // Build the overlapping row map
    int NumPrimaryEquations = 2;
    int NumLocalEquations = (NumProc>1)?(2*NumPrimaryEquations):NumPrimaryEquations;
    int start   = NumPrimaryEquations * MyPID;
    int g_total = NumPrimaryEquations * NumProc;
    Teuchos::Array<int> elementList(NumLocalEquations);
    for(int i=0; i<NumPrimaryEquations; i++) {
      elementList[i] = start + i;
      if(NumProc>1)
	elementList[i+NumPrimaryEquations] = (start + NumPrimaryEquations + i) % g_total;
    }
    RCP<const map_type > RowMap = rcp(new map_type(g_total, elementList(), 0, comm));
    
    // Create the not range map since this can't be allowed to overlap
    RCP<const map_type > RangeMap = rcp(new map_type(Teuchos::OrdinalTraits<int>::invalid(), NumPrimaryEquations, 0, comm));
    
    // Build exporter
    RCP<const export_type> Exporter = rcp(new export_type(RowMap,RangeMap));
    
    // Diagnostics on the maps
    printf("[%d] RangeMap  : ",MyPID);
    for(size_t i=0; i<RangeMap->getNodeNumElements(); i++) 
      printf("%d ",RangeMap->getGlobalElement(i));
    printf("\n");
    printf("[%d] RowMap    : ",MyPID);
    for(size_t i=0; i<RowMap->getNodeNumElements(); i++) 
      printf("%d ",RowMap->getGlobalElement(i));
    printf("\n");
    
    
    // Diagnostics on exporter
    RCP<const export_type> E = Exporter;
    printf("[%d] SameIDs   : ",MyPID);
    for(size_t i=0; i< E->getNumSameIDs(); i++)
      printf("%d ",(int) i);
    printf("\n");
    printf("[%d] PermuteIDs: ",MyPID);
    for(size_t i=0; i<E->getNumPermuteIDs(); i++)
      printf("(%d->%d) ",(int)E->getPermuteFromLIDs()[i],(int)E->getPermuteToLIDs()[i]);
    printf("\n");
    printf("[%d] RemoteIDs : ",MyPID);
    for(size_t i=0; i<E->getNumRemoteIDs(); i++)
      printf("%d ",(int)E->getRemoteLIDs()[i]);
    printf("\n");
    printf("[%d] ExportIDs : ",MyPID);
    for(size_t i=0; i<E->getNumExportIDs(); i++)
      printf("%d ",(int)E->getExportLIDs()[i]);
    printf("\n");
    fflush(stdout);


    if(!MyPID) printf("CMS: Note that neither processor thinks it has any remotes --- thus off-processor entries will be ignored\n");

    Tpetra::Vector<int,int,int,Node> v1(RowMap), v2(RangeMap);
    v1.putScalar(2);
    v2.putScalar(0);
    
    
    printf("[%d] v1   : ",MyPID);
    for(size_t i=0; i< RowMap->getNodeNumElements(); i++)
      printf("%d ",v1.get1dView()[i]);
    printf("\n");
    
    v2.doExport(v1,*Exporter,Tpetra::ADD);
    
    printf("[%d] v2   : ",MyPID);
    for(size_t i=0; i< RangeMap->getNodeNumElements(); i++)
      printf("%d ",v2.get1dView()[i]);
    printf("\n");
    
    if(!MyPID) printf("CMS: Note that the correct answer for v2 is supposed to be 4 in each entry, not 2\n");

    success=true;
    for(size_t i=0; i< RangeMap->getNodeNumElements(); i++)
      if(v2.get1dView()[i]!=4) success=false;

    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


#define UNIT_TEST_GROUP_NODE( NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Export, Bug6126, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
  TPETRA_INSTANTIATE_N( UNIT_TEST_GROUP_NODE )
#else
  TPETRA_INSTANTIATE_N_NOGPU( UNIT_TEST_GROUP_NODE )
#endif
}



#endif // DO_COMPILATION
