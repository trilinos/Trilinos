// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Apply_Helpers.hpp"
#include "TpetraUtils_MatrixGenerator.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_KokkosCounter.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <type_traits> // std::is_same

#include <Teuchos_StackedTimer.hpp>

namespace {

  // no ScalarTraits<>::eps() for integer types


  using std::endl;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createContigMapWithNode;


  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MatvecH2D, LO, GO, Scalar, Node )
  {

    // Test inspired by the reproduce in this issue: https://github.com/trilinos/Trilinos/issues/12560
    // So we can isolate H2D/D2H transfers in the matvec

    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;

    // This code is left in in case people want to debug future issues using the Kokkos profiling 
    // hooks in Tpetra
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("TransferPerf"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
    int numRanks = comm->getSize();

    const size_t nsize=3;

   /* Create the identity matrix, three rows per proc */
    RCP<MAT> A1 = Tpetra::Utils::MatrixGenerator<MAT>::generate_miniFE_matrix(nsize, comm);
    if(!A1->isFillComplete()) A1->fillComplete();
    auto map = A1->getRowMap();

    RCP<MV> X = rcp(new MV(map,1));  
    RCP<MV> Y = rcp(new MV(map,1));
    X->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    Y->putScalar(Teuchos::ScalarTraits<Scalar>::zero());

    Teuchos::Array<Mag> normX(1);

    const Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), beta = Teuchos::ScalarTraits<Scalar>::zero();

    {
      Tpetra::Details::ProfilingRegion r ("Sacrificial Matvec");
      A1->apply(*X,*Y,Teuchos::NO_TRANS, alpha, beta);
    }
    

    // Check to make sure we have the right number of H2D/D2H transfers
    using namespace Tpetra::Details;
    DeepCopyCounter::reset();
    DeepCopyCounter::start();
    size_t iter_num = 10;
      
    {
      Tpetra::Details::ProfilingRegion r ("Matvec Loop");
      
      for (size_t i = 0; i < iter_num; i ++) {
        A1->apply(*X,*Y,Teuchos::NO_TRANS, alpha, beta);
        X->update(-1.0, *Y, 1.0); 
        X->norm2(normX());
      }
    }
    DeepCopyCounter::stop();


    if (numRanks ==1 || Node::is_cpu || Tpetra::Details::Behavior::assumeMpiIsGPUAware()) {
      // No H2D/D2H if MPI is GPU aware, not a CPU, or if there's no MPI
      TEST_EQUALITY(DeepCopyCounter::get_count_different_space(),0);
    }
    else {
      // If there are multiple ranks w/o GPU-aware MPI, we currently have 1 copy per iter if the memory space
      // is a shared memory space, and 3 copies otherwise.
#ifdef KOKKOS_HAS_SHARED_SPACE
      constexpr size_t num_copies_per_iter = std::is_same_v<typename Node::memory_space, Kokkos::SharedSpace> ? 1 : 3;
#else
      constexpr size_t num_copies_per_iter = 3;
#endif
      TEST_EQUALITY(DeepCopyCounter::get_count_different_space(), num_copies_per_iter * iter_num);
    }
      

    stacked_timer->stopBaseTimer();
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(out, comm, options);

  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, MatvecH2D, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
