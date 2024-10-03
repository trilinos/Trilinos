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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MatvecFence, LO, GO, Scalar, Node )
  {

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
    FenceCounter::reset();
    FenceCounter::start();
    size_t iter_num = 10;

    {
      Tpetra::Details::ProfilingRegion r ("Matvec Loop");

      for (size_t i = 0; i < iter_num; i ++) {
        A1->apply(*X,*Y,Teuchos::NO_TRANS, alpha, beta);
        X->update(-1.0, *Y, 1.0);
        X->norm2(normX());
      }
    }
    FenceCounter::stop();

    auto exec_space = typename Node::execution_space();
    size_t expectedGlobalCount = 0;
    size_t expectedInstanceCount = 0;
    if (numRanks == 1) {
      expectedGlobalCount = 0;
      if (Node::is_cpu) {
        expectedInstanceCount = iter_num;
      } else {
        expectedInstanceCount = 2*iter_num;
      }
    } else {
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_SERIAL)
      // OpenMP in a Serial+OpenMP build
      if (Node::is_cpu && !Node::is_serial)
#else
      // OpenMP in an OpenMP build, Serial in a Serial or Serial+Cuda build
      if (Node::is_cpu)
#endif
      {
        if (Tpetra::Details::Behavior::assumeMpiIsGPUAware()) {
          expectedGlobalCount = iter_num;
        } else {
          expectedGlobalCount = 0;
        }
        if (Tpetra::Details::Behavior::debug()) {
          expectedInstanceCount = 3*iter_num;
        } else {
          expectedInstanceCount = 2*iter_num;
        }
      }
#ifdef HAVE_TPETRA_INST_OPENMP
      // Serial in a OpenMP+Serial build
      else if (Node::is_serial) {
        // Did not test the case of Serial node in build with Serial and OpenMP and GPU-aware
        expectedGlobalCount = iter_num;
        if (Tpetra::Details::Behavior::debug()) {
#if KOKKOS_VERSION >= 40499
          expectedInstanceCount = 3*iter_num;
#else
          expectedInstanceCount = 4*iter_num;
#endif
        } else {
#if KOKKOS_VERSION >= 40499
          expectedInstanceCount = 2*iter_num;
#else
          expectedInstanceCount = 3*iter_num;
#endif
        }
      }
#endif
      else {
        if (Tpetra::Details::Behavior::assumeMpiIsGPUAware()) {
          expectedGlobalCount = iter_num;
        } else {
          expectedGlobalCount = 6 * iter_num;
        }
#ifdef HAVE_TPETRA_INST_HIP
        if constexpr (std::is_same_v<typename Node::execution_space, Kokkos::HIP>) {
          if (Tpetra::Details::Behavior::debug()) {
            expectedInstanceCount = 4*iter_num;
          } else {
            expectedInstanceCount = 2*iter_num;
          }
        } else
#endif
        {
          if (Tpetra::Details::Behavior::debug()) {
            expectedInstanceCount = 5*iter_num;
          } else {
            expectedInstanceCount = 3*iter_num;
          }
        }
      }
    }

    TEST_EQUALITY(FenceCounter::get_count_global(exec_space.name()),   expectedGlobalCount);
    TEST_EQUALITY(FenceCounter::get_count_instance(exec_space.name()), expectedInstanceCount);

    stacked_timer->stopBaseTimer();
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(out, comm, options);

  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, MatvecFence, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
