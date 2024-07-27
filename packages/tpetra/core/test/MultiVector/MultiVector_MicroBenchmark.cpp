// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "KokkosBlas.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::RCP;
  using std::endl;

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, MicroBenchmark_Update, SC, LO, GO, NT )
  {
    typedef Tpetra::Map<LO, GO, NT> map_type;
    typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;
    
    typename MV::impl_scalar_type ONE = STS::one();
    typename MV::impl_scalar_type ZERO = STS::zero();

    // Problem parameters
    size_t vector_size = 125000;
    size_t num_repeats = 2000;

    auto comm = getDefaultComm ();
    const int numProcs = comm->getSize ();
    if(numProcs != 1) return;

    const LO lclNumRows = vector_size;
    const GO indexBase = 0;
    RCP<const map_type> map (new map_type (lclNumRows, lclNumRows,
                                           indexBase, comm));
    const LO numVecs = 1;  // FIXME: This needs to be 1 for the lambda loop to work
    MV X (map, numVecs), Y(map,numVecs);
    X.putScalar(ZERO);
    Y.putScalar(ONE);

    Kokkos::fence();


    typedef typename MV::dual_view_type::t_dev dev_view_type;
    //    typedef typename MV::dual_view_type::t_host host_view_type;
    typedef typename dev_view_type::memory_space cur_memory_space;
    typedef typename dev_view_type::execution_space cur_exec_space;

    // Run Tpetra Update Loop
    {
      ::Tpetra::Details::ProfilingRegion region ("Tpetra Update Loop");
      for(size_t i=0; i<num_repeats; i++) {
        X.update(ONE,Y,ONE);
        Kokkos::fence();
      }
    }

    X.putScalar(ZERO);

    // Run KokkosBlas Update Loop
    {
      auto X_lcl =
           X.template getLocalView<cur_memory_space>(Tpetra::Access::ReadOnly);
      auto Y_lcl = 
           Y.template getLocalView<cur_memory_space>(Tpetra::Access::ReadWrite);

      ::Tpetra::Details::ProfilingRegion region ("KokkosBlas Update Loop");
      for(size_t i=0; i<num_repeats; i++) {
        KokkosBlas::axpby(ONE, X_lcl, ONE, Y_lcl);
        Kokkos::fence();
      }
    }

    X.putScalar(ZERO);
    Kokkos::fence();

    // Run raw lambda loop
    Kokkos::RangePolicy<cur_exec_space> policy (0, vector_size);
    {
      auto X_lcl =
           X.template getLocalView<cur_memory_space>(Tpetra::Access::ReadWrite);
      auto Y_lcl = 
           Y.template getLocalView<cur_memory_space>(Tpetra::Access::ReadOnly);
      ::Tpetra::Details::ProfilingRegion region ("Raw Lambda Update Loop");
      for(size_t rep=0; rep<num_repeats; rep++) {
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const size_t &i) {
            // This is what we should do in general
            //            X_lcl(i,0) = ONE * X_lcl(i,0) + ONE * Y_lcl(i,0);
            // This is the special case code that actually executes inside KokkosKernels for alpha=1
            X_lcl(i,0) += ONE * Y_lcl(i,0);
          });
        Kokkos::fence();
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, MicroBenchmark_Update, SC, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

