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
#include "Tpetra_Details_gemm.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_BLAS.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include <typeinfo>

namespace {
  using std::endl;
  typedef int LO;

  template<class ValueType,
           const bool isInteger = std::is_integral<ValueType>::value>
  struct MachinePrecision
  {
    typedef typename Kokkos::ArithTraits<ValueType>::mag_type mag_type;

    static mag_type machinePrecision ();
  };

  template<class ValueType>
  struct MachinePrecision<ValueType, true>
  {
    typedef typename Kokkos::ArithTraits<ValueType>::mag_type mag_type;

    static mag_type machinePrecision () { return 0; }
  };

  template<class ValueType>
  struct MachinePrecision<ValueType, false>
  {
    typedef typename Kokkos::ArithTraits<ValueType>::mag_type mag_type;

    static mag_type machinePrecision () {
      return Kokkos::ArithTraits<ValueType>::eps ();
    }
  };

  template<class ValueType>
  typename MachinePrecision<ValueType>::mag_type
  machinePrecision ()
  {
    return MachinePrecision<ValueType>::machinePrecision ();
  }

  template<class DeviceType>
  using PseudorandomPoolType = Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space>;

  template<class DeviceType>
  PseudorandomPoolType<DeviceType>
  preparePseudorandomNumberGenerator ()
  {
    // Seed the pseudorandom number generator using the same seed each
    // time.  This test doesn't care about MPI processes, so we don't
    // need to decorrelate different process' pseudorandom streams.
    // The point is really just to get some input numbers for tests.
    constexpr unsigned int seed = 42;
    return PseudorandomPoolType<DeviceType> (seed);
  }

  template<class EntryType, class CoeffType, class DeviceType>
  void
  testGemmVsTeuchosBlasForOneTransComb (Teuchos::FancyOStream& out,
                                        bool& success,
                                        PseudorandomPoolType<DeviceType>& randPool,
                                        const LO m,
                                        const LO n,
                                        const LO k,
                                        const char trans_A,
                                        const bool trans_A_is_trans,
                                        const bool trans_A_is_conj,
                                        const char trans_B,
                                        const bool trans_B_is_trans,
                                        const bool trans_B_is_conj)
  {
    // Convert types like std::complex into their equivalent
    // Kokkos-friendly types.
    typedef typename Kokkos::ArithTraits<EntryType>::val_type entry_type;
    typedef typename Kokkos::ArithTraits<CoeffType>::val_type coeff_type;
    // Teuchos::BLAS only works with LayoutLeft matrices.
    typedef Kokkos::LayoutLeft layout_type;
    typedef Kokkos::View<entry_type**, layout_type, DeviceType> mat_type;
    typedef typename mat_type::HostMirror::execution_space host_exec_space;
    typedef Kokkos::Device<host_exec_space, Kokkos::HostSpace> host_dev_type;
    typedef Kokkos::View<entry_type**, layout_type, host_dev_type> host_mat_type;
    typedef PseudorandomPoolType<DeviceType> pool_type;
    typedef typename pool_type::generator_type generator_type;
    typedef typename Kokkos::ArithTraits<entry_type>::mag_type mag_type;

    Teuchos::OSTab tab0 (out);
    out << "Test instance:" << endl;
    Teuchos::OSTab tab1 (out);
    out << "entry_type: "
        << typeid (entry_type).name () << endl
        << "coeff_type: "
        << typeid (coeff_type).name () << endl
        << "execution_space: "
        << typeid (typename DeviceType::execution_space).name () << endl
        << "memory_space: "
        << typeid (typename DeviceType::memory_space).name () << endl
        << "m: " << m << endl
        << "n: " << n << endl
        << "k: " << k << endl
        << "trans_A: " << trans_A << endl
        << "trans_B: " << trans_B << endl;

    mat_type A_orig ("A_orig",
                     trans_A_is_trans ? k : m,
                     trans_A_is_trans ? m : k);
    mat_type B_orig ("B_orig",
                     trans_B_is_trans ? n : k,
                     trans_B_is_trans ? k : n);
    mat_type C_orig ("C_orig", m, n);

    {
      typedef Kokkos::ArithTraits<entry_type> KATE;
      const entry_type maxVal =
        Kokkos::rand<generator_type, entry_type>::max ();
      const entry_type minVal =
        KATE::is_signed ? entry_type (-maxVal) : KATE::zero ();
      Kokkos::fill_random (A_orig, randPool, minVal, maxVal);
      Kokkos::fill_random (B_orig, randPool, minVal, maxVal);
      Kokkos::fill_random (C_orig, randPool, minVal, maxVal);
    }

    mat_type A ("A", A_orig.dimension_0 (), A_orig.dimension_1 ());
    mat_type B ("B", B_orig.dimension_0 (), B_orig.dimension_1 ());
    mat_type C ("C", C_orig.dimension_0 (), C_orig.dimension_1 ());
    auto A_host = Kokkos::create_mirror_view (A);
    auto B_host = Kokkos::create_mirror_view (B);
    auto C_host = Kokkos::create_mirror_view (C);

    // It's important that these be separate allocations from A, B,
    // and C, in case the routine to test is buggy and modifies the
    // wrong thing.  The host mirror of a View may be the View itself,
    // thus the same allocation.
    host_mat_type A2 ("A2", A.dimension_0 (), A.dimension_1 ());
    host_mat_type B2 ("B2", B.dimension_0 (), B.dimension_1 ());
    host_mat_type C2 ("C2", C.dimension_0 (), C.dimension_1 ());

    // Use the host versions of A, B, and C to compute max norms.
    // We'll need these for relative error bounds.
    mag_type A_norm = Kokkos::ArithTraits<mag_type>::zero ();
    mag_type B_norm = Kokkos::ArithTraits<mag_type>::zero ();
    mag_type C_norm = Kokkos::ArithTraits<mag_type>::zero ();
    {
      Kokkos::deep_copy (A2, A_orig);
      for (LO i = 0; i < static_cast<LO> (A2.dimension_0 ()); ++i) {
        for (LO j = 0; j < static_cast<LO> (A2.dimension_1 ()); ++j) {
          const mag_type curAbs =
            Kokkos::ArithTraits<entry_type>::abs (A2(i,j));
          A_norm = (curAbs > A_norm) ? curAbs : A_norm;
        }
      }
      Kokkos::deep_copy (B2, B_orig);
      for (LO i = 0; i < static_cast<LO> (B2.dimension_0 ()); ++i) {
        for (LO j = 0; j < static_cast<LO> (B2.dimension_1 ()); ++j) {
          const mag_type curAbs =
            Kokkos::ArithTraits<entry_type>::abs (B2(i,j));
          B_norm = (curAbs > B_norm) ? curAbs : B_norm;
        }
      }
      Kokkos::deep_copy (C2, A_orig);
      for (LO i = 0; i < static_cast<LO> (C2.dimension_0 ()); ++i) {
        for (LO j = 0; j < static_cast<LO> (C2.dimension_1 ()); ++j) {
          const mag_type curAbs =
            Kokkos::ArithTraits<entry_type>::abs (C2(i,j));
          C_norm = (curAbs > C_norm) ? curAbs : C_norm;
        }
      }
    }

    const LO A2_stride = Tpetra::Details::Blas::getStride2DView (A2);
    const LO B2_stride = Tpetra::Details::Blas::getStride2DView (B2);
    const LO C2_stride = Tpetra::Details::Blas::getStride2DView (C2);

    typedef Kokkos::ArithTraits<coeff_type> KAT;
    const coeff_type zero = KAT::zero ();
    const coeff_type one = KAT::one ();
    const coeff_type two = one + one;
    constexpr int numCoeffOpts = 5;
    // coeff_type may not have constexpr constructors, so we can't
    // necessarily use a constexpr array here.
    const coeff_type alphaOpts[numCoeffOpts] = {zero, one, -one, two, -two};
    const coeff_type betaOpts[numCoeffOpts] = {zero, one, -one, two, -two};

    const mag_type eps = machinePrecision<entry_type> ();

    // Inf-norm error bound depends on the number of terms in the sum.
    out << "A_norm: " << A_norm << endl
        << "B_norm: " << B_norm << endl
        << "C_norm: " << C_norm << endl;
    // Add a little "fudge factor."  2 is enough for real, and 4 is
    // enough for complex.
    const mag_type fudgeFactor = Kokkos::ArithTraits<entry_type>::is_complex ?
      static_cast<mag_type> (8) :
      static_cast<mag_type> (4);

    const mag_type boundOrig =
      eps * fudgeFactor * (A_norm * B_norm * static_cast<mag_type> (k));
    const mag_type bound = boundOrig < (fudgeFactor*eps) ? (fudgeFactor*eps) : boundOrig;
    out << "bound: " << bound << endl;

    // This needs to be EntryType and not
    // ArithTraits<EntryType>::val_type, because the BLAS is
    // specialized for std::complex and not Kokkos::complex.
    Teuchos::BLAS<int, EntryType> teuchosBlas;

    for (coeff_type alpha : alphaOpts) {
      for (coeff_type beta : betaOpts) {
        // Refresh test problem from original version.
        Kokkos::deep_copy (A, A_orig);
        Kokkos::deep_copy (B, B_orig);
        Kokkos::deep_copy (C, C_orig);
        Kokkos::deep_copy (A2, A_orig);
        Kokkos::deep_copy (B2, B_orig);
        Kokkos::deep_copy (C2, C_orig);

        {
          ::Tpetra::Details::Blas::gemm (trans_A, trans_B, alpha, A, B, beta, C);
          const Teuchos::ETransp teuchosTransA = trans_A_is_trans ?
            (trans_A_is_conj ? Teuchos::CONJ_TRANS : Teuchos::TRANS) :
            Teuchos::NO_TRANS;
          const Teuchos::ETransp teuchosTransB = trans_B_is_trans ?
            (trans_B_is_conj ? Teuchos::CONJ_TRANS : Teuchos::TRANS) :
            Teuchos::NO_TRANS;
          teuchosBlas.GEMM (teuchosTransA, teuchosTransB, m, n, k, alpha,
                            const_cast<const EntryType*> (reinterpret_cast<EntryType*> (A2.data ())),
                            A2_stride,
                            const_cast<const EntryType*> (reinterpret_cast<EntryType*> (B2.data ())),
                            B2_stride,
                            beta,
                            reinterpret_cast<EntryType*> (C2.data ()),
                            C2_stride);
          Kokkos::deep_copy (C_host, C);
        }

        {
          mag_type maxErr = Kokkos::ArithTraits<mag_type>::zero ();
          for (LO i = 0; i < static_cast<LO> (C_host.dimension_0 ()); ++i) {
            for (LO j = 0; j < static_cast<LO> (C_host.dimension_1 ()); ++j) {
              const mag_type curErr =
                Kokkos::ArithTraits<entry_type>::abs (C_host(i,j) - C2(i,j));
              maxErr = (curErr > maxErr) ? curErr : maxErr;
            }
          }
          const bool wrong = (maxErr > bound);
          if (wrong) {
            TEST_ASSERT( false );
            out << "Comparing Tpetra::Details::Blas::gemm against "
              "Teuchos::BLAS::GEMM: Max error " << maxErr << " > error bound "
                << bound << endl;
          }
        }

        // Repeat test for the "default" Tpetra::Details::Blas::gemm implementation.
        {
          Kokkos::deep_copy (A, A_orig);
          Kokkos::deep_copy (B, B_orig);
          Kokkos::deep_copy (C, C_orig);
          ::Tpetra::Details::Blas::Default::gemm (trans_A, trans_B, alpha,
                                                  A, B, beta, C);
          Kokkos::deep_copy (C_host, C);
        }

        {
          mag_type maxErr = Kokkos::ArithTraits<mag_type>::zero ();
          for (LO i = 0; i < static_cast<LO> (C_host.dimension_0 ()); ++i) {
            for (LO j = 0; j < static_cast<LO> (C_host.dimension_1 ()); ++j) {
              const mag_type curErr =
                Kokkos::ArithTraits<entry_type>::abs (C_host(i,j) - C2(i,j));
              maxErr = (curErr > maxErr) ? curErr : maxErr;
            }
          }
          const bool wrong = (maxErr > bound);
          if (wrong) {
            TEST_ASSERT( false );
            out << "Comparing Tpetra::Details::Blas::Default::gemm against "
              "Teuchos::BLAS::GEMM: Max error " << maxErr << " > error bound "
                << bound << endl;
          }
        }

        if (! success) {
          out << "At least one test FAILED; abandoning the others." << endl;
          return;
        }
      } // beta
    } // alpha
  }


  template<class EntryType, class CoeffType, class DeviceType>
  void
  testGemmVsTeuchosBlas (Teuchos::FancyOStream& out,
                         bool& success,
                         PseudorandomPoolType<DeviceType>& randPool,
                         const LO m,
                         const LO n,
                         const LO k)
  {
    constexpr int numTransOpts = 6;
    constexpr char transOpts[numTransOpts] =
      {'N',   'n',   'T',   't',   'C',  'c'};
    constexpr bool transOptIsTrans[numTransOpts] =
      {false, false, true,  true,  true, true};
    constexpr bool transOptIsConj[numTransOpts] =
      {false, false, false, false, true, true};

    for (int transInd_A = 0; transInd_A < numTransOpts; ++transInd_A) {
      const char trans_A = transOpts[transInd_A];
      const bool isTrans_A = transOptIsTrans[transInd_A];
      const bool isConj_A = transOptIsConj[transInd_A];

      for (int transInd_B = 0; transInd_B < numTransOpts; ++transInd_B) {
        const char trans_B = transOpts[transInd_B];
        const bool isTrans_B = transOptIsTrans[transInd_B];
        const bool isConj_B = transOptIsConj[transInd_B];

        testGemmVsTeuchosBlasForOneTransComb<EntryType, CoeffType, DeviceType>
          (out, success, randPool, m, n, k,
           trans_A, isTrans_A, isConj_A,
           trans_B, isTrans_B, isConj_B);
        if (! success) {
          out << "At least one test FAILED; abandoning the others." << endl;
          return;
        }
      }
    }
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Blas, Gemm, SCALAR )
  {
    typedef SCALAR entry_type;
    typedef SCALAR coeff_type;
    typedef Tpetra::Map<> map_type;
    typedef map_type::device_type device_type;

    Teuchos::OSTab tab0 (out);
    out << "Test \"Tpetra::Details::Blas::gemm\"" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    // Creating a Map instance takes care of Kokkos initialization and
    // finalization automatically.
    Tpetra::Map<> map (comm->getSize (), 1, 0, comm);

    auto randPool = preparePseudorandomNumberGenerator<device_type> ();
    const LO m_vals[] = {1, 2, 5, 13};
    const LO n_vals[] = {1, 2, 5, 13};
    const LO k_vals[] = {1, 2, 5, 13};
    for (LO m : m_vals) {
      for (LO n : n_vals) {
        for (LO k : k_vals) {
          testGemmVsTeuchosBlas<entry_type, coeff_type, device_type> (out,
                                                                      success,
                                                                      randPool,
                                                                      m, n, k);
          if (! success) {
            out << "At least one test FAILED; abandoning the others." << endl;
            return;
          }
        }
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Blas, Gemm, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_GROUP )

} // namespace (anonymous)


