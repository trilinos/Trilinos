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
#include "Teuchos_BLAS.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

// kokkos kernels
#include "KokkosBlas.hpp"
#include "KokkosKernels_Utils.hpp"

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
  testGemvVsTeuchosBlas (Teuchos::FancyOStream& out,
                         bool& success,
                         PseudorandomPoolType<DeviceType>& randPool,
                         const LO numRows,
                         const LO numCols)
  {
    // Convert types like std::complex into their equivalent
    // Kokkos-friendly types.
    typedef typename Kokkos::ArithTraits<EntryType>::val_type entry_type;
    typedef typename Kokkos::ArithTraits<CoeffType>::val_type coeff_type;
    // Teuchos::BLAS only works with LayoutLeft matrices.
    typedef Kokkos::LayoutLeft layout_type;
    typedef Kokkos::View<entry_type**, layout_type, DeviceType> mat_type;
    typedef Kokkos::View<entry_type*, layout_type, DeviceType> vec_type;
    typedef typename mat_type::HostMirror::execution_space host_exec_space;
    typedef Kokkos::Device<host_exec_space, Kokkos::HostSpace> host_dev_type;
    typedef Kokkos::View<entry_type**, layout_type, host_dev_type> host_mat_type;
    typedef Kokkos::View<entry_type*, layout_type, host_dev_type> host_vec_type;
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
        << "numRows: " << numRows << endl
        << "numCols: " << numCols << endl;

    mat_type A_orig ("A_orig", numRows, numCols);
    vec_type x_orig ("x_orig", numCols);
    vec_type y_orig ("y_orig", numRows);

    {
      typedef Kokkos::ArithTraits<entry_type> KATE;
      const entry_type maxVal =
        Kokkos::rand<generator_type, entry_type>::max ();
      const entry_type minVal =
        KATE::is_signed ? entry_type (-maxVal) : KATE::zero ();
      Kokkos::fill_random (A_orig, randPool, minVal, maxVal);
      Kokkos::fill_random (x_orig, randPool, minVal, maxVal);
      Kokkos::fill_random (y_orig, randPool, minVal, maxVal);
    }

    mat_type A ("A", numRows, numCols);
    vec_type x ("x", numCols);
    vec_type y ("y", numRows);
    auto x_host = Kokkos::create_mirror_view (x);
    auto y_host = Kokkos::create_mirror_view (y);

    // It's important that these be separate allocations from A, x,
    // and y.  The host mirror of a View may be the View itself, thus
    // the same allocation.
    host_mat_type A2 ("A2", numRows, numCols);
    host_vec_type x2 ("x2", numCols);
    host_vec_type y2 ("y2", numRows);

    // Use the host versions of A, x, and y to compute max norms.
    // We'll need these for relative error bounds.
    mag_type A_norm = Kokkos::ArithTraits<mag_type>::zero ();
    mag_type x_norm = Kokkos::ArithTraits<mag_type>::zero ();
    mag_type y_norm = Kokkos::ArithTraits<mag_type>::zero ();
    {
      Kokkos::deep_copy (A2, A_orig);
      for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
          const mag_type curAbs =
            Kokkos::ArithTraits<entry_type>::abs (A2(i,j));
          A_norm = (curAbs > A_norm) ? curAbs : A_norm;
        }
      }
      Kokkos::deep_copy (x2, x_orig);
      for (int i = 0; i < static_cast<int> (x2.extent (0)); ++i) {
        const mag_type curAbs =
          Kokkos::ArithTraits<entry_type>::abs (x2(i));
        x_norm = (curAbs > x_norm) ? curAbs : x_norm;
      }
      Kokkos::deep_copy (y2, y_orig);
      for (int i = 0; i < static_cast<int> (y2.extent (0)); ++i) {
        const mag_type curAbs =
          Kokkos::ArithTraits<entry_type>::abs (y2(i));
        y_norm = (curAbs > y_norm) ? curAbs : y_norm;
      }
    }

    const int A2_stride = A2.stride(1);
    const int x2_inc = x2.stride(0);
    const int y2_inc = y2.stride(0);

    constexpr int numTransOpts = 6;
    constexpr char transOpts[numTransOpts] =
      {'N',   'n',   'T',   't',   'C',  'c'};
    constexpr bool transOptIsTrans[numTransOpts] =
      {false, false, true,  true,  true, true};
    constexpr bool transOptIsConj[numTransOpts] =
      {false, false, false, false, true, true};

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
        << "x_norm: " << x_norm << endl
        << "y_norm: " << y_norm << endl;
    // Add a little "fudge factor."  2 is enough for real, and 4 is
    // enough for complex.
    const mag_type fudgeFactor = Kokkos::ArithTraits<entry_type>::is_complex ?
      static_cast<mag_type> (4) :
      static_cast<mag_type> (2);
    const mag_type nonTransFactor = A_norm * x_norm > fudgeFactor ?
      A_norm * x_norm :
      fudgeFactor;
    const mag_type nonTransBound = nonTransFactor *
      static_cast<mag_type> (numCols) * eps;
    const mag_type transFactor = A_norm * y_norm > fudgeFactor ?
      A_norm * y_norm :
      fudgeFactor;
    const mag_type transBound = transFactor *
      static_cast<mag_type> (numRows) * eps;
    out << "nonTransBound: " << nonTransBound << endl
        << "transBound: " << transBound << endl;

    // This needs to be EntryType and not
    // ArithTraits<EntryType>::val_type, because the BLAS is
    // specialized for std::complex and not Kokkos::complex.
    Teuchos::BLAS<int, EntryType> teuchosBlas;

    for (coeff_type alpha : alphaOpts) {
      for (coeff_type beta : betaOpts) {
        for (int transInd = 0; transInd < numTransOpts; ++transInd) {
          const char trans = transOpts[transInd];
          const bool isTrans = transOptIsTrans[transInd];
          const bool isConj = transOptIsConj[transInd];

          // Refresh test problem from original version.
          Kokkos::deep_copy (A, A_orig);
          Kokkos::deep_copy (x, x_orig);
          Kokkos::deep_copy (y, y_orig);

          Kokkos::deep_copy (A2, A_orig);
          Kokkos::deep_copy (x2, x_orig);
          Kokkos::deep_copy (y2, y_orig);

          const Teuchos::ETransp teuchosTrans = isTrans ?
            (isConj ? Teuchos::CONJ_TRANS : Teuchos::TRANS) :
            Teuchos::NO_TRANS;

          // Whether (x,y) are (input,output) or vice versa depends on
          // whether we are exercising the transpose.
          if (isTrans) {
            KokkosBlas::gemv(&trans, alpha, A, y, beta, x); 
            teuchosBlas.GEMV (teuchosTrans, numRows, numCols, alpha,
                              reinterpret_cast<const EntryType*> (A2.data ()),
                              A2_stride,
                              reinterpret_cast<const EntryType*> (y2.data ()),
                              y2_inc, beta,
                              reinterpret_cast<EntryType*> (x2.data ()),
                              x2_inc);
            Kokkos::deep_copy (x_host, x);

            mag_type maxErr = Kokkos::ArithTraits<mag_type>::zero ();
            for (LO i = 0; i < static_cast<LO> (x.extent (0)); ++i) {
              const mag_type curErr =
                Kokkos::ArithTraits<entry_type>::abs (x_host(i) - x2(i));
              maxErr = (curErr > maxErr) ? curErr : maxErr;
            }
            const bool wrong = (maxErr > transBound);
            TEST_ASSERT( ! wrong );
            if (wrong) {
              out << "Max error " << maxErr
                  << " > error bound " << transBound << endl;
            }
          }
          else {
            KokkosBlas::gemv(&trans, alpha, A, x, beta, y); 
            teuchosBlas.GEMV (teuchosTrans, numRows, numCols, alpha,
                              reinterpret_cast<const EntryType*> (A2.data ()),
                              A2_stride,
                              reinterpret_cast<const EntryType*> (x2.data ()),
                              x2_inc, beta,
                              reinterpret_cast<EntryType*> (y2.data ()),
                              y2_inc);
            Kokkos::deep_copy (y_host, y);

            mag_type maxErr = Kokkos::ArithTraits<mag_type>::zero ();
            for (LO i = 0; i < static_cast<LO> (y.extent (0)); ++i) {
              const mag_type curErr =
                Kokkos::ArithTraits<entry_type>::abs (y_host(i) - y2(i));
              maxErr = (curErr > maxErr) ? curErr : maxErr;
            }
            const bool wrong = (maxErr > nonTransBound);
            TEST_ASSERT( ! wrong );
            if (wrong) {
              out << "Max error " << maxErr
                  << " > error bound " << nonTransBound << endl;
            }
          }
        }
      }
    }
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Blas, Gemv, SCALAR )
  {
    typedef SCALAR entry_type;
    typedef SCALAR coeff_type;
    typedef Tpetra::Map<> map_type;
    typedef map_type::device_type device_type;

    Teuchos::OSTab tab0 (out);
    out << "Test \"KokkosBlas::gemv\"" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    // Creating a Map instance takes care of Kokkos initialization and
    // finalization automatically.
    Tpetra::Map<> map (comm->getSize (), 1, 0, comm);

    auto randPool = preparePseudorandomNumberGenerator<device_type> ();
    //const LO numRowsVals[] = {0, 1, 2, 5, 13};
    //const LO numColsVals[] = {0, 1, 2, 5, 13};
    const LO numRowsVals[] = {1, 2, 5, 13};
    const LO numColsVals[] = {1, 2, 5, 13};
    for (LO numRows : numRowsVals) {
      for (LO numCols : numColsVals) {
        testGemvVsTeuchosBlas<entry_type, coeff_type, device_type> (out,
                                                                    success,
                                                                    randPool,
                                                                    numRows,
                                                                    numCols);
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Blas, Gemv, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_GROUP )

} // namespace (anonymous)


