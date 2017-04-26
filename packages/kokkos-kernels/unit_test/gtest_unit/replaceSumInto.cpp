/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <gtest/gtest.h>
#include "KokkosKernels_Test_Macros.hpp"


// mfh 21 Jun 2016: CUDA 7.5 with GCC 4.8.4 gives me funny build
// errors if I put this functor in an anonymous namespace.  If I name
// the namespace, it builds just fine.
namespace KokkosSparseTest {
  template<class CrsMatrixType>
  class ModifyEvenNumberedRows {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::value_type value_type;

    ModifyEvenNumberedRows (const CrsMatrixType& A,
                            const bool replace,
                            const bool sorted,
                            const bool atomic) :
      A_ (A), replace_ (replace), sorted_ (sorted)
    {}

    KOKKOS_FUNCTION void
    operator () (const ordinal_type& lclRow) const
    {
      if (lclRow % static_cast<ordinal_type> (2) == 0) {
        const ordinal_type ncol = 1;
        ordinal_type cols[1];
        value_type vals[1];

        const value_type ONE = Kokkos::Details::ArithTraits<value_type>::one ();
        const value_type THREE = ONE + ONE + ONE;

        cols[0] = lclRow;
        vals[0] = replace_ ? THREE : ONE;

        if (replace_) {
          A_.replaceValues (lclRow, cols, ncol, vals, sorted_, atomic_);
        }
        else { // sumInto
          A_.sumIntoValues (lclRow, cols, ncol, vals, sorted_, atomic_);
        }
      }
    }

  private:
    CrsMatrixType A_;
    bool replace_;
    bool sorted_;
    bool atomic_;
  };
} // namespace KokkosSparseTest

namespace { // (anonymous)
  using std::endl;

  template<class CrsMatrixType>
  void
  modifyEvenNumberedRows (const CrsMatrixType& A,
                          const bool replace,
                          const bool sorted,
                          const bool atomic)
  {
    typedef typename CrsMatrixType::device_type::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, typename CrsMatrixType::ordinal_type> policy_type;

    ::KokkosSparseTest::ModifyEvenNumberedRows<CrsMatrixType> functor (A, replace, sorted, atomic);
    Kokkos::parallel_for (policy_type (0, A.numRows ()), functor);
  }

  template<class CrsMatrixType>
  bool
  checkWhetherEvenNumberedRowsWereModified (const CrsMatrixType& A,
                                            const bool replace,
                                            const bool /* sorted */,
                                            const bool /* atomic */ )
  {
    typedef typename CrsMatrixType::value_type SC;
    typedef typename CrsMatrixType::ordinal_type LO;

    const SC ONE = Kokkos::Details::ArithTraits<SC>::one ();
    const SC TWO = ONE + ONE;
    const SC THREE = ONE + ONE + ONE;

    typename CrsMatrixType::values_type val = A.values;
    typename CrsMatrixType::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
    Kokkos::deep_copy (val_h, val);

    const LO numRows = A.numRows ();
    bool success = true;
    for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
      if (lclRow % 2 == 0) {
        if (replace && val_h(lclRow) != THREE) {
          success = false;
          break;
        }
        else if (! replace && val_h(lclRow) != TWO) {
          success = false;
          break;
        }
      }
      else if (val_h(lclRow) != ONE) {
        success = false;
        break;
      }
    }
    return success;
  }

  template<class CrsMatrixType>
  void
  testOneCase (bool& success,
               //Teuchos::FancyOStream& out,
          	   std::ostream &out,
               const CrsMatrixType& A,
               const bool replace,
               const bool sorted,
               const bool atomic)
  {
    using Kokkos::Details::ArithTraits;
    typedef typename CrsMatrixType::value_type value_type;

    //Teuchos::OSTab tab0 (out);
    out << "replace: " << (replace ? "true" : "false")
        << ", sorted: " << (sorted ? "true" : "false")
        << ", atomic: " << (atomic ? "true" : "false")
        << endl;

    modifyEvenNumberedRows (A, replace, sorted, atomic);
    const bool lclSuccess =
      checkWhetherEvenNumberedRowsWereModified (A, replace, sorted, atomic);
    KK_TEST_ASSERT( lclSuccess ); // this modifies 'success' and prints to 'out'
    // Restore original values.
    Kokkos::deep_copy (A.values, ArithTraits<value_type>::one ());
  }

  // Test findRelOffset with various array data types and for various cases.
  //
  // This takes the same arguments as if it were declared via the
  // TEUCHOS_UNIT_TEST macro.
  void generalTest (bool& success, std::ostream &out)
		  	  	  	  	  	  	  //Teuchos::FancyOStream& out)
  {
    typedef double SC;
    typedef int LO;
    typedef Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space> DT;
    typedef KokkosSparse::CrsMatrix<SC, LO, DT> matrix_type;

    //Teuchos::OSTab tab0 (out);
    out << "Test KokkosSparse::CrsMatrix::{replace,sumInto}Values*" << endl;
    //Teuchos::OSTab tab1 (out);

    out << "Create a diagonal matrix as a test problem" << endl;

    const LO numRows = 10;
    matrix_type::size_type numEnt = 0; // to be updated below
    matrix_type::row_map_type::non_const_type ptr ("ptr", numRows+1);
    {
      matrix_type::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
      ptr_h[0] = 0;
      for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
        ptr_h[lclRow+1] = ptr_h[lclRow] + 1; // 1 entry in each row
      }
      numEnt = ptr_h[numRows];
      Kokkos::deep_copy (ptr, ptr_h);
    }

    matrix_type::index_type::non_const_type ind ("ind", numEnt);
    {
      matrix_type::index_type::HostMirror ind_h = Kokkos::create_mirror_view (ind);
      for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
        ind_h[lclRow] = lclRow; // diagonal matrix
      }
      Kokkos::deep_copy (ind, ind_h);
    }

    matrix_type::values_type val ("val", numEnt);
    {
      matrix_type::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
      for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
        val_h[lclRow] = 1.0; // diagonal matrix
      }
      Kokkos::deep_copy (val, val_h);
    }

    const LO numCols = numRows; // square diagonal matrix
    matrix_type A ("A", numRows, numCols, numEnt, val, ptr, ind);

    for (int replaceInt = 0; replaceInt < 2; ++replaceInt) {
      const bool replace = replaceInt != 0;
      for (int sortedInt = 0; sortedInt < 2; ++sortedInt) {
        const bool sorted = sortedInt != 0;
        for (int atomicInt = 0; atomicInt < 2; ++atomicInt) {
          const bool atomic = atomicInt != 0;
          testOneCase (success, out, A, replace, sorted, atomic);
        }
      }
    }
  }

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::endl;

  /*
  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  Teuchos::FancyOStream& out = *outPtr;
  */
  std::ostream &out = std::cout;
  out << "Call Kokkos::initialize" << endl;
  Kokkos::initialize (argc, argv);

  bool success = true;
  out << "Run test" << endl;
  generalTest (success, out);

  out << "Call Kokkos::finalize" << endl;
  Kokkos::finalize ();

  if (success) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


