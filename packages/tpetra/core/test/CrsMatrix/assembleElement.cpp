// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "Tpetra_Details_crsMatrixAssembleElement.hpp"
#include <typeinfo> // methods on std::type_info
#include <utility> // std::pair

// CUDA 7.5 at least doesn't like functors in anonymous namespaces.
namespace TpetraTest {

  // Kokkos::parallel_for functor for testing
  // Tpetra::Details::crsMatrixAssembleElement_sortedLinear.
  // MUST be run over an index "range" of a single index, i = 0.
  template<class SparseMatrixType,
           class VectorViewType,
           class RhsViewType,
           class LhsViewType>
  class TestCrsMatrixAssembleElementSortedLinear {
  public:
    typedef typename SparseMatrixType::ordinal_type ordinal_type;
    typedef typename SparseMatrixType::device_type device_type;

    TestCrsMatrixAssembleElementSortedLinear (const SparseMatrixType& A,
                                              const VectorViewType& x,
                                              const Kokkos::View<ordinal_type*, device_type>& lids,
                                              const Kokkos::View<ordinal_type*, device_type>& sortPerm,
                                              const RhsViewType& rhs,
                                              const LhsViewType& lhs,
                                              const bool forceAtomic,
                                              const bool checkInputIndices) :
      A_ (A),
      x_ (x),
      lids_ (lids),
      sortPerm_ (sortPerm),
      rhs_ (rhs),
      lhs_ (lhs),
      forceAtomic_ (forceAtomic),
      checkInputIndices_ (checkInputIndices),
      result_ ("result")
    {
      static_assert (Kokkos::is_view<VectorViewType>::value,
                     "VectorViewType must be a Kokkos::View specialization.");
      static_assert (Kokkos::is_view<RhsViewType>::value,
                     "RhsViewType must be a Kokkos::View specialization.");
      static_assert (Kokkos::is_view<LhsViewType>::value,
                     "LhsViewType must be a Kokkos::View specialization.");
      static_assert (static_cast<int> (RhsViewType::rank) == 1,
                     "RhsViewType must be a rank-1 Kokkos::View.");
      static_assert (static_cast<int> (LhsViewType::rank) == 2,
                     "LhsViewType must be a rank-2 Kokkos::View.");
      static_assert (std::is_integral<ordinal_type>::value,
                     "SparseMatrixType::ordinal_type must be a built-in integer type.");
    }

  // Only meant to be called for i = 0.
  KOKKOS_FUNCTION void
  operator() (const ordinal_type& i) const
  {
    using ::Tpetra::Details::crsMatrixAssembleElement_sortedLinear;

    if (i == 0) {
      const ordinal_type retval =
        crsMatrixAssembleElement_sortedLinear (A_, x_, lids_.data (),
                                               sortPerm_.data (),
                                               rhs_, lhs_, forceAtomic_,
                                               checkInputIndices_);
      result_() = retval;
    }
  }

  ordinal_type
  getReturnValue () const
  {
    auto result_h = Kokkos::create_mirror_view (result_);
    Kokkos::deep_copy (result_h, result_);
    return result_h();
  }

  private:
    SparseMatrixType A_;
    VectorViewType x_;
    Kokkos::View<ordinal_type*, device_type> lids_;
    Kokkos::View<ordinal_type*, device_type> sortPerm_;
    typename RhsViewType::const_type rhs_;
    typename LhsViewType::const_type lhs_;
    bool forceAtomic_;
    bool checkInputIndices_;
    Kokkos::View<ordinal_type, device_type> result_;
  };
} // namespace TpetraTest

namespace { // (anonymous)

  // Call Tpetra::Details::crsMatrixAssembleElement_sortedLinear on
  // device with the given input, and return the function's return
  // value.
  template<class SparseMatrixType,
           class VectorViewType,
           class RhsViewType,
           class LhsViewType>
  //std::pair<typename SparseMatrixType::ordinal_type, std::pair<bool, bool> >
  typename SparseMatrixType::ordinal_type
  testCrsMatrixAssembleElementSortedLinear (const SparseMatrixType& A,
                                            const VectorViewType& x,
                                            const Kokkos::View<typename SparseMatrixType::ordinal_type*, typename SparseMatrixType::device_type>& lids,
                                            const Kokkos::View<typename SparseMatrixType::ordinal_type*, typename SparseMatrixType::device_type>& sortPerm,
                                            const RhsViewType& rhs,
                                            const LhsViewType& lhs,
                                            const typename RhsViewType::const_type& expectedVectorValues,
                                            const typename SparseMatrixType::values_type::const_type& expectedMatrixValues,
                                            const bool forceAtomic =
#ifdef KOKKOS_ENABLE_SERIAL
                                            ! std::is_same<typename SparseMatrixType::device_type::execution_space, Kokkos::Serial>::value,
#else // NOT KOKKOS_ENABLE_SERIAL
                                            false,
#endif // KOKKOS_ENABLE_SERIAL
                                            const bool checkInputIndices = true)
  {
    static_assert (Kokkos::is_view<VectorViewType>::value,
                   "VectorViewType must be a Kokkos::View specialization.");
    static_assert (Kokkos::is_view<RhsViewType>::value,
                   "RhsViewType must be a Kokkos::View specialization.");
    static_assert (Kokkos::is_view<LhsViewType>::value,
                   "LhsViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (RhsViewType::rank) == 1,
                   "RhsViewType must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (LhsViewType::rank) == 2,
                   "LhsViewType must be a rank-2 Kokkos::View.");
    typedef TpetraTest::TestCrsMatrixAssembleElementSortedLinear<SparseMatrixType,
      VectorViewType, RhsViewType, LhsViewType> functor_type;
    typedef typename SparseMatrixType::ordinal_type ordinal_type;
    static_assert (std::is_integral<ordinal_type>::value,
                   "SparseMatrixType::ordinal_type must be a built-in integer type.");
    typedef typename SparseMatrixType::device_type device_type;
    typedef typename device_type::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, ordinal_type> range_type;

    TEUCHOS_TEST_FOR_EXCEPTION
      (expectedVectorValues.extent (0) != x.extent (0),
       std::invalid_argument,
       "expectedVectorValues.extent(0) = " << expectedVectorValues.extent (0)
       << " != x.extent(0) = " << x.extent (0) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (expectedMatrixValues.extent (0) != A.values.extent (0),
       std::invalid_argument,
       "expectedMatrixValues.extent(0) = " << expectedMatrixValues.extent (0)
       << " != A.values.extent(0) = " << A.values.extent (0) << ".");

    functor_type functor (A, x, lids, sortPerm, rhs, lhs, forceAtomic, checkInputIndices);
    // It's a "parallel" loop with one loop iteration.  The point is
    // to run on device.
    Kokkos::parallel_for (range_type (0, 1), functor);
    Kokkos::fence (); // make sure that the kernel finished

    return functor.getReturnValue ();
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsMatrix, assembleElement, ScalarType )
  {
    using std::endl;
    typedef typename Kokkos::ArithTraits<ScalarType>::val_type SC;
    typedef Tpetra::Map<>::local_ordinal_type LO;
    //typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef Tpetra::Map<>::device_type DT;
    typedef Tpetra::Map<> map_type;
    typedef KokkosSparse::CrsMatrix<SC, LO, DT, void> sparse_matrix_type;
    typedef typename sparse_matrix_type::size_type offset_type;
    typedef typename sparse_matrix_type::StaticCrsGraphType sparse_graph_type;

    out << "Test crsMatrixAssembleElement" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Create Map just to initialize Kokkos correctly" << endl;
    auto comm = Tpetra::getDefaultComm ();
    map_type mapToInitKokkos (Tpetra::global_size_t (100), 0, comm);

    const bool execSpaceInitd = Kokkos::is_initialized ();
    TEST_ASSERT( execSpaceInitd );
    if (! execSpaceInitd) {
      out << "Tpetra::Map failed to initialize Kokkos execution space \""
          << typeid (DT::execution_space).name () << "\"." << endl;
      return;
    }

    // Put the entire test (not counting the above initialization) in
    // its own scope.  That way, Tpetra::Map's destructor (actually,
    // the Node's destructor) will not finalize Kokkos until all the
    // Kokkos::View objects have fallen out of scope.
    {
      // Dimension of the elements to test.
      constexpr LO eltDim = 4;
      // Sparsity pattern of elements to use for scattering into the matrix.
      const LO eltSparsityPattern[eltDim] = {3, 5, 8, 11};

      // Values for the element to scatter into the matrix.  Choose
      // unique values, to make sure we get the indices right.
      const SC ONE = Kokkos::ArithTraits<SC>::one ();

      out << "Constructing element matrix" << endl;

      Kokkos::View<SC**, DT> lhs_d ("lhs_d", eltDim, eltDim);
      {
        auto lhs_h = Kokkos::create_mirror_view (lhs_d);

        SC curVal = ONE;
        for (LO i = 0; i < eltDim; ++i) {
          for (LO j = 0; j < eltDim; ++j) {
            lhs_h(i,j) = curVal;
            curVal = curVal + ONE;
          }
        }
        out << "Element matrix (lhs): [";
        for (LO i = 0; i < static_cast<LO> (lhs_h.extent (0)); ++i) {
          for (LO j = 0; j < static_cast<LO> (lhs_h.extent (1)); ++j) {
            constexpr int width = Kokkos::ArithTraits<SC>::is_complex ? 7 : 3;
            out << std::setw (width) << lhs_h(i,j);
            if (j + LO (1) < static_cast<LO> (lhs_h.extent (1))) {
              out << " ";
            }
          }
          if (i + LO (1) < static_cast<LO> (lhs_h.extent (0))) {
            out << endl;
          }
        }
        out << "]" << endl;

        Kokkos::deep_copy (lhs_d, lhs_h);
      }

      out << "Constructing element vector" << endl;
      Kokkos::View<SC*, DT> rhs_d ("rhs_d", eltDim);
      {
        auto rhs_h = Kokkos::create_mirror_view (rhs_d);
        SC curVal = ONE;
        for (LO i = 0; i < static_cast<LO> (rhs_h.extent (0)); ++i) {
          rhs_h(i) = -curVal;
          curVal = curVal + ONE;
        }
        Kokkos::deep_copy (rhs_d, rhs_h);
      }

      // Number of rows in the sparse matrix A, and number of entries in
      // the dense vector b.  This is a constexpr so we can easily
      // construct arrays with it.
      constexpr LO numRows = 13;
      // Number of columns in the sparse matrix A
      const LO numCols = numRows;

      // When defining the matrix sparsity pattern, make sure that some
      // rows not in the above element sparsity pattern just happen to
      // have some column indices that match those in the element
      // sparsity pattern.
      const LO matSparsityPattern[] = {
        0, 4, 5, 7, 10,            // row 0 (matches none of the indices)
        // row 1 is empty
        2, 3, 7, 9, 11,            // row 2 (matches some of the indices)
        0, 3, 5, 7, 8, 10, 11, 12, // row 3 (all indices, plus extra)
        5, 6, 8,                   // row 4 (matches some of the indices)
        1, 3, 5, 6, 8, 11,         // row 5 (all indices, plus extra)
        2, 10, 12,                 // row 6 (matches none of the indices)
        3, 5, 8, 11,               // row 7 (happens to match all the indices)
        3, 5, 8, 11,               // row 8 (exactly the indices, no more)
        // row 9 is empty
        // row 10 is empty (two consecutive empty rows)
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, // row 11 (all indices, plus extra)
        5, 10                      // row 12 (matches some of the indices)
      };

      // Number of entries in each row of the sparse matrix.
      const offset_type numEntPerRow[] = {
        5,
        0,
        5,
        8,
        3,
        6,
        3,
        4,
        4,
        0,
        0,
        13,
        2
      };

#if 0
      // Pattern of entries in the matrix that an assembleElement
      // operation should have changed.
      const LO matChangedPattern[] = {
        0, 0, 0, 0, 0,             // row 0 (matches none of the indices)
        // row 1 is empty
        0, 0, 0, 0, 0,             // row 2 (matches some of the indices)
        0, 1, 1, 0, 1, 0, 1, 0,    // row 3 (all indices, plus extra)
        0, 0, 0,                   // row 4 (matches some of the indices)
        0, 1, 1, 0, 1, 1,         // row 5 (all indices, plus extra)
        0, 0, 0,                   // row 6 (matches none of the indices)
        0, 0, 0, 0,                // row 7 (happens to match all the indices)
        1, 1, 1, 1,               // row 8 (exactly the indices, no more)
        // row 9 is empty
        // row 10 is empty (two consecutive empty rows)
        0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, // row 11 (all indices, plus extra)
        0, 0                      // row 12 (matches some of the indices)
      };
#endif // 0

      // Total number of entries in the sparse matrix.
      const offset_type numEnt = 53;

      // Expected values in the matrix after the assembleElement operation.
      out << "Make the array of expected matrix entries" << endl;
      typedef typename sparse_matrix_type::values_type::non_const_type values_type;
      values_type expectedMatrixValues ("expectedMatrixValues", numEnt);
      {
        Teuchos::OSTab tabX (out);
        out << "Allocated device view" << endl;
        auto expectedMatrixValues_h = Kokkos::create_mirror_view (expectedMatrixValues);
        out << "Created mirror view" << endl;
        Kokkos::deep_copy (expectedMatrixValues_h, Kokkos::ArithTraits<SC>::zero ());
        out << "Filled mirror view with 0s" << endl;

        SC curVal = ONE;
        expectedMatrixValues_h(11) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(12) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(14) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(16) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(22) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(23) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(25) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(26) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(34) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(35) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(36) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(37) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(41) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(43) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(46) = curVal;
        curVal += ONE;
        expectedMatrixValues_h(49) = curVal;
        //curVal += ONE;
        Kokkos::deep_copy (expectedMatrixValues, expectedMatrixValues_h);
        out << "Copied to device" << endl;
      }

      // Expected values in the right-hand-side vector after the
      // assembleElement operation.
      out << "Make the array of expected right-hand-side vector entries" << endl;
      Kokkos::View<SC*, DT> expectedVectorValues ("expectedVectorValues", numRows);
      {
        auto expectedVectorValues_h = Kokkos::create_mirror_view (expectedVectorValues);
        Kokkos::deep_copy (expectedVectorValues_h, Kokkos::ArithTraits<SC>::zero ());
        SC curVal = ONE;
        expectedVectorValues_h(3) = -curVal;
        curVal += ONE;
        expectedVectorValues_h(5) = -curVal;
        curVal += ONE;
        expectedVectorValues_h(8) = -curVal;
        curVal += ONE;
        expectedVectorValues_h(11) = -curVal;
        //curVal += ONE;
        Kokkos::deep_copy (expectedVectorValues, expectedVectorValues_h);
      }

      out << "Create the sparse matrix A, and fill its values with zeros" << endl;
      typename sparse_graph_type::row_map_type::non_const_type
        A_ptr ("A_ptr", numRows+1);
      {
        auto A_ptr_host = Kokkos::create_mirror_view (A_ptr);
        A_ptr_host[0] = 0;
        for (offset_type i = 1; i <= static_cast<offset_type> (numRows); ++i) {
          A_ptr_host[i] = A_ptr_host[i-1] + numEntPerRow[i-1];
        }
        TEST_EQUALITY( A_ptr_host[numRows], numEnt );
        if (A_ptr_host[numRows] != numEnt) {
          out << "The sparse matrix for the test was not constructed "
              << "correctly, since the last entry of the offsets array, "
              << "A_ptr_host[numRows=" << numRows << "] != numEnt = " << numEnt
              << ".  Please go back and fix it." << endl;
          return; // no sense in continuing the test at this point
        }
        Kokkos::deep_copy (A_ptr, A_ptr_host);
      }
      typename sparse_graph_type::entries_type::non_const_type
        A_ind ("A_ind", numEnt);
      {
        auto A_ind_host = Kokkos::create_mirror_view (A_ind);
        for (size_t k = 0; k < static_cast<size_t> (numEnt); ++k) {
          A_ind_host[k] = matSparsityPattern[k];
        }
        Kokkos::deep_copy (A_ind, A_ind_host);
      }
      sparse_graph_type A_graph (A_ind, A_ptr);
      sparse_matrix_type A ("A", A_graph, numCols);
      Kokkos::deep_copy (A.values, Kokkos::ArithTraits<SC>::zero ());

      out << "The sparse matrix, as copied to device:" << endl;
      {
        Teuchos::OSTab tabX (out);

        auto ptr = Kokkos::create_mirror_view (A.graph.row_map);
        Kokkos::deep_copy (ptr, A.graph.row_map);
        out << "ptr: [";
        for (LO k = 0; k < static_cast<LO> (ptr.extent (0)); ++k) {
          out << ptr(k);
          if (k + LO (1) < static_cast<LO> (ptr.extent (0))) {
            out << ",";
          }
        }
        out << "]" << endl;
        auto ind = Kokkos::create_mirror_view (A.graph.entries);
        Kokkos::deep_copy (ind, A.graph.entries);
        out << "ind: [";
        for (offset_type k = 0; k < static_cast<offset_type> (ind.extent (0)); ++k) {
          out << ind(k);
          if (k + offset_type (1) < static_cast<offset_type> (ind.extent (0))) {
            out << ",";
          }
        }
        out << "]" << endl;
        auto val = Kokkos::create_mirror_view (A.values);
        Kokkos::deep_copy (val, A.values);
        out << "val: [";
        for (offset_type k = 0; k < static_cast<offset_type> (val.extent (0)); ++k) {
          out << val(k);
          if (k + offset_type (1) < static_cast<offset_type> (val.extent (0))) {
            out << ",";
          }
        }
        out << "]" << endl;
      }

      out << "Create the \"right-hand side\" vector b, and fill it with zeros" << endl;
      Kokkos::View<SC*, DT> b ("b", numRows);
      Kokkos::deep_copy (b, Kokkos::ArithTraits<SC>::zero ());

      out << "Make the list of indices input/output array" << endl;
      Kokkos::View<LO*, DT> lids_d ("lids", eltDim);
      {
        Kokkos::View<const LO*, typename Kokkos::View<LO*, DT>::array_layout,
          Kokkos::HostSpace, Kokkos::MemoryUnmanaged> lids_h (eltSparsityPattern, eltDim);
        //typename Kokkos::View<const LO*, DT>::HostMirror lids_h (eltSparsityPattern, eltDim);
        Kokkos::deep_copy (lids_d, lids_h);
      }

      out << "Make the sort permutation output array" << endl;
      Kokkos::View<LO*, DT> sortPerm_d ("sortPerm", eltDim);

      out << "Call the function to test" << endl;
      auto retval =
        testCrsMatrixAssembleElementSortedLinear (A, b, lids_d, sortPerm_d,
                                                  rhs_d, lhs_d,
                                                  expectedVectorValues,
                                                  expectedMatrixValues);
      const LO numEntFound = retval; // retval.first;
      TEST_EQUALITY( numEntFound, eltDim*eltDim );
      // TEST_ASSERT( retval.second.first );
      // TEST_ASSERT( retval.second.second );
      out << "Function returned numEntFound=" << numEntFound << endl;

      {
        auto A_val_h = Kokkos::create_mirror_view (A.values);
        Kokkos::deep_copy (A_val_h, A.values);
        auto val_h = Kokkos::create_mirror_view (expectedMatrixValues);
        Kokkos::deep_copy (val_h, expectedMatrixValues);
        TEST_EQUALITY( A_val_h.extent (0), val_h.extent (0) );
        if (A_val_h.extent (0) == val_h.extent (0)) {
          bool same = true;
          const offset_type len = A_val_h.extent (0);
          for (offset_type k = 0; k < len; ++k) {
            if (A_val_h(k) != val_h(k)) {
              same = false;
              break;
            }
          }
          TEST_ASSERT( same );
        }
        out << "A.values            : [";
        for (offset_type k = 0; k < numEnt; ++k) {
          constexpr int width = Kokkos::ArithTraits<SC>::is_complex ? 7 : 3;
          out << std::setw (width) << A_val_h(k);
          if (k + offset_type (1) < numEnt) {
            out << ",";
          }
        }
        out << "]" << endl;
        out << "expectedMatrixValues: [";
        for (offset_type k = 0; k < numEnt; ++k) {
          constexpr int width = Kokkos::ArithTraits<SC>::is_complex ? 7 : 3;
          out << std::setw (width) << val_h(k);
          if (k + offset_type (1) < numEnt) {
            out << ",";
          }
        }
        out << "]" << endl;
      }
      {
        auto b_h = Kokkos::create_mirror_view (b);
        Kokkos::deep_copy (b_h, b);
        auto b_exp_h = Kokkos::create_mirror_view (expectedVectorValues);
        Kokkos::deep_copy (b_exp_h, expectedVectorValues);
        TEST_EQUALITY( b_h.extent (0), b_exp_h.extent (0) );
        if (b_h.extent (0), b_exp_h.extent (0)) {
          bool same = true;
          const offset_type len = b_h.extent (0);
          for (offset_type k = 0; k < len; ++k) {
            if (b_h(k) != b_exp_h(k)) {
              same = false;
              break;
            }
          }
          TEST_ASSERT( same );
        }
        out << "Actual output b  : [";
        for (offset_type k = 0;
             k < static_cast<offset_type> (b_h.extent (0)); ++k) {
          out << b_h(k);
          if (k + offset_type (1) <
              static_cast<offset_type> (b_h.extent (0))) {
            out << ",";
          }
        }
        out << "]" << endl;
        out << "Expected output b: [";
        for (offset_type k = 0;
             k < static_cast<offset_type> (b_exp_h.extent (0)); ++k) {
          out << b_exp_h(k);
          if (k + offset_type (1) <
              static_cast<offset_type> (b_exp_h.extent (0))) {
            out << ",";
          }
        }
        out << "]" << endl;
      }
    }
    // The above curly brace invokes the destructors of all remaining
    // Kokkos::View objects.  At this point, all Kokkos::View objects'
    // reference counts should have gone to zero, and all of them
    // should have been deallocated.  Thus, it is safe to call
    // Kokkos::finalize now.  Tpetra::Map's destructor (actually, the
    // Node's destructor) should have handled that for us.
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, assembleElement, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

// FIXME_SYCL
#ifndef KOKKOS_ENABLE_SYCL
  UNIT_TEST_GROUP( double )
#endif
  //TPETRA_INSTANTIATE_S_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)

