/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#ifndef KOKKOS_SPARSE_IMPL_MKL_HPP_
#define KOKKOS_SPARSE_IMPL_MKL_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include "Kokkos_Sparse_impl_copyIntegers.hpp"

#ifdef HAVE_TPETRAKERNELS_MKL
#  include "mkl_spblas.h"
#endif // HAVE_TPETRAKERNELS_MKL


namespace KokkosSparse {
namespace Impl {
namespace Mkl {

#ifndef HAVE_TPETRAKERNELS_MKL
typedef int MKL_INT;
typedef int sparse_diag_type_t;
typedef int sparse_fill_mode_t;
typedef int sparse_index_base_t;
typedef int sparse_layout_t;
typedef int sparse_matrix_t;
typedef int sparse_matrix_type_t;
typedef int sparse_memory_usage_t;
typedef int sparse_operation_t;
typedef int sparse_status_t;

// Fall-back typedefs, so this file builds even if MKL is not enabled.
struct matrix_descr {
  sparse_matrix_type_t type;
  sparse_fill_mode_t mode;
  sparse_diag_type_t diag;
};
#endif // NOT HAVE_TPETRAKERNELS_MKL

/// \brief Does MKL support operations with sparse matrices whose
///   entries have type ValueType?
/// \tparam ValueType Type of the entries in the sparse matrix.
///
/// \c value is true if and only if MKL supports sparse operations
/// with matrices whose value_type is \c ValueType.
template<class ValueType>
struct SupportsValueType {
  /// \brief Whether MKL supports operations with sparse matrices
  ///   whose entries have type \c ValueType.
  static const bool value =
#ifdef HAVE_TPETRAKERNELS_MKL
    std::is_same<ValueType, double>::value ||
    std::is_same<ValueType, float>::value ||
    std::is_same<ValueType, ::Kokkos::complex<double> >::value ||
    std::is_same<ValueType, ::Kokkos::complex<float> >::value;
#else
    false;
#endif // HAVE_TPETRAKERNELS_MKL
};

//! Value returned by an MKL operation indicating success.
constexpr sparse_status_t tplStatusSuccessful () {
#ifdef HAVE_TPETRAKERNELS_MKL
  return SPARSE_STATUS_SUCCESS;
#else
  return 0; // the usual *nix convention for success
#endif // HAVE_TPETRAKERNELS_MKL
}

/// \brief Value returned by an MKL operation indicating that the
///   operation is not supported.
constexpr sparse_status_t tplStatusNotSupported () {
#ifdef HAVE_TPETRAKERNELS_MKL
  return SPARSE_STATUS_NOT_SUPPORTED;
#else
  return -1; // must not be the same as "success" (see above)
#endif // HAVE_TPETRAKERNELS_MKL
}

/// \brief Value returned by an MKL operation indicating that the
///   matrix is not initialized.
constexpr sparse_status_t tplStatusNotInitialized () {
#ifdef HAVE_TPETRAKERNELS_MKL
  return SPARSE_STATUS_NOT_INITIALIZED;
#else
  return -2; // must not be the same as "success" (see above)
#endif // HAVE_TPETRAKERNELS_MKL
}

/// \brief Convert "status" return value of an MKL operation into a
///   human-readable string.
std::string tplStatusToString (const sparse_status_t status) {
#ifdef HAVE_TPETRAKERNELS_MKL
  // All strings but the last (default case) come directly from
  // Intel's MKL documentation:
  // https://software.intel.com/en-us/node/590115
  // [last accessed 26 Aug 2016].
  if (status == SPARSE_STATUS_SUCCESS) {
    return "The operation was successful.";
  }
  else if (status == SPARSE_STATUS_NOT_INITIALIZED) {
    return "The routine encountered an empty handle or matrix array.";
  }
  else if (status == SPARSE_STATUS_ALLOC_FAILED) {
    return "Internal memory allocation failed.";
  }
  else if (status == SPARSE_STATUS_INVALID_VALUE) {
    return "The input parameters contain an invalid value.";
  }
  else if (status == SPARSE_STATUS_EXECUTION_FAILED) {
    return "Execution failed.";
  }
  else if (status == SPARSE_STATUS_INTERNAL_ERROR) {
    return "An error in algorithm implementation occurred.";
  }
  else if (status == SPARSE_STATUS_NOT_SUPPORTED) {
    return "The requested operation is not supported.";
  }
  else {
    return "Invalid status value.";
  }
#else
  return "Invalid status value.";
#endif // HAVE_TPETRAKERNELS_MKL
}

/// \brief Operations on a "raw" MKL sparse matrix handle.
/// \tparam ValueType Type of the entries in the sparse matrix.
template<class ValueType>
struct RawTplMatrixHandle {
  //! Type of the entries of the sparse matrix, as Trilinos stores them.
  typedef ValueType value_type;
  /// \brief Type of the entries of the sparse matrix, as MKL stores them.
  ///
  /// This may differ from value_type, especially if value_type is
  /// complex.  However, the two types must always be bitwise
  /// equivalent.
  typedef ValueType internal_value_type;
  typedef sparse_matrix_t handle_type;
  typedef sparse_status_t status_type;

  static status_type
  create (handle_type* /* handle */, // output argument
          sparse_index_base_t /* indexing */,
          MKL_INT /* numRows */,
          MKL_INT /* numCols */,
          MKL_INT* /* rowBeg */,
          MKL_INT* /* rowEnd */,
          MKL_INT* /* colInd */,
          value_type* /* values */)
  {
    // This ValueType is not supported.  Exploit existing return value.
    return tplStatusNotSupported ();
  }

  static status_type destroy (handle_type /* handle */) {
    // This ValueType is not supported.  Exploit existing return value.
    return tplStatusNotSupported ();
  }
};

//! Full specialization of RawTplMatrixHandle for ValueType = double.
template<>
struct RawTplMatrixHandle<double> {
  typedef double value_type;
  typedef double internal_value_type;
  typedef sparse_matrix_t handle_type;
  typedef sparse_status_t status_type;

  static status_type
  create (handle_type* handle, // output argument
          sparse_index_base_t indexing,
          MKL_INT numRows,
          MKL_INT numCols,
          MKL_INT* rowBeg,
          MKL_INT* rowEnd,
          MKL_INT* colInd,
          value_type* values)
  {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_d_create_csr (handle, indexing, numRows, numCols,
                                    rowBeg, rowEnd, colInd, values);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }

  static status_type destroy (handle_type handle) {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_destroy (handle);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }
};

//! Full specialization of RawTplMatrixHandle for ValueType = float.
template<>
struct RawTplMatrixHandle<float> {
  typedef float value_type;
  typedef float internal_value_type;
  typedef sparse_matrix_t handle_type;
  typedef sparse_status_t status_type;

  static status_type
  create (handle_type* handle, // output argument
          sparse_index_base_t indexing,
          MKL_INT numRows,
          MKL_INT numCols,
          MKL_INT* rowBeg,
          MKL_INT* rowEnd,
          MKL_INT* colInd,
          value_type* values)
  {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_s_create_csr (handle, indexing, numRows, numCols,
                                    rowBeg, rowEnd, colInd, values);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }

  static status_type destroy (handle_type handle) {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_destroy (handle);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }
};

//! Full specialization of RawTplMatrixHandle for complex (double).
template<>
struct RawTplMatrixHandle< ::Kokkos::complex<double> > {
  typedef ::Kokkos::complex<double> value_type;
#ifdef HAVE_TPETRAKERNELS_MKL
  typedef MKL_Complex16 internal_value_type;
#else
  typedef value_type internal_value_type;
#endif // HAVE_TPETRAKERNELS_MKL
  typedef sparse_matrix_t handle_type;
  typedef sparse_status_t status_type;

  static status_type
  create (handle_type* handle, // output argument
          sparse_index_base_t indexing,
          MKL_INT numRows,
          MKL_INT numCols,
          MKL_INT* rowBeg,
          MKL_INT* rowEnd,
          MKL_INT* colInd,
          value_type* values)
  {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_z_create_csr (handle, indexing, numRows, numCols,
                                    rowBeg, rowEnd, colInd,
                                    reinterpret_cast<internal_value_type*> (values));
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }

  static status_type destroy (handle_type handle) {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_destroy (handle);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }
};

//! Full specialization of RawTplMatrixHandle for complex (float).
template<>
struct RawTplMatrixHandle< ::Kokkos::complex<float> > {
  typedef ::Kokkos::complex<float> value_type;
#ifdef HAVE_TPETRAKERNELS_MKL
  typedef MKL_Complex8 internal_value_type;
#else
  typedef value_type internal_value_type;
#endif // HAVE_TPETRAKERNELS_MKL
  typedef sparse_matrix_t handle_type;
  typedef sparse_status_t status_type;

  static status_type
  create (handle_type* handle, // output argument
          sparse_index_base_t indexing,
          MKL_INT numRows,
          MKL_INT numCols,
          MKL_INT* rowBeg,
          MKL_INT* rowEnd,
          MKL_INT* colInd,
          value_type* values)
  {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_c_create_csr (handle, indexing, numRows, numCols,
                                    rowBeg, rowEnd, colInd,
                                    reinterpret_cast<internal_value_type*> (values));
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }

  static status_type destroy (handle_type handle) {
#ifdef HAVE_TPETRAKERNELS_MKL
    return mkl_sparse_destroy (handle);
#else
    return tplStatusNotSupported ();
#endif // HAVE_TPETRAKERNELS_MKL
  }
};


// "Deleter" for std::shared_ptr of the "raw" TPL handle.  This thing
// "frees" the handle once its reference count goes to zero.  "Frees"
// means telling MKL to destroy the (raw) handle.
template<class MatrixType>
struct RawTplMatrixHandleDeleter {
  void
  operator() (typename RawTplMatrixHandle<MatrixType>::handle_type* handle) const
  {
    if (handle != NULL) {
      // This is a destructor, so we shouldn't throw if not successful.
      (void) RawTplMatrixHandle<MatrixType>::destroy (*handle);
      delete handle;
    }
  }
};

/// \brief Wrapped MKL handle for a sparse matrix
/// \note Tpetra, Ifpack2, etc. developers: USE THIS CLASS!
/// \tparam MatrixType KokkosSparse::CrsMatrix specialization
///
/// Tpetra, Ifpack2, etc. developers may use this class to create and
/// access the MKL handle for a KokkosSparse::CrsMatrix.  If they need
/// to share this handle, they should wrap it in std::shared_ptr or
/// Teuchos::RCP.  Direct copying doesn't make sense.
///
/// You must call setMatrix under any of the following conditions:
/// <ul>
/// <li> If any of the matrix's pointers change </li>
/// <li> If the graph structure of the matrix changes </li>
/// <li> If any values in the matrix change </li>
/// </ul>
/// (MKL technically lets users change the values without recreating
/// its handle, but we assume that this is not allowed either.  This
/// is because it is likely cheaper for users to change the matrix
/// through Trilinos' data structure.)
///
/// It is better to keep the WrappedTplMatrixHandle instance around
/// than to create a new one, because the instance can reuse internal
/// index storage.  (It needs internal storage because MKL may store
/// row offsets using a different integer type than Trilinos uses.)
template<class MatrixType>
class WrappedTplMatrixHandle {
private:
  typedef typename std::remove_const<typename MatrixType::value_type>::type value_type;
  typedef typename RawTplMatrixHandle<MatrixType>::handle_type handle_type;
  typedef typename RawTplMatrixHandle<MatrixType>::status_type status_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

public:
  WrappedTplMatrixHandle () :
    status_ (tplStatusNotInitialized ())
  {}

  WrappedTplMatrixHandle (const MatrixType& A,
                          const bool reuseGraph) :
    status_ (tplStatusNotInitialized ())
  {
    setMatrix (A, reuseGraph);
  }

  //! Free any state that needs freeing.
  void
  reset ()
  {
    if (status_ == tplStatusSuccessful ()) {
      // This is a destructor of sorts, so don't attempt to detect
      // errors on failed destruction of a member.
      (void) RawTplMatrixHandle<value_type>::destroy (handle_);
    }
    status_ = tplStatusNotInitialized ();
    ptrBeg_ = decltype (ptrBeg_) ();
    ptrEnd_ = decltype (ptrEnd_) ();
    colInd_ = decltype (colInd_) ();
  }

  /// \brief Set up the MKL handle for the given sparse matrix A.
  ///
  /// \param A [in] The sparse matrix for which to set up MKL
  /// \param reuseGraph [in] Whether to reuse any previously allocated
  ///   graph information
  ///
  /// This method does not satisfy the strong exception guarantee.
  /// This is an optimization to improve reuse of previously allocated
  /// memory.  However, if this method throws, it will first free or
  /// invalidate any of its state.
  void
  setMatrix (const MatrixType& A, const bool reuseGraph)
  {
#ifdef HAVE_TPETRAKERNELS_MKL
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
#else
    sparse_index_base_t indexing = 0;
#endif // HAVE_TPETRAKERNELS_MKL

    // Get the matrix's dimensions.  Make sure that we can cast them
    // to MKL_INT without overflow.
    MKL_INT numRows = 0;
    MKL_INT numCols = 0;
    if (sizeof (MKL_INT) > sizeof (ordinal_type) ||
        (sizeof (MKL_INT) == sizeof (ordinal_type) &&
         std::numeric_limits<MKL_INT>::is_signed &&
         ! std::numeric_limits<ordinal_type>::is_signed)) {
      // Assigning ordinal_type to MKL_INT could overflow.
      // In that case, we need to check.

      const ordinal_type numRowsOT = A.numRows ();
      const ordinal_type numColsOT = A.numCols ();
      if (numRowsOT > static_cast<ordinal_type> (std::numeric_limits<MKL_INT>::max ()) ||
          numColsOT > static_cast<ordinal_type> (std::numeric_limits<MKL_INT>::max ())) {
        std::ostringstream os;
        os << "A's dimensions do not fit in MKL_INT.  A is " << numRowsOT
           << " x " << numColsOT << ", but the maximum allowed MKL_INT value "
          "is " << std::numeric_limits<MKL_INT>::max () << ".";
        throw std::runtime_error (os.str ());
      }
      numRows = static_cast<MKL_INT> (numRowsOT);
      numCols = static_cast<MKL_INT> (numColsOT);
    }
    else { // assigning ordinal_type to MKL_INT won't overflow
      numRows = A.numRows ();
      numCols = A.numCols ();
    }

    // It may be that the type of the matrix's row offsets is not the
    // same as MKL_INT.  In that case, attempt to copy from the given
    // type into MKL_INT.  Throw on overflow.
    //
    // This code will attempt to reuse output storage if it has
    // already been allocated.  However, it will always copy the data,
    // unless you tell it to "reuse the graph."  This means that it
    // won't copy if the output allocation already exists.
    {
      constexpr bool sameOffsetType =
        std::is_same<MKL_INT, typename std::decay<decltype (A.graph.row_map[0]) >::type>::value;
      decltype (ptrBeg_) ptrBegEmpty;
      // If the offset types are the same, attempt to assign the input
      // array of row offsets to ptrBeg_.  Otherwise, "assign" an
      // empty array.
      decltype (ptrBeg_) tmpPtrBeg =
        Kokkos::Impl::if_c<sameOffsetType,
        decltype (A.graph.row_map),
        decltype (ptrBeg_) >::select (A.graph.row_map, ptrBegEmpty);
      if (tmpPtrBeg.dimension_0 () == A.graph.row_map.dimension_0 ()) {
        ptrBeg_nc_ = decltype (ptrBeg_nc_) (); // clear out existing storage
        ptrBeg_ = tmpPtrBeg; // done
      }
      else {
        if (ptrBeg_nc_.dimension_0 () != A.graph.row_map.dimension_0 ()) {
          ptrBeg_nc_ = decltype (ptrBeg_nc_) ("ptr", A.graph.row_map.dimension_0 ());
          copyIntegers (ptrBeg_nc_, A.graph.row_map);
          ptrBeg_ = ptrBeg_nc_;
        }
        else if (! reuseGraph) {
          copyIntegers (ptrBeg_nc_, A.graph.row_map);
          ptrBeg_ = ptrBeg_nc_;
        }
      }
    }

    if (A.graph.row_map.dimension_0 () != 0) {
      Kokkos::pair<MKL_INT, MKL_INT> range (1, A.graph.row_map.dimension_0 ());
      // NOTE (mfh 26 Aug 2016) ptrBeg should really have length
      // numRows, not numRows+1.  However, MKL doesn't care, and it's
      // easier (given that we use 3-array CSR in Trilinos) to make
      // ptrBeg have the same length as our row offsets array, namely
      // numRows+1.  Thus, we can make ptrEnd an offset view of
      // ptrBeg.
      ptrEnd_ = Kokkos::subview (ptrBeg_, range);
    }
    else {
      ptrEnd_ = decltype (ptrEnd_) (); // empty View
    }
    // MKL offers a C interface, and therefore does not accept const
    // qualifiers.  However, it won't change our data, so it's safe to
    // cast away const.  We do so here and below.
    MKL_INT* rowBegRaw = const_cast<MKL_INT*> (ptrBeg_.ptr_on_device ());
    MKL_INT* rowEndRaw = const_cast<MKL_INT*> (ptrEnd_.ptr_on_device ());

    {
      constexpr bool sameOrdinalType =
        std::is_same<MKL_INT, typename std::decay<decltype (A.graph.entries[0]) >::type>::value;
      decltype (colInd_) colIndEmpty;
      // If the column index ("ordinal") types are the same, attempt
      // to assign the input array of column indices to colInd_.
      // Otherwise, "assign" an empty array.
      decltype (colInd_) tmpColInd =
        Kokkos::Impl::if_c<sameOrdinalType,
        decltype (A.graph.entries),
        decltype (colInd_) >::select (A.graph.entries, colIndEmpty);
      if (tmpColInd.dimension_0 () == A.graph.entries.dimension_0 ()) {
        colInd_nc_ = decltype (colInd_nc_) (); // clear out existing storage
        colInd_ = tmpColInd; // done
      }
      else {
        if (colInd_nc_.dimension_0 () != A.graph.entries.dimension_0 ()) {
          colInd_nc_ = decltype (colInd_nc_) ("ind", A.graph.entries.dimension_0 ());
          copyIntegers (colInd_nc_, A.graph.entries);
          colInd_ = colInd_nc_;
        }
        else if (! reuseGraph) {
          copyIntegers (colInd_nc_, A.graph.entries);
          colInd_ = colInd_nc_;
        }
      }
    }

    MKL_INT* colIndRaw = const_cast<MKL_INT*> (colInd_.ptr_on_device ());
    value_type* valRaw = const_cast<value_type*> (A.values.ptr_on_device ());

    // This unique_ptr won't attempt to free the handle through MKL's
    // interface.  It will just call 'delete' on the allocated struct.
    // That's OK; that's actually what we want.  If this operation
    // succeeds, we'll set the custom deallocator below.  We use
    // unique_ptr because it's cheaper than shared_ptr, and we don't
    // need shared_ptr's generality here.

    handle_type newHandle;
    const status_type newStatus =
      RawTplMatrixHandle<value_type>::create (&newHandle, indexing,
                                              numRows, numCols, rowBegRaw,
                                              rowEndRaw, colIndRaw, valRaw);

    if (newStatus == tplStatusSuccessful ()) {
      handle_ = newHandle; // shallow copy of the handle
      status_ = newStatus;
    }
    else {
      reset (); // on failure, don't leak memory
      const std::string statusStr = tplStatusToString (newStatus);
      std::ostringstream os;
      os << "MKL failed to create the sparse matrix handle.  "
        "It reports the following error message: " << statusStr;
      throw std::runtime_error (os.str ());
    }
  }

  /// \note Get MKL's sparse matrix handle.
  ///
  /// \warning Do NOT call mkl_sparse_destroy on the return value!!!
  handle_type getHandle () const {
    if (status_ == tplStatusSuccessful ()) {
      return handle_;
    }
    else {
      throw std::runtime_error ("Sparse matrix handle is not ready!");
    }
  }

private:
  /// \brief "Raw" MKL handle to the sparse matrix.
  ///
  /// This handle is valid if and only if
  /// status_ == tplStatusSuccessful().
  handle_type handle_;
  //! Status returned from creating TPL handle.
  status_type status_;

  //! Offsets of row starts; may be a deep copy if MKL_INT != ordinal_type.
  ::Kokkos::View<const MKL_INT*,
                 typename MatrixType::row_map_type::array_layout,
                 typename MatrixType::row_map_type::device_type> ptrBeg_;
  //! Nonconst version of ptrBeg_; only used if row offsets not MKL_INT.
  ::Kokkos::View<MKL_INT*,
                 typename MatrixType::row_map_type::array_layout,
                 typename MatrixType::row_map_type::device_type> ptrBeg_nc_;

  //! Offsets of row ends; may be a deep copy if MKL_INT != ordinal_type.
  ::Kokkos::View<const MKL_INT*,
                 typename MatrixType::row_map_type::array_layout,
                 typename MatrixType::row_map_type::device_type> ptrEnd_;

  //! Column indices; may be a deep copy if MKL_INT != ordinal_type.
  ::Kokkos::View<const MKL_INT*,
                 typename MatrixType::StaticCrsGraphType::entries_type::array_layout,
                 typename MatrixType::StaticCrsGraphType::entries_type::device_type> colInd_;
  //! Nonconst version of colInd_; only used if column indices not MKL_INT.
  ::Kokkos::View<MKL_INT*,
                 typename MatrixType::StaticCrsGraphType::entries_type::array_layout,
                 typename MatrixType::StaticCrsGraphType::entries_type::device_type> colInd_nc_;

  /// \brief Forbid copy construction syntactically.
  ///
  /// Trilinos' sparse matrices may share this object by holding an
  /// std::shared_ptr or Teuchos::RCP to it.
  WrappedTplMatrixHandle (const WrappedTplMatrixHandle<MatrixType>& rhs);

  /// \brief Forbid assignment (operator=) syntactically.
  ///
  /// Trilinos' sparse matrices may share this object by holding an
  /// std::shared_ptr or Teuchos::RCP to it.
  WrappedTplMatrixHandle<MatrixType>&
  operator= (const WrappedTplMatrixHandle<MatrixType>& rhs);
};

} // namespace Mkl
} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_IMPL_MKL_HPP_
