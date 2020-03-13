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
// ************************************************************************
// @HEADER

#include "Tpetra_Details_LocalCuSparseCrsMatrixOperator.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Tpetra_Details_CuSparseVector_fwd.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include "Teuchos_Assert.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosBlas1_scal.hpp"

namespace Tpetra {
namespace Details {

template<class Scalar>
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
LocalCuSparseCrsMatrixOperatorRealBase(
  const execution_space& execSpace,
  const std::shared_ptr<local_matrix_type>& A)
  : base_type(A),
    handle_(getCuSparseHandle(execSpace)),
    matrix_(nullptr, &Impl::deleteCuSparseMatrix)
{}

template<class Scalar>
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
~LocalCuSparseCrsMatrixOperatorRealBase() {}

template<class Scalar>
void
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
resumeFill()
{
  matrix_.reset();
  base_type::resumeFill();
}

template<class Scalar>
void
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
setMinMaxNumberOfEntriesPerRow(const LO minNumEntPerRow,
                               const LO maxNumEntPerRow)
{
  minNumEntPerRow_ = minNumEntPerRow;
  maxNumEntPerRow_ = maxNumEntPerRow;
  base_type::setMinMaxNumberOfEntriesPerRow(minNumEntPerRow,
                                            maxNumEntPerRow);
}

template<class Scalar>
void
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
fillComplete()
{
  using std::endl;
  const char funcName[] = "Tpetra::Details::"
    "LocalCuSparseCrsMatrixOperatorRealBase::fillComplete";

  const bool debug = Details::Behavior::debug("cuSPARSE");
  const bool verbose = Details::Behavior::verbose("cuSPARSE");
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    std::ostringstream os;
    os << funcName << ": ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start" << endl;
    std::cerr << os.str();
  }
  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On entry to "
       << funcName << ", CUDA is in an erroneous state \""
       << cudaGetErrorName(lastErr) << "\".");
  }

  if (! this->isFillComplete()) {
    auto A = this->getLocalMatrix();
    const int64_t numRows = static_cast<int64_t>(A.numRows());
    const int64_t numCols = static_cast<int64_t>(A.numCols());
    const int64_t numEnt = numRows == 0 || numCols == 0 ?
      int64_t(0) : static_cast<int64_t>(A.nnz());

    if (numEnt == 0) {
      matrix_.reset();
    }
    else {
      using LO = DefaultTypes::local_ordinal_type;

      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;

      // NOTE (mfh 09 Mar 2020) I'm assuming that if ptr_ has the
      // right length, then we don't need to copy it again.  This ties
      // in with Tpetra::CrsMatrix's assumption that the matrix's
      // graph structure is fixed after first fillComplete.

      const auto newPtrLen = A.graph.row_map.extent(0);
      if (debug) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (A.graph.entries.extent(0) != A.values.extent(0),
           std::logic_error, "A.graph.entries.extent(0)="
           << A.graph.entries.extent(0) << " != A.values.extent(0)="
           << A.values.extent(0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION
          (newPtrLen == 0 && A.graph.entries.extent(0) != 0,
           std::logic_error, "A.graph.row_map.extent(0)=0, but "
           "A.graph.entries.extent(0)=" << A.graph.entries.extent(0)
           << " != 0.");
      }

      if (newPtrLen != ptr_.extent(0)) {
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "Reallocate ptr_: old=" << ptr_.extent(0)
             << ", new=" << newPtrLen << endl;
          std::cerr << os.str();
        }
        // Free memory before (re)allocating; this may reduce the
        // high-water memory mark.
        ptr_ = Kokkos::View<LO*, device_type>();
        ptr_ = Kokkos::View<LO*, device_type>(
          view_alloc("Tpetra::CrsMatrix cuSPARSE row offsets",
                     WithoutInitializing),
          newPtrLen);
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "Call copyOffsets" << endl;
          std::cerr << os.str();
        }
        Details::copyOffsets(ptr_, A.graph.row_map);
        if (debug) {
          const cudaError_t lastErr = cudaGetLastError();
          TEUCHOS_TEST_FOR_EXCEPTION
            (lastErr != cudaSuccess, std::runtime_error, "In "
             << funcName << ", after calling copyOffsets, CUDA is in "
             "an erroneous state \"" << cudaGetErrorName(lastErr)
             << "\".");
        }
      }

      if (debug && newPtrLen != 0) {
        LO nnz = 0;
        Kokkos::deep_copy(nnz, Kokkos::subview(ptr_, numRows));
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "ptr_[numRows=" << numRows << "]="
             << nnz << endl;
          std::cerr << os.str();
        }
        TEUCHOS_TEST_FOR_EXCEPTION
          (nnz != LO(A.graph.entries.extent(0)), std::logic_error,
           "ptr_[numRows=" << numRows << "]=" << nnz << " != "
           "A.graph.entries.extent(0)=" << A.graph.entries.extent(0)
           << ".");
      }

      LO* ptr = ptr_.data();
      LO* ind = const_cast<LO*>(A.graph.entries.data());
      Scalar* val = const_cast<Scalar*>(A.values.data());

      // FIXME (mfh 09 Mar 2020) Replace this with the logic to pick the
      // right algorithm.  For now, we just want to test merge path.
      const CuSparseMatrixVectorMultiplyAlgorithm alg =
        CuSparseMatrixVectorMultiplyAlgorithm::LOAD_BALANCED;
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Call getCuSparseMatrix" << endl;
        std::cerr << os.str();
      }
      matrix_ = getCuSparseMatrix(numRows, numCols, numEnt,
                                  ptr, ind, val, alg);
    }
    base_type::fillComplete();
  }

  if (debug) {
    cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On exit of "
       << funcName << ", before fence, CUDA is in an erroneous "
       "state \"" << cudaGetErrorName(lastErr) << "\".");
    Kokkos::fence();
    lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On exit of "
       << funcName << ", after fence, CUDA is in an erroneous "
       "state \"" << cudaGetErrorName(lastErr) << "\".");
  }
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str();
  }
}

template<class Scalar>
void
LocalCuSparseCrsMatrixOperatorRealBase<Scalar>::
apply(Kokkos::View<const Scalar**, array_layout,
        device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
      Kokkos::View<Scalar**, array_layout,
        device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
      const Teuchos::ETransp mode,
      const Scalar alpha,
      const Scalar beta) const
{
  using std::endl;
  const char tfecfFuncName[] = "apply";

  const bool debug = Details::Behavior::debug("cuSPARSE");
  const bool verbose = Details::Behavior::verbose("cuSPARSE");
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    std::ostringstream os;
    os << "Tpetra::Details::"
      "LocalCuSparseCrsMatrixOperatorRealBase::apply: ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start" << endl;
    std::cerr << os.str();
  }
  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lastErr != cudaSuccess, std::runtime_error, ": On entry, CUDA "
       "is in an erroneous state \"" << cudaGetErrorName(lastErr)
       << "\".");
  }

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (! this->isFillComplete(), std::runtime_error, ": This object is "
     "not fill complete.  You must call fillComplete on this object "
     "without an intervening resumeFill call, before you may call "
     "apply.");
  // if (! this->isFillComplete()) {
  //   using this_type = LocalCuSparseCrsMatrixOperatorRealBase<Scalar>;
  //   this_type& thisRef = const_cast<this_type&>(*this);
  //   thisRef.fillComplete(); // evil hack
  // }

  using size_type = Kokkos::Cuda::memory_space::size_type;
  const size_type numCols = X.extent(1);
  TEUCHOS_ASSERT( Y.extent(1) == numCols );

  if (matrix_.get() == nullptr) {
    if (beta == Scalar{}) {
      Kokkos::deep_copy(Y, Scalar{});
    }
    else { // beta != 0
      KokkosBlas::scal(Y, beta, Y);
    }
  }
  else {
    for (size_type col = 0; col < numCols; ++col) {
      auto X_col = Kokkos::subview(X, Kokkos::ALL(), col);
      auto Y_col = Kokkos::subview(Y, Kokkos::ALL(), col);
      auto X_col_cuda = getCuSparseVector(
        const_cast<Scalar*>(X_col.data()), X_col.extent(0));
      auto Y_col_cuda = getCuSparseVector(
        Y_col.data(), Y_col.extent(0));
      cuSparseMatrixVectorMultiply(*handle_, mode, alpha, *matrix_,
                                   *X_col_cuda, beta, *Y_col_cuda);
    }
  }

  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On exit of "
       "Tpetra::Details::LocalCuSparseCrsMatrixOperatorRealBase::"
       "apply, CUDA is in an erroneous state \""
       << cudaGetErrorName(lastErr) << "\".");
  }
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str();
  }
}

// full specializations begin here

#ifdef HAVE_TPETRA_INST_FLOAT
LocalCuSparseCrsMatrixOperator<
  float, float,
  Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>>::
LocalCuSparseCrsMatrixOperator(
  const execution_space& execSpace,
  const std::shared_ptr<local_matrix_type>& A)
  : base_type(execSpace, A)
{}
#endif // HAVE_TPETRA_INST_FLOAT

#ifdef HAVE_TPETRA_INST_DOUBLE
LocalCuSparseCrsMatrixOperator<
  double, double,
  Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>>::
LocalCuSparseCrsMatrixOperator(
  const execution_space& execSpace,
  const std::shared_ptr<local_matrix_type>& A)
  : base_type(execSpace, A)
{}
#endif // HAVE_TPETRA_INST_DOUBLE

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE
