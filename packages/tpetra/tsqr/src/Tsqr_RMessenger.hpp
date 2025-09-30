// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_RMessenger_hpp
#define __TSQR_RMessenger_hpp

#include "Tsqr_MatView.hpp"
#include "Tsqr_MessengerBase.hpp"
#include "Teuchos_RCP.hpp"
#include <algorithm>
#include <vector>

namespace TSQR {

/// \class RMessenger
/// \brief Send, receive, and broadcast square R factors.
///
/// Object that handles sending, receiving, and broadcasting square
/// upper triangular matrices containing data of type Scalar, and
/// indexed by indices of type Ordinal.
template <class Ordinal, class Scalar>
class RMessenger {
 public:
  typedef Scalar scalar_type;
  typedef Ordinal ordinal_type;
  typedef MessengerBase<Scalar> messenger_type;
  typedef Teuchos::RCP<messenger_type> messenger_ptr;

  RMessenger() = delete;

  /// \brief Constructor
  ///
  /// \param messenger [in/out] Pointer to the communicator wrapper.
  RMessenger(const messenger_ptr& messenger)
    : messenger_(messenger) {}

  template <class ConstMatrixViewType>
  void
  send(const ConstMatrixViewType& R, const int destProc) {
    pack(R);
    messenger_->send(&buffer_[0], buffer_.size(), destProc, 0);
  }

  template <class MatrixViewType>
  void
  recv(MatrixViewType& R, const int srcProc) {
    const typename MatrixViewType::ordinal_type ncols = R.extent(1);
    const Ordinal buflen                              = buffer_length(ncols);
    buffer_.resize(buflen);
    messenger_->recv(&buffer_[0], buflen, srcProc, 0);
    unpack(R);
  }

  template <class MatrixViewType>
  void
  broadcast(MatrixViewType& R, const int rootProc) {
    const int myRank = messenger_->rank();
    if (myRank == rootProc)
      pack(R);
    messenger_->broadcast(buffer_.data(), buffer_length(R.extent(1)), rootProc);
    if (myRank != rootProc)
      unpack(R);
  }

  //! Copy constructor
  RMessenger(const RMessenger& rhs)
    : messenger_(rhs.messenger_)
    , buffer_(0)  // don't need to copy the buffer
  {}

  //! Assignment operator
  RMessenger& operator=(const RMessenger& rhs) {
    if (this != &rhs) {
      this->messenger_ = rhs.messenger_;
      // Don't need to do anything to this->buffer_; the various
      // operations such as pack() will resize it as necessary.
    }
    return *this;
  }

 private:
  messenger_ptr messenger_;
  std::vector<Scalar> buffer_;

  /// \brief Buffer length as a function of R factor dimension.
  ///
  /// \param ncols [in] Number of columns (and number of rows)
  ///   in the R factor input.
  Ordinal buffer_length(const Ordinal ncols) const {
    return (ncols * (ncols + Ordinal(1))) / Ordinal(2);
  }

  template <class ConstMatrixViewType>
  void
  pack(const ConstMatrixViewType& R) {
    using view_scalar_type  = typename ConstMatrixViewType::non_const_value_type;
    using view_ordinal_type = typename ConstMatrixViewType::ordinal_type;

    const view_ordinal_type ncols = R.extent(1);
    const Ordinal buf_length      = buffer_length(ncols);
    buffer_.resize(buf_length);
    auto iter = buffer_.begin();
    for (view_ordinal_type j = 0; j < ncols; ++j) {
      const view_scalar_type* const R_j = &R(0, j);
      std::copy(R_j, R_j + (j + 1), iter);
      iter += (j + 1);
    }
  }

  template <class MatrixViewType>
  void
  unpack(MatrixViewType& R) {
    typedef typename MatrixViewType::ordinal_type view_ordinal_type;
    typedef typename std::vector<Scalar>::const_iterator const_iter_type;

    const view_ordinal_type ncols = R.extent(1);
    const_iter_type iter          = buffer_.begin();
    for (view_ordinal_type j = 0; j < ncols; ++j) {
      std::copy(iter, iter + (j + 1), &R(0, j));
      iter += (j + 1);
    }
  }
};

/// \brief Distribute a stack of R factors.
///
/// \param R_stack [in] nprocs*ncols by ncols stack of square upper
///   triangular matrices.  The whole stack is stored in
///   column-major order.
///
/// \param R_local [out] ncols by ncols upper triangular matrix,
///   stored in column-major order (in unpacked form).
///
/// \param messenger [in/out] Object that handles communication
///
template <class MatrixViewType, class ConstMatrixViewType>
void scatterStack(const ConstMatrixViewType& R_stack,
                  MatrixViewType& R_local,
                  const Teuchos::RCP<MessengerBase<typename MatrixViewType::non_const_value_type>>& messenger) {
  using ordinal_type    = typename MatrixViewType::ordinal_type;
  using scalar_type     = typename MatrixViewType::non_const_value_type;
  using const_view_type = MatView<ordinal_type, const scalar_type>;

  const int nprocs  = messenger->size();
  const int my_rank = messenger->rank();

  if (my_rank == 0) {
    const ordinal_type ncols = R_stack.extent(1);
    // Copy data from top ncols x ncols block of R_stack into R_local.
    const_view_type R_stack_view_first(ncols, ncols, R_stack.data(),
                                       R_stack.stride(1));
    deep_copy(R_local, R_stack_view_first);

    // Loop through all other processors, sending each the next
    // ncols x ncols block of R_stack.
    RMessenger<ordinal_type, scalar_type> sender(messenger);
    for (int destProc = 1; destProc < nprocs; ++destProc) {
      auto R_ptr = R_stack.data() + destProc * ncols;
      const_view_type R_stack_view_cur(ncols, ncols, R_ptr, R_stack.stride(1));
      sender.send(R_stack_view_cur, destProc);
    }
  } else {
    const int srcProc = 0;
    deep_copy(R_local, scalar_type{});
    RMessenger<ordinal_type, scalar_type> receiver(messenger);
    receiver.recv(R_local, srcProc);
  }
}

template <class MatrixViewType, class ConstMatrixViewType>
void gatherStack(MatrixViewType& R_stack,
                 ConstMatrixViewType& R_local,
                 const Teuchos::RCP<MessengerBase<typename MatrixViewType::non_const_value_type>>& messenger) {
  using ordinal_type  = typename MatrixViewType::ordinal_type;
  using scalar_type   = typename MatrixViewType::non_const_value_type;
  using mat_view_type = MatView<ordinal_type, scalar_type>;

  const int nprocs  = messenger->size();
  const int my_rank = messenger->rank();

  if (my_rank == 0) {
    const ordinal_type ncols = R_stack.extent(1);
    // Copy data from R_local into top ncols x ncols block of R_stack.
    mat_view_type R_stack_view_first(ncols, ncols, R_stack.data(),
                                     R_stack.stride(1));
    deep_copy(R_stack_view_first, R_local);

    // Loop through all other processors, fetching their matrix data.
    RMessenger<ordinal_type, scalar_type> receiver(messenger);
    for (int srcProc = 1; srcProc < nprocs; ++srcProc) {
      auto R_ptr = R_stack.data() + srcProc * ncols;
      mat_view_type R_stack_view_cur(ncols, ncols, R_ptr,
                                     R_stack.stride(1));
      // Fill (the lower triangle) with zeros, since
      // RMessenger::recv() only writes to the upper triangle.
      deep_copy(R_stack_view_cur, scalar_type{});
      receiver.recv(R_stack_view_cur, srcProc);
    }
  } else {
    // We only read R_stack on Proc 0, not on this proc.
    // Send data from R_local to Proc 0.
    const int destProc = 0;
    RMessenger<ordinal_type, scalar_type> sender(messenger);
    sender.send(R_local, destProc);
  }
  messenger->barrier();
}

}  // namespace TSQR

#endif  // __TSQR_RMessenger_hpp
