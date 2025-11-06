// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GETRF_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_GETRF_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Iamax.hpp>
#include <KokkosBatched_Laswp.hpp>

namespace KokkosBatched {
namespace Impl {

struct Stack {
 private:
  constexpr static int STACK_SIZE = 48;

  // (state, m_start, n_start, piv_start, m_size, n_size, piv_size)
  int m_stack[7][STACK_SIZE];
  int m_top;

 public:
  KOKKOS_FUNCTION
  Stack() : m_top(-1) {}  // Initialize top to -1, indicating the stack is empty

  KOKKOS_INLINE_FUNCTION
  void push(int values[]) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (m_top >= STACK_SIZE - 1) {
      Kokkos::printf("Stack overflow: Cannot push, the stack is full.\n");
      return;
    }
#endif
    ++m_top;
    for (int i = 0; i < 7; i++) {
      // Increment top and add value
      m_stack[i][m_top] = values[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void pop(int values[]) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (m_top < 0) {
      // Check if the stack is empty
      Kokkos::printf("Stack underflow: Cannot pop, the stack is empty.");
      return;
    }
#endif
    for (int i = 0; i < 7; i++) {
      // Return the top value and decrement top
      values[i] = m_stack[i][m_top];
    }
    m_top--;
  }

  KOKKOS_INLINE_FUNCTION
  bool isEmpty() const { return m_top == -1; }
};

// Host only implementation with recursive algorithm
template <typename AlgoType>
struct SerialGetrfInternalHost {
  template <typename AViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &ipiv);
};

template <>
template <typename AViewType, typename PivViewType>
KOKKOS_INLINE_FUNCTION int SerialGetrfInternalHost<Algo::Getrf::Unblocked>::invoke(const AViewType &A,
                                                                                   const PivViewType &ipiv) {
  using ScalarType = typename AViewType::non_const_value_type;

  const int m = A.extent(0), n = A.extent(1);

  // Quick return if possible
  if (m <= 0 || n <= 0) return 0;

  int info = 0;

  // Use unblocked code for one row case
  // Just need to handle ipiv and info
  if (m == 1) {
    ipiv(0) = 0;
    if (A(0, 0) == 0) return 1;

    return 0;
  } else if (n == 1) {
    // Use unblocked code for one column case
    // Compute machine safe minimum
    auto col_A = Kokkos::subview(A, Kokkos::ALL, 0);

    int i   = SerialIamax::invoke(col_A);
    ipiv(0) = i;

    if (A(i, 0) == 0) return 1;

    // Apply the interchange
    if (i != 0) {
      Kokkos::kokkos_swap(A(i, 0), A(0, 0));
    }

    // Compute elements
    const ScalarType alpha          = 1.0 / A(0, 0);
    auto sub_col_A                  = Kokkos::subview(A, Kokkos::pair<int, int>(1, m), 0);
    [[maybe_unused]] auto info_scal = KokkosBlas::SerialScale::invoke(alpha, sub_col_A);

    return 0;
  } else {
    // Use recursive code
    auto n1 = Kokkos::min(m, n) / 2;

    // Factor A0 = [[A00],
    //              [A10]]

    // split A into two submatrices A = [A0, A1]
    auto A0    = Kokkos::subview(A, Kokkos::ALL, Kokkos::pair<int, int>(0, n1));
    auto A1    = Kokkos::subview(A, Kokkos::ALL, Kokkos::pair<int, int>(n1, n));
    auto ipiv0 = Kokkos::subview(ipiv, Kokkos::pair<int, int>(0, n1));
    auto iinfo = invoke(A0, ipiv0);

    if (info == 0 && iinfo > 0) info = iinfo;

    // Apply interchanges to A1 = [[A01],
    //                             [A11]]

    [[maybe_unused]] auto info_laswp = KokkosBatched::SerialLaswp<Direct::Forward>::invoke(ipiv0, A1);

    // split A into four submatrices
    // A = [[A00, A01],
    //      [A10, A11]]
    auto A00 = Kokkos::subview(A, Kokkos::pair<int, int>(0, n1), Kokkos::pair<int, int>(0, n1));
    auto A01 = Kokkos::subview(A, Kokkos::pair<int, int>(0, n1), Kokkos::pair<int, int>(n1, n));
    auto A10 = Kokkos::subview(A, Kokkos::pair<int, int>(n1, m), Kokkos::pair<int, int>(0, n1));
    auto A11 = Kokkos::subview(A, Kokkos::pair<int, int>(n1, m), Kokkos::pair<int, int>(n1, n));

    // Solve A00 * X = A01
    [[maybe_unused]] auto info_trsm = KokkosBatched::SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit,
                                                                Algo::Trsm::Unblocked>::invoke(1.0, A00, A01);

    // Update A11 = A11 - A10 * A01
    [[maybe_unused]] auto info_gemm =
        KokkosBatched::SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(-1.0, A10, A01,
                                                                                                         1.0, A11);

    // Factor A11
    auto ipiv1 = Kokkos::subview(ipiv, Kokkos::pair<int, int>(n1, Kokkos::min(m, n)));
    iinfo      = invoke(A11, ipiv1);

    if (info == 0 && iinfo > 0) info = iinfo + n1;

    // Apply interchanges to A10
    info_laswp = KokkosBatched::SerialLaswp<Direct::Forward>::invoke(ipiv1, A10);

    // Pivot indices
    for (int i = n1; i < Kokkos::min(m, n); i++) {
      ipiv(i) += n1;
    }

    return info;
  }
}

// Device only implementation with recursive algorithm
template <typename AlgoType>
struct SerialGetrfInternalDevice {
  template <typename AViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &ipiv);
};

template <>
template <typename AViewType, typename PivViewType>
KOKKOS_INLINE_FUNCTION int SerialGetrfInternalDevice<Algo::Getrf::Unblocked>::invoke(const AViewType &A,
                                                                                     const PivViewType &ipiv) {
  using ScalarType = typename AViewType::non_const_value_type;

  const int m = A.extent(0), n = A.extent(1), init_piv_size = ipiv.extent(0);

  Stack stack;
  int initial[7] = {0, 0, 0, 0, m, n, init_piv_size};
  stack.push(initial);

  // Quick return if possible
  if (m <= 0 || n <= 0) return 0;

  while (!stack.isEmpty()) {
    // Firstly, make a subview based on the current state
    int current[7];
    stack.pop(current);

    int state = current[0], m_start = current[1], n_start = current[2], piv_start = current[3], m_size = current[4],
        n_size = current[5], piv_size = current[6];

    // Quick return if possible
    if (m_size <= 0 || n_size <= 0) continue;

    auto A_current = Kokkos::subview(A, Kokkos::pair<int, int>(m_start, m_start + m_size),
                                     Kokkos::pair<int, int>(n_start, n_start + n_size));

    auto ipiv_current = Kokkos::subview(ipiv, Kokkos::pair<int, int>(piv_start, piv_start + piv_size));
    auto n1           = Kokkos::min(m_size, n_size) / 2;

    // split A into two submatrices A = [A0, A1]
    auto A0    = Kokkos::subview(A_current, Kokkos::ALL, Kokkos::pair<int, int>(0, n1));
    auto A1    = Kokkos::subview(A_current, Kokkos::ALL, Kokkos::pair<int, int>(n1, n_size));
    auto ipiv0 = Kokkos::subview(ipiv_current, Kokkos::pair<int, int>(0, n1));
    auto ipiv1 = Kokkos::subview(ipiv_current, Kokkos::pair<int, int>(n1, Kokkos::min(m_size, n_size)));

    // split A into four submatrices
    // A = [[A00, A01],
    //      [A10, A11]]
    auto A00 = Kokkos::subview(A_current, Kokkos::pair<int, int>(0, n1), Kokkos::pair<int, int>(0, n1));
    auto A01 = Kokkos::subview(A_current, Kokkos::pair<int, int>(0, n1), Kokkos::pair<int, int>(n1, n_size));
    auto A10 = Kokkos::subview(A_current, Kokkos::pair<int, int>(n1, m_size), Kokkos::pair<int, int>(0, n1));
    auto A11 = Kokkos::subview(A_current, Kokkos::pair<int, int>(n1, m_size), Kokkos::pair<int, int>(n1, n_size));

    if (state == 0) {
      // start state
      if (m_size == 1) {
        ipiv_current(0) = 0;
        if (A_current(0, 0) == 0) return 1;
        continue;
      } else if (n_size == 1) {
        // Use unblocked code for one column case
        // Compute machine safe minimum
        auto col_A = Kokkos::subview(A_current, Kokkos::ALL, 0);

        int i           = SerialIamax::invoke(col_A);
        ipiv_current(0) = i;

        if (A_current(i, 0) == 0) return 1;

        // Apply the interchange
        if (i != 0) {
          Kokkos::kokkos_swap(A_current(i, 0), A_current(0, 0));
        }

        // Compute elements
        const ScalarType alpha          = 1.0 / A_current(0, 0);
        auto sub_col_A                  = Kokkos::subview(A_current, Kokkos::pair<int, int>(1, m_size), 0);
        [[maybe_unused]] auto info_scal = KokkosBlas::SerialScale::invoke(alpha, sub_col_A);
        continue;
      }

      // Push states onto the stack in reverse order of how they are executed
      // in the recursive version
      int after_second[7] = {2, m_start, n_start, piv_start, m_size, n_size, piv_size};
      int second[7]       = {0,
                             m_start + n1,
                             n_start + n1,
                             piv_start + n1,
                             m_size - n1,
                             n_size - n1,
                             static_cast<int>(Kokkos::min(m_size, n_size)) - n1};
      int after_first[7]  = {1, m_start, n_start, piv_start, m_size, n_size, piv_size};
      int first[7]        = {0, m_start, n_start, piv_start, m_size, n1, n1};

      stack.push(after_second);
      stack.push(second);
      stack.push(after_first);
      stack.push(first);

    } else if (state == 1) {
      // after first recursive call
      // Factor A0 = [[A00],
      //              [A10]]

      // Apply interchanges to A1 = [[A01],
      //                             [A11]]
      KokkosBatched::SerialLaswp<Direct::Forward>::invoke(ipiv0, A1);

      // Solve A00 * X = A01
      [[maybe_unused]] auto info_trsm =
          KokkosBatched::SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit,
                                    Algo::Trsm::Unblocked>::invoke(1.0, A00, A01);

      // Update A11 = A11 - A10 * A01
      [[maybe_unused]] auto info_gemm =
          KokkosBatched::SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
              -1.0, A10, A01, 1.0, A11);

    } else if (state == 2) {
      // after second recursive call
      // Apply interchanges to A10
      KokkosBatched::SerialLaswp<Direct::Forward>::invoke(ipiv1, A10);

      // Pivot indices
      for (int i = n1; i < Kokkos::min(m_size, n_size); i++) {
        ipiv_current(i) += n1;
      }
    }
  }
  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GETRF_SERIAL_INTERNAL_HPP_
