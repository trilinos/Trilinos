// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_PARTITION_HPP__
#define __TACHO_PARTITION_HPP__

/// \file Tacho_Partition.hpp
/// \brief Matrix partitioning utilities.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Part_2x2(const MatView A, MatView &ATL, MatView &ATR,
                                     /**************/ MatView &ABL, MatView &ABR, const ordinal_type bm,
                                     const ordinal_type bn, const int quadrant) {
  ordinal_type bmm, bnn;

  switch (quadrant) {
  case Partition::TopLeft:
    bmm = min(bm, A.extent(0));
    bnn = min(bn, A.extent(1));

    ATL.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  case Partition::TopRight:
  case Partition::BottomLeft:
    TACHO_TEST_FOR_ABORT(true, MSG_NOT_IMPLEMENTED);
    break;
  case Partition::BottomRight:
    bmm = A.extent(0) - min(bm, A.extent(0));
    bnn = A.extent(1) - min(bn, A.extent(1));

    ATL.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }

  ATR.set_view(A, A.offset_0(), ATL.extent(0), A.offset_1() + ATL.extent(1), A.extent(1) - ATL.extent(1));

  ABL.set_view(A, A.offset_0() + ATL.extent(0), A.extent(0) - ATL.extent(0), A.offset_1(), ATL.extent(1));

  ABR.set_view(A, A.offset_0() + ATL.extent(0), A.extent(0) - ATL.extent(0), A.offset_1() + ATL.extent(1),
               A.extent(1) - ATL.extent(1));
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Part_1x2(const MatView A, MatView &AL, MatView &AR, const ordinal_type bn, const int side) {
  ordinal_type bmm, bnn;

  switch (side) {
  case Partition::Left:
    bmm = A.extent(0);
    bnn = min(bn, A.extent(1));

    AL.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  case Partition::Right:
    bmm = A.extent(0);
    bnn = A.extent(1) - min(bn, A.extent(1));

    AL.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }

  AR.set_view(A, A.offset_0(), A.extent(0), A.offset_1() + AL.extent(1), A.extent(1) - AL.extent(1));
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Part_2x1(const MatView A, MatView &AT,
                                     /*************/ MatView &AB, const ordinal_type bm, const int side) {
  ordinal_type bmm, bnn;

  switch (side) {
  case Partition::Top:
    bmm = min(bm, A.extent(0));
    bnn = A.extent(1);

    AT.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  case Partition::Bottom:
    bmm = A.extent(0) - min(bm, A.extent(0));
    bnn = A.extent(1);

    AT.set_view(A, A.offset_0(), bmm, A.offset_1(), bnn);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }

  AB.set_view(A, A.offset_0() + AT.extent(0), A.extent(0) - AT.extent(0), A.offset_1(), A.extent(1));
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void
Part_2x2_to_3x3(const MatView ATL, const MatView ATR, MatView &A00, MatView &A01, MatView &A02,
                /***********************************/ MatView &A10, MatView &A11, MatView &A12, const MatView ABL,
                const MatView ABR, MatView &A20, MatView &A21, MatView &A22, const ordinal_type bm,
                const ordinal_type bn, const int quadrant) {
  switch (quadrant) {
  case Partition::TopLeft:
    Part_2x2(ATL, A00, A01,
             /**/ A10, A11, bm, bn, Partition::BottomRight);

    Part_2x1(ATR, A02,
             /**/ A12, bm, Partition::Bottom);

    Part_1x2(ABL, A20, A21, bn, Partition::Right);

    A22.set_view(ABR, ABR.offset_0(), ABR.extent(0), ABR.offset_1(), ABR.extent(1));
    break;
  case Partition::TopRight:
  case Partition::BottomLeft:
    TACHO_TEST_FOR_ABORT(true, MSG_NOT_IMPLEMENTED);
    break;
  case Partition::BottomRight:
    A00.set_view(ATL, ATL.offset_0(), ATL.extent(0), ATL.offset_1(), ATL.extent(1));

    Part_1x2(ATR, A01, A02, bn, Partition::Left);

    Part_2x1(ABL, A10,
             /**/ A20, bm, Partition::Top);

    Part_2x2(ABR, A11, A12,
             /**/ A21, A22, bm, bn, Partition::TopLeft);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Part_2x1_to_3x1(const MatView AT, MatView &A0,
                                            /***************/ MatView &A1, const MatView AB, MatView &A2,
                                            const ordinal_type bm, const int side) {
  switch (side) {
  case Partition::Top:
    Part_2x1(AT, A0,
             /**/ A1, bm, Partition::Bottom);

    A2.set_view(AB, AB.offset_0(), AB.extent(0), AB.offset_1(), AB.extent(1));
    break;
  case Partition::Bottom:
    A0.set_view(AT, AT.offset_0(), AT.extent(0), AT.offset_1(), AT.extent(1));

    Part_2x1(AB, A1,
             /**/ A2, bm, Partition::Top);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Part_1x2_to_1x3(const MatView AL, const MatView AR, MatView &A0, MatView &A1, MatView &A2,
                                            const ordinal_type bn, const int side) {
  switch (side) {
  case Partition::Left:
    Part_1x2(AL, A0, A1, bn, Partition::Right);

    A2.set_view(AR.BaseObaject(), AR.offset_0(), AR.extent(0), AR.offset_1(), AR.extent(1));
    break;
  case Partition::Right:
    A0.set_view(AL, AL.offset_0(), AL.extent(0), AL.offset_1(), AL.extent(1));

    Part_1x2(AR, A1, A2, bn, Partition::Left);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Merge_2x2(const MatView ATL, const MatView ATR, const MatView ABL, const MatView ABR,
                                      MatView &A) {
  A.set_view(ATL, ATL.offset_0(), ATL.extent(0) + ABR.extent(0), ATL.offset_1(), ATL.extent(1) + ABR.extent(1));
}

template <typename MatView> KOKKOS_INLINE_FUNCTION void Merge_1x2(const MatView AL, const MatView AR, MatView &A) {
  A.set_view(AL, AL.offset_0(), AL.extent(0), AL.offset_1(), AL.extent(1) + AR.extent(1));
}

template <typename MatView> KOKKOS_INLINE_FUNCTION void Merge_2x1(const MatView AT, const MatView AB, MatView &A) {
  A.set_view(AT, AT.offset_0(), AT.extent(0) + AB.extent(0), AT.offset_1(), AT.extent(1));
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Merge_3x3_to_2x2(const MatView A00, const MatView A01, const MatView A02, MatView &ATL,
                                             MatView &ATR, const MatView A10, const MatView A11, const MatView A12,
                                             const MatView A20, const MatView A21, const MatView A22, MatView &ABL,
                                             MatView &ABR, const int quadrant) {
  switch (quadrant) {
  case Partition::TopLeft:
    Merge_2x2(A00, A01, A10, A11, ATL);

    Merge_2x1(A02, A12, ATR);

    Merge_1x2(A20, A21, ABL);

    ABR.set_view(A22, A22.offset_0(), A22.extent(0), A22.offset_1(), A22.extent(1));
    break;
  case Partition::TopRight:
  case Partition::BottomLeft:
    TACHO_TEST_FOR_ABORT(true, MSG_NOT_IMPLEMENTED);
    break;
  case Partition::BottomRight:
    ATL.set_view(A00, A00.offset_0(), A00.extent(0), A00.offset_1(), A00.extent(1));

    Merge_1x2(A01, A02, ATR);

    Merge_2x1(A10, A20, ABL);

    Merge_2x2(A11, A12, A21, A22, ABR);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Merge_3x1_to_2x1(const MatView A0, MatView &AT, const MatView A1, const MatView A2,
                                             MatView &AB, const int side) {
  switch (side) {
  case Partition::Top:
    Merge_2x1(A0, A1, AT);

    AB.set_view(A2, A2.offset_0(), A2.extent(0), A2.offset_1(), A2.extent(1));
    break;
  case Partition::Bottom:
    AT.set_view(A0, A0.offset_0(), A0.extent(0), A0.offset_1(), A0.extent(1));

    Merge_2x1(A1, A2, AB);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

template <typename MatView>
KOKKOS_INLINE_FUNCTION void Merge_1x3_to_1x2(const MatView A0, const MatView A1, const MatView A2, MatView &AL,
                                             MatView &AR, const int side) {
  switch (side) {
  case Partition::Left:
    Merge_1x2(A0, A1, AL);

    AR.set_view(A2, A2.offset_0(), A2.extent(0), A2.offset_1(), A2.extent(1));
    break;
  case Partition::Right:
    AL.set_view(A0, A0.offset_0(), A0.extent(0), A0.offset_1(), A0.extent(1));

    Merge_1x2(A1, A2, AR);
    break;
  default:
    TACHO_TEST_FOR_ABORT(true, MSG_INVALID_INPUT);
    break;
  }
}

} // namespace Tacho

#endif
