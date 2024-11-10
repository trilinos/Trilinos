// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZOPERATORS_FWD_HPP
#define _FROSCH_SCHWARZOPERATORS_FWD_HPP

namespace FROSch {

    template <class SC,
              class LO,
              class GO,
              class NO>
    class AlgebraicOverlappingOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class CoarseOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class GDSWCoarseOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class HarmonicCoarseOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class IPOUHarmonicCoarseOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class MultiplicativeOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class OverlappingOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class RGDSWCoarseOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class SchwarzOperator;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class SumOperator;

}

#endif
