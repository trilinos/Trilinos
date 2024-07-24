// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZPRECONDITIONERS_FWD_HPP
#define _FROSCH_SCHWARZPRECONDITIONERS_FWD_HPP

namespace FROSch {

    template <class SC,
              class LO,
              class GO,
              class NO>
    class AlgebraicOverlappingPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class GDSWPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class OneLevelPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class RGDSWPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class SchwarzPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class TwoLevelBlockPreconditioner;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class TwoLevelPreconditioner;

}

#endif
