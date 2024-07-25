// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __SHYLU_FASTUTIL_HPP__
#define __SHYLU_FASTUTIL_HPP__

namespace FastILU {

  enum class SpTRSV {
    Standard,
    StandardHost,
    Fast
  };

}

#endif

