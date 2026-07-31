// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Each tu_*.cc translation unit in this directory includes exactly one
// public MiniTensor header and nothing else. Compiling and linking them
// into this executable proves every public header is self-contained.
// Reaching main means compilation succeeded, so the test passes.
int
main()
{
  std::cout << "Self-containedness check: End Result: TEST PASSED" << std::endl;
  return 0;
}
