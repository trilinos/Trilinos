// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GlobalMPISession_HPP
#define ROL_GlobalMPISession_HPP


namespace ROL {

struct GlobalMPISession {
  GlobalMPISession(...) {}
  static void abort() {}
  static void barrier() {}
  static int getNProc() { return 1; }
  static int getRank() { return 1; }
  static int sum(int localVal) { return localVal; }
};  

} 

#endif // ROL_GlobalMPISession_HPP

