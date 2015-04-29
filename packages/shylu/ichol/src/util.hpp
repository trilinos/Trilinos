#pragma once
#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <stdio.h>
#include <string.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>

#include <cmath>
#include <complex>

/// \file util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This provides utility functions for implementing mini-app for incomplete
/// sparse matrix factorization with task-data parallelism e.g., parameter
/// classes, error handling, ostream << overloading.
///
/// Note: The reference of the "static const int" members in the enum-like 
/// classes should not be used as function arguments but their values only. 


using namespace std;

namespace Example {
  
#undef CHKERR
#define CHKERR(ierr)                                                    \
  if (ierr != 0) { cout << endl << ">> Error in " << __FILE__ << ", " << __LINE__ << endl; }

#define MSG_NOT_YET_IMPLEMENTED ">> Not yet implemented"
#define MSG_INVALID_INPUT(what) ">> Invaid input argument: " #what
#define MSG_INVALID_TEMPLATE_ARGS ">> Invaid template arguments"
#define ERROR(msg)                                                      \
  { cout << endl << ">> Error in " << __FILE__ << ", " << __LINE__ << endl << msg << endl; }
  
  /// \class Partition
  /// \brief Matrix partition parameters.
  class Partition { 
  public:
    static const int Top         = 101;
    static const int Bottom      = 102;

    static const int Left        = 201;
    static const int Right       = 202;

    static const int TopLeft     = 401;
    static const int TopRight    = 402;
    static const int BottomLeft  = 403;
    static const int BottomRight = 404;
  };

  /// \class Uplo
  /// \brief Matrix upper/lower parameters.
  class Uplo {
  public:
    static const int Upper = 501;
    static const int Lower = 502;
  };

  /// \class Side
  /// \brief Matrix left/right parameters.
  class Side {
  public:
    static const int Left  = 601;
    static const int Right = 602;
  };

  /// \class Diag
  /// \brief Matrix unit/non-unit diag parameters.
  class Diag {
  public:
    static const int Unit    = 701;
    static const int NonUnit = 702;
  };

  /// \class Trans
  /// \brief Matrix upper/lower parameters.
  class Trans {
  public:
    static const int Transpose     = 801;
    static const int ConjTranspose = 802;
    static const int NoTranspose   = 803;
  };

  /// \class Loop
  /// \brief outer/innner parameters
  class Loop {
    static const int Outer = 901;
    static const int Inner = 902;
    static const int Fused = 903;
  };

  /// \class AlgoChol
  /// \brief Algorithmic variants in sparse factorization and sparse BLAS operations. 
  class AlgoChol {
  public:
    // One side factorization on flat matrices
    static const int LeftUnblocked      = 1001;
    static const int LeftUnblockedOpt1  = 1002;
    static const int LeftBlocked        = 1101;

    static const int RightUnblocked     = 1201;
    static const int RightUnblockedOpt1 = 1202;
    static const int RightBlocked       = 1301;

    // One side factorization on hier matrices
    static const int LeftByBlocks       = 1401;
    static const int RightByBlocks      = 1501;
  };

  // aliasing name space
  typedef AlgoChol AlgoIChol;

  class AlgoGemm {
  public:
    // One side factorization on flat matrices
    static const int ForFactorLeftBlocked      = 2001;
    static const int ForFactorRightBlocked     = 2002;

    // B and C are dense matrix and used for solve phase
    static const int ForSolveDenseMatrix = 2003;
  };

  typedef AlgoGemm AlgoHerk;
  typedef AlgoGemm AlgoTrsm;

  /// \brief Interface for overloaded stream operators.
  template<typename T> 
  inline 
  ostream& operator<<(ostream &os, const auto_ptr<T> &p) {
    return p->showMe(os);
  }

  /// \class Disp
  /// \brief Interface for the stream operator.
  class Disp {
    friend ostream& operator<<(ostream &os, const Disp &disp);
  public:
    Disp() { }
    virtual ostream& showMe(ostream &os) const {
      return os;
    }
  };

  /// \brief Implementation of the overloaded stream operator.
  inline 
  ostream& operator<<(ostream &os, const Disp &disp) {
    return disp.showMe(os);
  }  

}

#endif
