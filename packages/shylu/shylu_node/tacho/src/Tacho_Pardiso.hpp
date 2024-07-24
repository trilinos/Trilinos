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
#ifndef __TACHO_EXAMPLE_PARDISO_HPP__
#define __TACHO_EXAMPLE_PARDISO_HPP__

#if defined(__INTEL_MKL__)
using namespace std;

#include "Tacho_Util.hpp"
#include "mkl_pardiso.h"

namespace Tacho {

class Pardiso {
public:
  enum Phase {
    Analyze = 11,
    AnalyzeFactorize = 12,
    AnalyzeFactorizeSolve = 13,
    Factorize = 22,
    FactorizeSolve = 23,
    Solve = 33,
    ReleaseInternal = 0, // release internal memory for LU matrix number MNUM
    ReleaseAll = -1      // release all internal memory for all matrices
  };

private:
  int _mtype, // matrix type 2 - spd, 4 - hpd
      _phase;

  int _maxfct, // maxfct - maximum number of factors (=1)
      _mnum,   // mnum   - actual matrix for the solution phase (=1)
      _msglvl; // msglvl - print out level 0: nothing, 1: statistics

  // parameters
  int _iparm[64];
  double _dparm[64];

  // internal data address pointers: NEVER TOUCH
  void *_pt[64];

  // Translate Fortran index to C-index
  int Fort(const int i) const { return i - 1; }

  // Matrix type
  //  1 - real structrue sym
  //  2 - real sym pos def
  // -2 - real sym indef
  //  3 - complex structrue sym
  //  4 - complex hermitian pos def
  // -4 - complex and hermitian indef
  //  6 - complex and sym
  // 11 - real and nonsym
  // 13 - complex and nonsym
  template <typename ValueType, int Algo = 2> void setMatrixType();

public:
  Pardiso() : _mtype(0), _phase(0) {}

  void setDefaultParameters() {
    // initialize arrays
    for (int i = 0; i < 64; ++i) {
      _iparm[i] = 0;
      _dparm[i] = 0.0;
      _pt[i] = 0;
    }

    _maxfct = 1;
    _mnum = 1;
    _msglvl = 1;
  }

  void setParameter(const int id, const int value) { _iparm[Fort(id)] = value; }
  void setParameter(const int id, const double value) { _dparm[Fort(id)] = value; }

  template <typename ValueType, int Algo> int init() {
    int ierr = 0;

    setMatrixType<ValueType, Algo>();
    setDefaultParameters();

    // load default param
    pardisoinit(_pt, &_mtype, _iparm);

    // overload default parameters

    // setParameter( 1, 1); // default param: 0 - default, 1 - user provided
    // setParameter( 2, 0); // reordering: 0 - mindegreem, 2 - nd, 3 - parallel nd
    setParameter(27, 1); // mat check: 0 - no, 1 - check
    setParameter(35, 1); // row and col index: 0 - fortran, 1 - CXX

    return ierr;
  }

  ostream &showErrorCode(ostream &os) const {
    os << "   0 No error" << std::endl
       << "-  1 Input inconsistent" << std::endl
       << "-  2 Not enough memory" << std::endl
       << "-  3 Reordering problem" << std::endl
       << "-  4 Zero pivot" << std::endl
       << "-  5 Unclassified (internal) error" << std::endl
       << "-  6 Preordering fail (matrix types 11,13)" << std::endl
       << "-  7 Diagonal matrix problem" << std::endl
       << "-  8 32-bit integer overflow" << std::endl
       << "- 10 No license file pardiso.lic found" << std::endl
       << "- 11 License expired " << std::endl
       << "- 12 Wrong user name or host" << std::endl
       << "-100 over, Krylov fail" << std::endl;
    return os;
  }

  ostream &showStat(ostream &os, const Phase phase) const {
    switch (phase) {
    case Analyze:
      os << "- Phase: Analyze -" << std::endl
         << "Number of perturbed pivots          = " << _iparm[Fort(14)] << std::endl
         << "Number of peak memory symbolic      = " << _iparm[Fort(15)] << std::endl
         << "Number of permenant memory symbolic = " << _iparm[Fort(16)] << " KB " << std::endl;
      break;
    case Factorize:
      os << "- Phase: Factorize -" << std::endl
         << "Peak memory used in factorization   = " << max(_iparm[Fort(15)], _iparm[Fort(16)] + _iparm[Fort(17)])
         << " KB " << std::endl
         << "Memory numerical factorization      = " << _iparm[Fort(17)] << " KB " << std::endl
         << "Number of nonzeros in factors       = " << _iparm[Fort(18)] << std::endl
         << "Number of factorization MFLOP       = " << _iparm[Fort(19)] << std::endl
         << "MFLOPs                              = " << _iparm[Fort(19)] << std::endl;
      break;
    case Solve:
      os << "- Phase: Solve -" << std::endl << "Number of iterative refinements     = " << _iparm[Fort(7)] << std::endl;
      break;
    case AnalyzeFactorize:
      showStat(os, Analyze);
      showStat(os, Factorize);
      break;
    case AnalyzeFactorizeSolve:
      showStat(os, Analyze);
      showStat(os, Factorize);
      showStat(os, Solve);
      break;
    case FactorizeSolve:
      showStat(os, Factorize);
      showStat(os, Solve);
      break;
    default:
      os << "- Phase: " << phase << " -" << std::endl << "Nothing serious in this phase" << std::endl;
      break;
    }
    return os;
  }

private:
  int _ndof, *_ia, *_ja, *_pivot, _nrhs;
  double *_a, *_b, *_x;

public:
  void setProblem(const int ndof, double *a, int *ia, int *ja, int *pivot, const int nrhs, double *b, double *x) {
    _ndof = ndof;
    _a = a;
    _ia = ia;
    _ja = ja;
    _pivot = pivot;
    _nrhs = nrhs;
    _b = b;
    _x = x;
  }

  int run(int phase, int msglvl = 1) {
    int ierr = 0;

    _msglvl = msglvl;
    pardiso(_pt, &_maxfct, &_mnum, &_mtype, &phase, &_ndof, _a, _ia, _ja, _pivot, &_nrhs, _iparm, &_msglvl, _b, _x,
            &ierr);

    return ierr;
  }
};

// Pardiso mtype
//  1 - real    structure sym
//  2 - real    sym posdef
// -2 - real    sym indef
//  3 - complex structure sym
//  4 - complex her posdef
// -4 - complex her indef
//  6 - complex structure sym
// 11 - real    non sym
// 13 - complex non sym
//
// pardiso does not like 2 4; for not we use 1 and 3 or 11 and 13.

template <> void Pardiso::setMatrixType<float, 2>() {
  _mtype = 2;
  setParameter(28, 1);
} // SPD

template <> void Pardiso::setMatrixType<double, 2>() {
  _mtype = 2;
  setParameter(28, 0);
} // SPD

template <> void Pardiso::setMatrixType<complex<float>, 2>() {
  _mtype = 3;
  setParameter(28, 1);
} // HPD

template <> void Pardiso::setMatrixType<complex<double>, 2>() {
  _mtype = 3;
  setParameter(28, 0);
} // HPD
} // namespace Tacho

#endif
#endif
