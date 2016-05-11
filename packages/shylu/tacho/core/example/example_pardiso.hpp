#pragma once
#ifndef __EXAMPLE_PARDISO_HPP__
#define __EXAMPLE_PARDISO_HPP__

#ifdef HAVE_SHYLUTACHO_MKL
using namespace std;

#include "util.hpp"
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
      ReleaseInternal = 0,  // release internal memory for LU matrix number MNUM
      ReleaseAll = -1 // release all internal memory for all matrices
    };

  private:
    int _mtype,  // matrix type 2 - spd, 4 - hpd
      _phase;

    int _maxfct, // maxfct - maximum number of factors (=1)
      _mnum,     // mnum   - actual matrix for the solution phase (=1)
      _msglvl;   // msglvl - print out level 0: nothing, 1: statistics

    // parameters
    int _iparm[64];
    double _dparm[64];

    // internal data address pointers: NEVER TOUCH
    void *_pt[64];

    // Translate Fortran index to C-index
    int Fort(const int i) const { return i-1; }

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
    template<typename ValueType, int Algo>
    void setMatrixType();

  public:
    Pardiso()
      : _mtype(0),
        _phase(0) {
    }

    void setDefaultParameters() {
      // initialize arrays
      for (int i=0;i<64;++i) {
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

    template<typename ValueType, int Algo>
    int init() {
      int ierr = 0;

      setMatrixType<ValueType,Algo>();
      setDefaultParameters();

      // load default param
      pardisoinit(_pt,  &_mtype, _iparm);

      // overload default parameters

      // setParameter( 1, 1); // default param: 0 - default, 1 - user provided
      // setParameter( 2, 0); // reordering: 0 - mindegreem, 2 - nd, 3 - parallel nd
      // setParameter(27, 1); // mat check: 0 - no, 1 - check
      setParameter(35, 1); // row and col index: 0 - fortran, 1 - CXX

      return ierr;
    }

    ostream& showErrorCode(ostream &os) const {
      os << "   0 No error" << endl
         << "-  1 Input inconsistent" << endl
         << "-  2 Not enough memory" << endl
         << "-  3 Reordering problem" << endl
         << "-  4 Zero pivot" << endl
         << "-  5 Unclassified (internal) error" << endl
         << "-  6 Preordering fail (matrix types 11,13)" << endl
         << "-  7 Diagonal matrix problem" << endl
         << "-  8 32-bit integer overflow" << endl
         << "- 10 No license file pardiso.lic found" << endl
         << "- 11 License expired " << endl
         << "- 12 Wrong user name or host" << endl
         << "-100 over, Krylov fail" << endl ;

      return os;
    }

    ostream& showStat(ostream &os, const Phase phase) const {
      switch(phase) {
      case Analyze:
        os << "- Phase: Analyze -" << endl
           << "Number of perturbed pivots          = " << _iparm[Fort(14)] << endl
           << "Number of peak memory symbolic      = " << _iparm[Fort(15)] << endl
           << "Number of permenant memory symbolic = " << _iparm[Fort(16)] << endl;
        break;
      case Factorize:
        os << "- Phase: Factorize -" << endl
           << "Memory numerical factorization      = " << _iparm[Fort(17)] << endl
           << "Number of nonzeros in factors       = " << _iparm[Fort(18)] << endl
           << "Number of factorization MFLOP       = " << _iparm[Fort(19)] << endl
           << "MFLOPs                              = " << _iparm[Fort(19)] << endl;
        break;
      case Solve:
        os << "- Phase: Solve -" << endl
           << "Number of iterative refinements     = " << _iparm[Fort(7)] << endl;
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
        os << "- Phase: " << phase << " -" << endl
           << "Nothing serious in this phase" << endl;
        break;
      }

      return os;
    }

  private:
    int _ndof, *_ia, *_ja, *_pivot, _nrhs;
    double *_a, *_b, *_x;

  public:
    void setProblem(const int ndof,
                    double *a,
                    int *ia,
                    int *ja,
                    int *pivot,
                    const int nrhs,
                    double *b,
                    double *x) {
      _ndof = ndof;
      _a = a;
      _ia = ia;
      _ja = ja;
      _pivot = pivot;
      _nrhs = nrhs;
      _b = b;
      _x = x;
    }

    int run(int phase) {
      int ierr = 0;

      pardiso(_pt, 
              &_maxfct, &_mnum, &_mtype,
              &phase,
              &_ndof, 
              _a, _ia, _ja,
              _pivot, &_nrhs,
              _iparm, &_msglvl,
              _b, _x,
              &ierr);

      return ierr;
    }

  };

  template<>
  void Pardiso::setMatrixType<double, AlgoChol::ExternalPardiso>() { _mtype = 2; } // SPD
  
  template<>
  void Pardiso::setMatrixType<complex<double>, AlgoChol::ExternalPardiso>() { _mtype = 4; } // HPD
}

#endif
#endif
