// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_SUPERLUDIST_CPP
#define AMESOS2_SUPERLUDIST_CPP

#include "Amesos2_Superludist_decl.hpp"


#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#  include "Amesos2_Superludist_def.hpp"
#  include "Amesos2_ExplicitInstantiationHelpers.hpp"

namespace Amesos2 {
#ifdef HAVE_AMESOS2_EPETRA
  AMESOS2_SOLVER_EPETRA_INST(Superludist);
#endif

#ifdef HAVE_TPETRA_INST_INT_INT
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,int);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,int);
#endif
#endif

#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,unsigned);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,unsigned);
#endif
#endif

#ifdef HAVE_TPETRA_INST_INT_LONG
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,long);
#endif
#endif

#ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,unsigned long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,unsigned long);
#endif
#endif


#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,long long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,long long);
#endif
#endif


/*
 * #ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superludist,float,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superludist,float,int,unsigned int);
#endif
#endif

#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superludist,double,int,unsigned int);
#endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<float>,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<float>,int,unsigned int);
#endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superludist,std::complex<double>,int,unsigned int);
#endif
#endif
  */


}

#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include "TpetraCore_ETIHelperMacros.h"

/*#define AMESOS2_SUPERLUDIST_LOCAL_INSTANT(S,LO,GO,N)                        \
  template class Amesos2::Superludist<Tpetra::CrsMatrix<S, LO, GO, N>,      \
                                  Tpetra::MultiVector<S, LO, GO,  N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

#define NODETYPE Tpetra::Details::DefaultTypes::node_type
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SUPERLUDIST_LOCAL_INSTANT(float, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(float, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(float, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(float, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<float>, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<float>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<float>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<float>, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
#endif
#undef NODETYPE*/

#define AMESOS2_SUPERLUDIST_LOCAL_INSTANT(S,LO,GO,N)                        \
  template class Amesos2::Superludist<Tpetra::CrsMatrix<S, LO, GO, N>,      \
                                  Tpetra::MultiVector<S, LO, GO,  N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

#if defined(HAVE_TPETRA_INST_SERIAL) && !defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE)
#define NODETYPE Tpetra_KokkosCompat_KokkosSerialWrapperNode
#ifdef HAVE_TPETRA_INST_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned long, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned long, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_PTHREAD) && !defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE)
#define NODETYPE Tpetra_KokkosCompat_KokkosThreadsWrapperNode
#ifdef HAVE_TPETRA_INST_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned long, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned long, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_OPENMP) && !defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE)
#define NODETYPE Tpetra_KokkosCompat_KokkosOpenMPWrapperNode
#ifdef HAVE_TPETRA_INST_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned long, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    #ifdef HAVE_TPETRA_INST_INT_INT
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
      AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned long, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_CUDA) && !defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE)
#define NODETYPE Tpetra_KokkosCompat_KokkosCudaWrapperNode
#ifdef HAVE_TPETRA_INST_DOUBLE
  #ifdef HAVE_TPETRA_INST_INT_INT
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, int, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(double, int, unsigned long, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  #ifdef HAVE_TPETRA_INST_INT_INT
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED_LONG
    AMESOS2_SUPERLUDIST_LOCAL_INSTANT(std::complex<double>, int, unsigned long, NODETYPE)
  #endif
#endif
#undef NODETYPE
#endif

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#endif  // AMESOS2_SUPERLUDIST_CPP
