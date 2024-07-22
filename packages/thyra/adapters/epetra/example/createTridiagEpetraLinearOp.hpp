// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_CREATE_TRIDIAG_EPETRA_LINEAR_OP_HPP
#define THYRA_CREATE_TRIDIAG_EPETRA_LINEAR_OP_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

class Epetra_Operator;

/** \brief \brief This function generates a tridiagonal linear operator using Epetra.
 *
 * Specifically, this function returns a smart pointer to the matrix:
\f[

A=
\left[\begin{array}{rrrrrrrrrr}
2 a    & -1 \\
-1     &  2 a    & -1 \\
       & \ddots  & \ddots  & \ddots \\
       &         & -1      & 2 a     & -1 \\
       &         &         &  -1     & 2 a
\end{array}\right]
\f]
 *
 * where <tt>diagScale</tt> is \f$a\f$ and <tt>globalDim</tt> is the
 * glboal dimension of the matrix.
 */
Teuchos::RCP<Epetra_Operator>
createTridiagEpetraLinearOp(
  const int      globalDim
#ifdef HAVE_MPI
  ,MPI_Comm      mpiComm
#endif
  ,const double  diagScale
  ,const bool    verbose
  ,std::ostream  &out
  );

#endif // THYRA_CREATE_TRIDIAG_EPETRA_LINEAR_OP_HPP
