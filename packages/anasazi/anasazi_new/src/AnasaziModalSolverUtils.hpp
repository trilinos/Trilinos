// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_MODAL_SOLVER_UTILS_HPP
#define ANASAZI_MODAL_SOLVER_UTILS_HPP

/*!     \file AnasaziModalSolverUtils.hpp
        \brief Class which provides internal utilities for the Anasazi modal solvers.
*/

/*!    \class Anasazi::ModalSolverUtils
       \brief Anasazi's templated, static class providing utilities for
       the modal solvers.

       This class provides concrete, templated implementations of utilities necessary
       for the modal solvers (Davidson, LOBPCG, ...).  These utilities include
       sorting, orthogonalization, projecting/solving local eigensystems, and sanity
       checking.  These are internal utilties, so the user should not alter this class.

       \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_BLAS.hpp"

namespace Anasazi {

  template<class STYPE, class MV, class OP>
  class ModalSolverUtils 
  {  
  public:
    
    //@{ \name Sorting Methods

    static int sortScalars(int n, STYPE *y, int *perm = 0) const;

    static int sortScalars_Vectors(int n, STYPE* lambda, STYPE* Q = 0, int ldQ = 0) const;

    //@} 

    //@{ \name Eigensolver Projection Methods

    static int massOrthonormalize(MV &X, MV &MX, const OP *M, const MV &Q, int howMany,
				  int type = 0, STYPE *WS = 0, STYPE kappa = 1.5625) const;

    static void localProjection(int numRow, int numCol, int length,
				STYPE *U, int ldU, STYPE *MatV, int ldV,
				STYPE *UtMatV, int ldUtMatV, STYPE *work) const;
    
    static int directSolver(int, STYPE*, int, STYPE*, int, int&, 
			    STYPE*, int, STYPE*, int, int = 0) const;

    //@}

    //@{ \name Sanity Checking Methods

    static STYPE errorOrthogonality(const MV &X, const MV &R, const OP *M = 0) const;
    
    static STYPE errorOrthonormality(const MV &X, const OP &M = 0) const;
    
    static STYPE errorEquality(const MV &X, const MV &MX, const OP &M = 0) const;
    
    static int inputArguments(const int &numEigen, const OP &K, const OP &M, const OP &P,
			      const MV &Q, const int &minSize) const;
    
    //@}
  };

  //-----------------------------------------------------------------------------
  // 
  //  SORTING METHODS
  //
  //-----------------------------------------------------------------------------
  
  template<STYPE, MV, OP>
  static int ModalSolverUtils<STYPE, MV, OP>::sortScalars( int n, STYPE *y, int *perm = 0) const 
  {
  // Sort a vector into increasing order of algebraic values
  //
  // Input:
  //
  // n    (integer ) = Size of the array (input)
  // y    (double* ) = Array of length n to be sorted (input/output)
  // perm (integer*) = Array of length n with the permutation (input/output)
  //                   Optional argument

  int i, j;
  int igap = n / 2;

  if (igap == 0) {
    if ((n > 0) && (perm != 0)) {
      perm[0] = 0;
    }
    return 0;
  }

  if (perm) {
    for (i = 0; i < n; ++i)
      perm[i] = i;
  }

  while (igap > 0) {
    for (i=igap; i<n; ++i) {
      for (j=i-igap; j>=0; j-=igap) {
        if (y[j] > y[j+igap]) {
          double tmpD = y[j];
          y[j] = y[j+igap];
          y[j+igap] = tmpD;
          if (perm) {
            int tmpI = perm[j];
            perm[j] = perm[j+igap];
            perm[j+igap] = tmpI;
          }
        }
        else {
          break;
        }
      }
    }
    igap = igap / 2;
  }

  return 0;
  }

  template<STYPE, MV, OP>
  static int ModalSolverUtils<STYPE, MV, OP>::sortScalars_Vectors( int n, STYPE *lambda, STYPE *Q, int ldQ) const
  {
  // This routines sorts the scalars (stored in lambda) in ascending order.
  // The associated vectors (stored in Q) are accordingly ordered.
  // One vector is of length ldQ.
  // Q must be of size ldQ * num.

    Teuchos::BLAS<int,STYPE> blas;
    STYPE tmp;
    int info = 0;
    int i, j;
    
    int igap = num / 2;
    
    if ((Q) && (ldQ > 0)) {
      while (igap > 0) {
	for (i=igap; i < num; ++i) {
	  for (j=i-igap; j>=0; j-=igap) {
	    if (lambda[j] > lambda[j+igap]) {
	      // Swap two scalars
	      tmp = lambda[j];
	      lambda[j] = lambda[j+igap];
	      lambda[j+igap] = tmp;
	      // Swap corresponding vectors
	      blas.SWAP( ldQ, Q + j*ldQ, 1, Q + (j+igap)*ldQ, 1 );
	    } 
	    else {
	      break;
	    }
	  }
	}
	igap = igap / 2;
      } // while (igap > 0)
    } // if ((Q) && (ldQ > 0))
    else {
      while (igap > 0) {
	for (i=igap; i < num; ++i) {
	  for (j=i-igap; j>=0; j-=igap) {
	    if (lambda[j] > lambda[j+igap]) {
	      // Swap two scalars
	      tmp = lambda[j];
	      lambda[j] = lambda[j+igap];
	      lambda[j+igap] = tmp;
	    } 
	    else {
	      break;
	    }
	  }
	}
	igap = igap / 2;
      } // while (igap > 0)
    } // if ((Q) && (ldQ > 0))
    
  return info;  
  }

  //-----------------------------------------------------------------------------
  // 
  //  EIGENSOLVER PROJECTION METHODS
  //
  //-----------------------------------------------------------------------------

  template<STYPE, MV, OP>
  static int ModalSolverUtils<STYPE, MV, OP>::massOrthonormalize(MV &X, MV &MX, const OP *M, const MV &Q, int howMany,
								 int type, STYPE *WS, STYPE kappa) const
  {
  }
  
  template<STYPE, MV, OP>
  static void ModalSolverUtils<STYPE, MV, OP>::localProjection(int numRow, int numCol, int length,
							       STYPE *U, int ldU, STYPE *MatV, int ldV,
							       STYPE *UtMatV, int ldUtMatV, STYPE *work) const
  {
  }
  
  template<STYPE, MV, OP>
  static int ModalSolverUtils<STYPE, MV, OP>::directSolver(int, STYPE*, int, STYPE*, int, int&, 
							   STYPE*, int, STYPE*, int, int) const
  {
  }
  
  //-----------------------------------------------------------------------------
  // 
  //  SANITY CHECKING METHODS
  //
  //-----------------------------------------------------------------------------

  template<STYPE, MV, OP>
  static STYPE ModalSolverUtils<STYPE, MV, OP>::errorOrthogonality(const MV &X, const MV &R, 
								   const OP *M = 0) const
  {
  }
  
  template<STYPE, MV, OP>
  static STYPE ModalSolverUtils<STYPE, MV, OP>::errorOrthonormality(const MV &X, const OP &M = 0) const
  {
  }
  
  template<STYPE, MV, OP>
  static STYPE ModalSolverUtils<STYPE, MV, OP>::errorEquality(const MV &X, const MV &MX, 
							      const OP &M = 0) const
  {
  }    
  
  template<STYPE, MV, OP>
  static int ModalSolverUtils<STYPE, MV, OP>::inputArguments(const int &numEigen, const OP &K, 
							     const OP &M, const OP &P,
							     const MV &Q, const int &minSize) const
  {
  }  

} // end namespace Anasazi

#endif // ANASAZI_MODAL_SOLVER_UTILS_HPP

