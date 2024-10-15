// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef _EPETRA_SLU_H_
#define _EPETRA_SLU_H_

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_CrsGraph;

class SLUData;

//! Epetra_SLU:  An object-oriented wrapper for Xiaoye Li's serial sparse solver package:  Superlu.
/*!  Epetra_SLU will solve a linear systems of equations: \f$ AX=B
  \f$, using Epetra objects and the Superludist solver library, where
  \f$A\f$ is an Epetra_RowMatrix and \f$X\f$ and \f$B\f$ are 
  Epetra_MultiVector objects.

  SuperLU execution can be tuned through a variety of parameters.
  Three parameters can be controlled within the call to the constructor:  fill_fac, panel_size and relax. 
  Seven parameters can be controlled within the call to Solve():   
  Verbose, Equil, Factor, perm_type, pivot_thresh, Refact, Trans

*/
class Epetra_SLU
{

 public:


  //@{ \name Constructor methods
  //! Epetra_SLU Constructor.
  /*! Creates an Epetra_SLU instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object.
  */
  Epetra_SLU( Epetra_LinearProblem * Problem,
              int fill_fac = -1,
              int panel_size = -1,
              int relax = -1 );

  //! Epetra_SLU Destructor.
  /*! Completely deletes a Epetra_SLU object.  
  */
  ~Epetra_SLU();
  //@}


  //@{ \name Solve method
  //!  All computation is performed during the call to Solve() 
  /*!  Factor controls whether or not the matrix should be factored prior to the solve.
       Default is true.
   */
  int Solve( bool Verbose = false,
             bool Equil = true,
             bool Factor = true,
             int perm_type = 2,
             double pivot_thresh = -1,
             bool Refact = true,
             bool Trans = false );
  //@}

 private:

  void Copy();

  Epetra_CrsMatrix * A_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  double * R_;
  double * C_;

  Epetra_CrsGraph * AG_;

  SLUData * data_;

  int * Indices_;
  double * Values_;

  int * TransNumNZ_;
  int TotNumNZ_;

  int ** TransIndices_;
  double ** TransValues_;

  int * RowIndices_;
  double * RowValues_;
  int * ColPointers_;

  int * perm_r_;
  int * perm_c_;

  int count_;

  int * etree_;

  double * ferr_;
  double * berr_;
};

#endif
