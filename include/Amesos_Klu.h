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

/*!
 * \file Amesos_Klu.h
 *
 * \class Amesos_Klu
 *
 * \brief Interface to KLU internal solver.
 *
 * \date Last updated on 24-May-05.
 */

#ifndef AMESOS_KLU_H
#define AMESOS_KLU_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

// class EpetraExt::MultiVector_Reindex ;
// class EpetraExt::CrsMatrix_Reindex ;
//! Amesos_Klu:  A serial, unblocked code ideal for getting started and for very sparse matrices, such as circuit matrces.

/*!

Class Amesos_Klu is an object-oriented wrapper for KLU. KLU, whose sources
are distributed
within Amesos, is a serial solver for sparse matrices. KLU will solve a
linear system of equations: \f$ A X = B\f$, where
<TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are
Epetra_MultiVector objects.

Amesos_Klu computes \f$ A^T X = B\f$
more efficiently than \f$ >A X = B\f$.  The
latter requires a matrix transpose -- which costs both time and space.

KLU is Tim Davis' implementation of Gilbert-Peierl's left-looking
sparse partial pivoting algorithm, with Eisenstat & Liu's symmetric
pruning.  Gilbert's version appears as \c [L,U,P]=lu(A) in MATLAB.
It doesn't exploit dense matrix kernels, but it is the only sparse
LU factorization algorithm known to be asymptotically optimal,
in the sense that it takes time proportional to the number of
floating-point operations.  It is the precursor to SuperLU,
thus the name ("clark Kent LU").  For very sparse matrices that
do not suffer much fill-in (such as most circuit matrices when
permuted properly) dense matrix kernels do not help, and the
asymptotic run-time is of practical importance.

The \c klu_btf code first permutes the matrix to upper block
triangular form (using two algorithms by Duff and Reid,
MC13 and MC21, in the ACM Collected Algorithms).  It then permutes
each block via a symmetric minimum degree ordering (AMD, by Amestoy,
Davis, and Duff).  This ordering phase can be done just once
for a sequence of matrices.  Next, it factorizes each reordered
block via the klu routine, which also attempts to preserve
diagonal pivoting, but allows for partial pivoting if the diagonal
is to small.

*/

// Amesos_Klu_Pimpl contains a pointer to two structures defined in
// klu.h:  trilinos_klu_symbolic and trilinos_klu_numeric.  This prevents Amesos_Klu.h
// from having to include klu.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Klu_Pimpl ;
class Amesos_StandardIndex ;
#endif

class Amesos_Klu: public Amesos_BaseSolver,
                  private Amesos_Time,
                  private Amesos_NoCopiable,
                  private Amesos_Utils,
                  private Amesos_Control,
                  private Amesos_Status {

public:

  //@{ \name Constructors and Destructors
  //! Amesos_Klu Constructor.
  /*! Creates an Amesos_Klu instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object.

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Klu(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Klu Destructor.
  ~Amesos_Klu(void);

  //@}
  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if KLU can handle this matrix shape
  /*! Returns true if the matrix shape is one that KLU can
    handle. KLU only works with square matrices.
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Klu
  /*!
    If SetUseTranspose() is set to true,
    \f$ A^T X = B\f$ is computed.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int SetParameters( Teuchos::ParameterList &ParameterList );

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //! Prints timing information
  void PrintTiming() const;

  //! Prints information about the factorization and solution phases.
  void PrintStatus() const;

  //! Extracts timing information and places in parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming( TimingParameterList ); }

private:

  //@}
  //@{ \name Utility methods

  /*
  CreateLocalMatrixAndExporters - Prepare to convert matrix and vectors to serial
    Preconditions:
      Problem_ must be set

    Postconditions:
      UseDataInPlace_ is set to 1 if the input matrix can be used in place, i.e.
        1)  is entirely stored on process 0
        2)  range map and domain map are same as the row map
      The following are only set if (! UseDataInPlace_ )"
        SerialMap_
        ImportToSerial_
        SerialCrsMatrixA_

      SerialMatrix_
   */
  int CreateLocalMatrixAndExporters() ;
  /*
    ExportToSerial
    Preconditions:
       UseDataInPlace_ must be set
       ImportToSerial and SerialCrsMatrixA_ must be set if UseDataInPlace_ != 1
    Postconditions
       SerialMatrix_ points to a serial version of the matrix
   */
  int ExportToSerial() ;
  /*
    ConvertToKluCRS - Convert matrix to form expected by Klu: Ai, Ap, Aval
    Preconditions:
      numentries_, RowMatrixA_, ImportToSerial_, StdIndexMatrix_, Reindex_
    Postconditions:
      SerialCrsMatrixA_
  */
  int ConvertToKluCRS(bool firsttime);

  /*
    PerformSymbolicFactorization - Call Klu to perform symbolic factorization
    Preconditions:
      UseDataInPlace_ must be set to 1 if the input matrix is entirely stored on process 0
      Ap, Ai and Aval point to a compressed row storage version of the input matrix A.
    Postconditions:
      Symbolic points to an KLU internal opaque object containing the
        symbolic factorization and accompanying information.
      SymbolicFactorizationOK_ = true;
    Note:  All action is performed on process 0

    Returns non-zero if the symbolic factorization failed
  */

  int PerformSymbolicFactorization();

  /*
    PerformNumericFactorization - Call Klu to perform numeric factorization
    Preconditions:
      UseDataInPlace_ must be set
      Ap, Ai and Aval point to a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an KLU internal opaque object containing the
        numeric factorization and accompanying information.
      NumericFactorizationOK_ = true;
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization();

  // @}

  int SerialXlda_ ;

#ifdef Bug_8212
  int *lose_this_;
#endif
  //
  //  PrivateKluData_ contains pointers to data needed by klu whose
  //  data structures are defined by klu.h
  //
  Teuchos::RCP<Amesos_Klu_Pimpl> PrivateKluData_;
  Teuchos::RCP<Amesos_StandardIndex> StdIndex_;
  Teuchos::RCP<Amesos_StandardIndex> StdIndexRange_;
  Teuchos::RCP<Amesos_StandardIndex> StdIndexDomain_;

  //! Ap, Ai, Aval form the compressed row storage used by Klu
  //! Ai and Aval can point directly into a matrix if it is StorageOptimized(), hence
  //! they may either be in vector form or may be a pointer into Epetra_CrsMatrix
  //! internals.  Ap must always be constructed.
  std::vector <int> Ap;
  std::vector <int> VecAi;
  std::vector <double> VecAval;
  double* Aval;
  int *Ai;

  //! 1 if Problem_->GetOperator() is stored entirely on process 0
  int UseDataInPlace_;
  //! Number of non-zero entries in Problem_->GetOperator()
  long long numentries_;
  //! Number of rows and columns in the Problem_->GetOperator()
  long long NumGlobalElements_;

  //! Operator converted to a RowMatrix
  Epetra_RowMatrix* RowMatrixA_;
  //! Operator converted to a CrsMatrix
  Epetra_CrsMatrix* CrsMatrixA_;
#if 0
  //! Points to an object which reindexes a MultiVector to a contiguous map
  Teuchos::RCP<EpetraExt::MultiVector_Reindex> VecTrans_;
  //! Points to an object which reindexes a CrsMatrix to a contiguous map
  Teuchos::RCP<EpetraExt::CrsMatrix_Reindex> MatTrans_;
  //! Points to a Contiguous Map
  Teuchos::RCP<Epetra_Map> ContiguousMap_;
#endif
  //! Points to a Serial Map (unused if UseDataInPlace_ == 1 )
  Teuchos::RCP<Epetra_Map> SerialMap_;
  //! Points to a Serial Copy of A (unused if UseDataInPlace_==1)
  Teuchos::RCP<Epetra_CrsMatrix> SerialCrsMatrixA_;
  //! Points to a Contiguous Copy of A
  Epetra_RowMatrix* StdIndexMatrix_ ;
  //! Points to a Serial Copy of A
  Epetra_RowMatrix* SerialMatrix_ ;

  //! If \c true, no checks are made and the matrix is assume to be distributed
  //  serially, StorageOptimized, the LHS and RHS are assumed to be available
  //  when SymbolicFactorization is called and not to change (address or number
  //  of vectors) thereafter.
  bool TrustMe_;
  //! Number of vectors in RHS and LHS
  int NumVectors_;
  //! Pointer to the actual values in the serial version of X and B
  double *SerialXBvalues_ ;
  double *SerialBvalues_ ;
  //! Serial versions of the LHS and RHS (may point to the original vector if serial)
  Teuchos::RCP<Epetra_MultiVector> SerialB_ ;
  Teuchos::RCP<Epetra_MultiVector> SerialX_ ;
  //! Serial versions of the LHS and RHS (if necessary)
  Teuchos::RCP<Epetra_MultiVector> SerialXextract_;
  Teuchos::RCP<Epetra_MultiVector> SerialBextract_;

  //! If \c true, the transpose of A is used.
  bool UseTranspose_;
  //! Pointer to the linear system problem.
  const Epetra_LinearProblem * Problem_;

  //! Only used for RowMatrices to extract copies.
  std::vector<int> ColIndicesV_;
  //! Only used for RowMatrices to extract copies.
  std::vector<double> RowValuesV_;
  //! Importer to process 0.
  Teuchos::RCP<Epetra_Import> ImportToSerial_;
  Teuchos::RCP<Epetra_Import> ImportRangeToSerial_;
  Teuchos::RCP<Epetra_Import> ImportDomainToSerial_;
  //! Quick access ids for the individual timings
  int MtxRedistTime_, MtxConvTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_, OverheadTime_;

};  // class Amesos_Klu

#endif /* AMESOS_KLU_H */
