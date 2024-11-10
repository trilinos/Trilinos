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
 * \file Amesos_Paraklete.h
 *
 * \class Amesos_Paraklete
 *
 * \brief Interface to PARAKLETE internal solver.Interface to PARAKLETE internal solver.
 *
 * \date Last updated on 24-May-05.
 */

#ifndef AMESOS_PARAKLETE_H
#define AMESOS_PARAKLETE_H

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
#ifdef HAVE_AMESOS_EPETRAEXT
#include "EpetraExt_Transpose_RowMatrix.h"
#endif


//! Amesos_Paraklete:  A serial, unblocked code ideal for getting started and for very sparse matrices, such as circuit matrces.

/*! 

Class Amesos_Paraklete is an object-oriented wrapper for PARAKLETE. PARAKLETE, whose sources
are distributed
within Amesos, is a serial solver for sparse matrices. PARAKLETE will solve a 
linear system of equations: \f$A X = B\f$, where
<TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
Epetra_MultiVector objects.

Amesos_Paraklete computes \f$A^T X = B\f$ 
more efficiently than \f$>A X = B\f$.  The
latter requires a matrix transpose -- which costs both time and space.

Paraklete is Tim Davis' parallel version of KLU a low overhead non-blocked code which 
solves very sparse matrices fast.

*/

// Amesos_Paraklete_Pimpl contains a pointer to structures defined in 
// paraklete.h.  This prevents Amesos_Paraklete.h 
// from having to include paraklete.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Paraklete_Pimpl ; 
class Amesos_StandardIndex ; 
#endif

class Amesos_Paraklete: public Amesos_BaseSolver,  
                  private Amesos_Time, 
                  private Amesos_NoCopiable, 
                  private Amesos_Utils, 
                  private Amesos_Control, 
                  private Amesos_Status { 

public: 

  //@{ \name Constructors and Destructors
  //! Amesos_Paraklete Constructor.
  /*! Creates an Amesos_Paraklete instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Paraklete(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Paraklete Destructor.
  ~Amesos_Paraklete(void);
  
  //@}
  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name 

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if PARAKLETE can handle this matrix shape 
  /*! Returns true if the matrix shape is one that PARAKLETE can
    handle. PARAKLETE only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose()
  /*! 
    If SetUseTranspose() is set to true, 
    \f$A^T X = B\f$ is computed.
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
 
  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }
 
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
    ConvertToParakleteCRS - Convert matrix to form expected by Paraklete: Ai, Ap, Aval
    Preconditions:
      numentries_, RowMatrixA_, ImportToSerial_, StdIndexMatrix_, Reindex_ 
    Postconditions:
      SerialCrsMatrixA_
  */
  int ConvertToParakleteCRS(bool firsttime);     

  /*
    PerformSymbolicFactorization - Call Paraklete to perform symbolic factorization
    Preconditions:
      UseDataInPlace_ must be set to 1 if the input matrix is entirely stored on process 0
      Ap, Ai and Aval point to a compressed row storage version of the input matrix A.
    Postconditions:
      Symbolic points to an PARAKLETE internal opaque object containing the
        symbolic factorization and accompanying information.  
      SymbolicFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
      
  int PerformSymbolicFactorization(); 

  /*
    PerformNumericFactorization - Call Paraklete to perform numeric factorization
    Preconditions:
      UseDataInPlace_ must be set 
      Ap, Ai and Aval point to a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an PARAKLETE internal opaque object containing the
        numeric factorization and accompanying information.  
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

  // @}

  bool IamInGroup_;  // True if this process is involved in the computation.  Set by SymbolicFactorization
  
  int SerialXlda_ ;
  //
  //  PrivateParakleteData_ contains pointers to data needed by paraklete whose
  //  data structures are defined by paraklete.h
  //
  Teuchos::RCP<Amesos_Paraklete_Pimpl> PrivateParakleteData_; 
  Teuchos::RCP<Amesos_StandardIndex> StdIndex_; 
  Teuchos::RCP<Amesos_StandardIndex> StdIndexRange_; 
  Teuchos::RCP<Amesos_StandardIndex> StdIndexDomain_; 
  MPI_Comm ParakleteComm_;

  //! Ap, Ai, Aval form the compressed row storage used by Paraklete
  //! Ai and Aval can point directly into a matrix if it is StorageOptimized(), hence
  //! they may either be in vector form or may be a pointer into Epetra_CrsMatrix 
  //! internals.  Ap must always be constructed.  
  std::vector <long> Ap;
  std::vector <long> Ai;
  std::vector <double> VecAval;
  double* Aval;

  //! 1 if Problem_->GetOperator() is stored entirely on process 0
  int UseDataInPlace_;
  //! Number of non-zero entries in Problem_->GetOperator()
  int numentries_;
  //! Number of rows and columns in the Problem_->GetOperator()
  int NumGlobalElements_;

  //! Operator converted to a RowMatrix
  Epetra_RowMatrix* RowMatrixA_;
  //! Operator converted to a CrsMatrix
  Epetra_CrsMatrix* CrsMatrixA_;
  //
  //  transposer_ transposes a CrsMatrix 
  //  Created in CreateLocalMatrixAndExporters
  //  Used in ExportToSerial()
  //
#ifdef HAVE_AMESOS_EPETRAEXT
  Teuchos::RCP<EpetraExt::RowMatrix_Transpose> transposer_;
#endif
  //! Points to a Serial Map (unused if UseDataInPlace_ == 1 )
  Teuchos::RCP<Epetra_Map> SerialMap_;
  //! Points to a Serial Copy of A (unused if UseDataInPlace_==1)
  Teuchos::RCP<Epetra_CrsMatrix> SerialCrsMatrixA_;
  //! Points to a Contiguous Copy of A 
  Epetra_RowMatrix* StdIndexMatrix_ ; 
  Epetra_MultiVector* StdIndexDomainVector_ ; 
  Epetra_MultiVector* StdIndexRangeVector_ ; 
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
  Epetra_MultiVector* SerialB_ ;
  Epetra_MultiVector* SerialX_ ;
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
 
  //! Quick access pointers to internal timing information
  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_, OverheadTime_;
 
};  // class Amesos_Paraklete  

#endif /* AMESOS_PARAKLETE_H */
