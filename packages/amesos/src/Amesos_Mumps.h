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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef _AMESOS_MUMPS_H_
#define _AMESOS_MUMPS_H_

class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;
#include "Epetra_SerialDenseVector.h"
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;
class Amesos_EpetraInterface;
class EpetraExt_Redistor;
#include "Epetra_Time.h"

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Amesos_EpetraBaseSolver.h"

extern "C" {
#include "dmumps_c.h"
}

//! Amesos_Mumps:  An object-oriented wrapper for Mumps.
/*!  Amesos_Mumps will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Mumps solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

  Mumps execution can be tuned through a variety of parameters.
  Amesos_Mumps.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.

  MUMPS will not perform column permutation on a matrix provided
  in distributed form.  Amesos_Mumps will match this, allowing
  column permutation only if the matrix is provided in serial form.
  This is unfortunate because it is an exception to the general rule
  that the capability (given adequate memory) of any class
  implementing the Amesos_BaseSolver base class does not depend on
  the distribution of the input matrix.  However, neither of the
  other options are attractive.  Coalescing the matrix to a single
  process independent of whether column permutation is requested
  unnecessarily limits the size problem that can be solved.
  Coalescing the matrix to a single process only when column
  permutation is requested would cause some problems to run out of memory
  when column permutation is requested.


  \Note This class should be used with MUMPS 4.3 or 4.3.1 (never tested
  with older versions of MUMPS, and developed with 4.3.1).

  \author Marzio Sala, 9214
  
*/
class Amesos_Mumps : public Amesos_EpetraBaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Mumps Constructor.
  /*! Creates an Amesos_Mumps instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Mumps(const Epetra_LinearProblem& LinearProblem);

  //! Amesos_Mumps Destructor.
  /*! Completely deletes an Amesos_Mumps object.  
  */
  ~Amesos_Mumps(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      </ul>

      postconditions:<ul>
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      <li>MDS will be modified to reflect the symbolic factorization 
      which has been performed.
      </ul>

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      </ul>

      postconditions:<ul>
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      <li>MDS will be modified to reflect the numeric factorization 
      which has been performed.
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  void Destroy();
  
  //@}

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Returns true if MUMPS can handle this matrix shape 
  /*! Returns true if the matrix shape is one that MUMPS can
    handle. MUMPS only works with square matrices.  
  */

  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};


  //! Returns the Schur complement matrix as an Epetra_CrsMatrix.
  /*! Returns the (dense) Schur complement matrix as an Epetra_CrsMatrix. This
      matrix is defined on all the processes in the Epetra Communicator. However,
      it has rows on the host process only.
      If \in flag : if \c true, MUMPS will compute the Schur complement matrix,
      with respect to the (global) rows defined in the integer array
      \c SchurComplementRows, of size \c NumSchurComplementRows.
      Those two arrays are defined on the host only.
  */
  int ComputeSchurComplement(bool flag,
			     int NumSchurComplementRows, int * SchurComplementRows);

  //! Returns the Schur complement in an Epetra_CrsMatrix on host only.
  /*! Returns the Schur complement in an Epetra_CrsMatrix on host only. Note that
      no checks are performed to see whether this action is legal or not (that is,
      if the call comes after the solver has been invocated).
      Epetra_CrsMatrix must be freed by the user!
  */
  Epetra_CrsMatrix * GetCrsSchurComplement();

  //! Returns the Schur complement as a SerialDenseMatrix (on host only).
  /*! Returns the Schur complement in an Epetra_SerialDenseMatrix on host only. Note that
      no checks are performed to see whether this action is legal or not (that is,
      if the call comes after the solver has been invocated).
      Epetra_SerialDenseMatrix must be freed by the user!
  */
  Epetra_SerialDenseMatrix * GetDenseSchurComplement();
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  int SetParameters(Teuchos::ParameterList &ParameterList );
  
  //@}

  //! Set prescaling.
  /*! Use double precision vectors of size N (global dimension of the matrix) as
      scaling for columns and rows. \c ColSca and \c RowSca must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(8) = -1.
  */
  int SetPrecscaling(double * ColSca, double * RowSca )
  {
    ColSca_ = ColSca;
    RowSca_ = RowSca;
    return 0;
  }

  int SetRowScaling(double * RowSca )
  {
    RowSca_ = RowSca;
    return 0;
  }

  int SetColScaling(double * ColSca )
  {
    ColSca_ = ColSca;
    return 0;
  }

  //! Set ordering.
  /*! Use integer vectors of size N (global dimension of the matrix) as
      given ordering. \c PermIn must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(7) = 1.
  */
  int SetOrdering(int * PermIn)
  {
    PermIn_ = PermIn;
    return 0;
  }

  int SetMaxis(int Maxis)
  {
    Maxis_ = Maxis;
    return 0;
  }

  int SetMaxs( int Maxs) 
  {
    Maxs_ = Maxs;
    return 0;
  }

  //! Get the pointer to the RINFO array (defined on all processes).
  double * GetRINFO() 
  {
    return (MDS.rinfo);
  }

  //! Get the pointer to the INFO array (defined on all processes).
  int * GetINFO() 
  {
    return (MDS.info);
  }

  //! Get the pointer to the RINFOG array (defined on host only).
  double * GetRINFOG()
  {
    return (MDS.rinfog);
  }

  //! Get the pointer to the INFOG array (defined on host only).
  int * GetINFOG()
  {
    return (MDS.infog);
  }

  //! Copy the input array into the internally stored ICNTL array.
  int SetICNTL(int * ictnl);

  //! Set ICNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetICNTL(int pos, int value);

  //! Copy the input array into the internally stored CNTL array.
  int SetCNTL(double * ctnl);

  //! Set CNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetCNTL(int pos, double value);

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();

  void SetUseMpiCommSelf() {
    UseMpiCommSelf_ = true;
  }

protected:
  
  /*
    ConvertToTriplet - Convert matrix to form expected by Mumps: Row, Col, Val
    Preconditions:
      numentries_, NumGloalElements_ and SerialMatrix_ must be set.
    Postconditions:
      Row, Col, Val are resized and populated with values that reflect
      the input matrix A.

  */
  int ConvertToTriplet();     
  int ConvertToTripletValues();
  /*
    PerformSymbolicFactorization - Call Mumps to perform symbolic factorization
    Preconditions:
      Row, Col and Val must be set

    Postconditions:
      MDS will be modified to reflect the symbolic factorization 
      which has been performed.
      SymbolicFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
      
  int PerformSymbolicFactorization(); 

  /*
    PerformNumericFactorization - Call Mumps to perform numeric factorization
    Preconditions:
      IsLocal must be set to 1 if the input matrix is entirely stored on process 0
      Row, Col and Val must be set
    Postconditions:
      MDS will be modified to reflect the numeric factorization 
      which has been performed.
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

  void CheckError();

  void CheckParameters();
  
  void SetICNTLandCNTL();

  void RedistributeMatrix(const int NumProcs);

  void RedistributeMatrixValues(const int NumProcs);
  
  bool SymbolicFactorizationOK_;   // True if SymbolicFactorization has been done
  bool NumericFactorizationOK_;    // True if NumericFactorization has been done
  bool IsConvertToTripletOK_;
  bool IsComputeSchurComplementOK_;
  bool UseMpiCommSelf_;
  
  DMUMPS_STRUC_C MDS ;             // Mumps data structure 

  //
  //  Row, Col, Val form the triplet representation used by Mumps
  //
  Epetra_IntSerialDenseVector * Row; // store COO format in epetra vectors
  Epetra_IntSerialDenseVector * Col;
  Epetra_SerialDenseVector    * Val;

  int MyPID;               //  Process number (i.e. Comm().MyPID() 
  
  int numentries_;         //  Number of non-zero entries in Problem_->GetOperator()
  int NumGlobalElements_;  //  Number of rows and columns in the Problem_->GetOperator()

  bool  KeepMatrixDistributed_;          // this governs the ICNTL(18) parameter.
                                         // If false, then matrix is redistributed
                                         // to proc 0 before converting it to
                                         // triplet format. Then, MUMPS will take care
                                         // of reditribution. If true, the input
                                         // distributed matrix is passed to MUMPS.

  int MaxProcs_;
  int MaxProcsInputMatrix_;
  
  const Epetra_Map * Map_;

  int NumMUMPSNonzeros_;                  // actual number of nonzeros in the matrix
  int NumMyMUMPSNonzeros_;                // actual number of nonzeros in the matrix
  
  bool UseTranspose_;
  bool AddDiagElement_;
  
  double AddToDiag_;
  
  bool PrintTiming_;
  bool PrintStatus_;
  bool ComputeVectorNorms_;
  bool ComputeTrueResidual_;
  
  double Threshold_;
  
  int icntl_[40];                         // to allow users overwrite default settings
  double cntl_[5];                        // as specified by Amesos
  double * RowSca_, * ColSca_;
  int * PermIn_;
  int Maxis_, Maxs_;  

  int NumSchurComplementRows_;            // Schur complement section
  int * SchurComplementRows_;

  Epetra_CrsMatrix * CrsSchurComplement_;
  Epetra_SerialDenseMatrix * DenseSchurComplement_;

  int verbose_;
  int debug_;
  
  EpetraExt_Redistor * Redistor_;
  
  Epetra_RowMatrix * OldMatrix_;

  Epetra_MultiVector * TargetVector_;

  // some timing internal to MUMPS
  double ConTime_;                        // time to convert to MUMPS format
  double SymTime_;                        // time for symbolic factorization
  double NumTime_;                        // time for numeric factorization
  double SolTime_;                        // time for solution
  double VecTime_;                        // time to redistribute vectors
  double MatTime_;                        // time to redistribute matrix
  
  int NumSymbolicFact_;
  int NumNumericFact_;
  int NumSolve_;
  

  Epetra_Time Time;
  
#ifdef EPETRA_MPI
  MPI_Comm MUMPSComm_;
#endif
  
};  // End of  class Amesos_Mumps

#endif /* _AMESOS_MUMPS_H_ */
