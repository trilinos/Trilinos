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

#ifndef AMESOS_MUMPS_H
#define AMESOS_MUMPS_H

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

#ifndef HAVE_AMESOS_SMUMPS
#define AMESOS_TYPE double
#else
#define AMESOS_TYPE float
#endif

extern "C" {
#ifndef HAVE_AMESOS_SMUMPS
#include "dmumps_c.h"
#else
#include "smumps_c.h"
#endif
}

//! Amesos_Mumps:  An object-oriented wrapper for CERFACS' MUMPS.
/*!  Amesos_Mumps is an interface to the CERFACS' sparse parallel direct
  solver MUMPS.  Given an Epetra_RowMatrix A, and two
  Epetra_MultiVectors X and B, the solution with Amesos_Mumps reads as
  follows:

  -# Epetra_LinearProblem Problem; Amesos_BaseSolver *
  Solver; Amesos Amesos_Factory;
  -# Solver = Amesos_Factory.Create("Amesos_Mumps", Problem);
  -# if( Solver == 0 ) cerr << "library not available" << endl;
  -# Problem.SetMatrix(&A);
  -# Solver->SymbolicFactorization();
  -# Solver->NumericFactorization();
  -# Problem.SetLHS(&X);
  -# Problem.SetLHS(&B);
  -# Solver->Solve();
  
  A number of parameters is available to tune the performances of
  MUMPS. We refer to the Amesos Reference Guide for a detailed overview
  of these parameters. Here, we just recall that it is possible to solve
  the linear system on a subset of the processes contained in the Comm
  object of the Epetra_LinearProblem.

  Amesos_Mumps accepts any Epetra_RowMatrix derived class. However,
  special functions are available for Epetra_CrsMatrix and
  Epetra_VbrMatrix objects.

  The single-precision version of MUMPS can be used by enabling the
  option \c --enable-amesos-smumps.  Note that this option overwrites
  --enable-amesos-mumps.  The choice between single-precision and
  double-precision must be done at configuration (and compilation)
  time.

  As Amesos is based on Epetra, and Epetra is only double-precision, we
  still require an Epetra_LinearProblem composed by a double-precision
  matrix, and two double-precision vectors. The solution vector is
  casted to \c double after solution. Single precision may be of
  interest if Amesos is used with ML, to solve the coarse problem (for
  which single-precision can be enough in term of numerical error, and
  usually save memory and CPU time).

  Amesos_Mumps is based on Amesos_EpetraBaseSolver, that is derived from
  Amesos_BaseSolver. The main redistribution utilities, as well as a
  getrow function, is obtained by EpetraBaseSolver.
  
  \warning This interface has been developed with MUMPS 4.3.1.

  \author Marzio Sala, 9214
  
*/
class Amesos_Mumps : public Amesos_EpetraBaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Mumps Constructor.
  /*! Creates an Amesos_Mumps instance, using an Epetra_LinearProblem,
    
  */
  Amesos_Mumps(const Epetra_LinearProblem& LinearProblem);

  //! Amesos_Mumps Destructor.
  /*! Deletes an Amesos_Mumps object.  
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

      At this point, the numerical values of A are not required.
      
    \return Integer error code, set to 0 if successful.
  */
  int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*! Performs the numeric factorization (and symbolic factorization
      if necessary) on the matrix A. It is supposed that the structure
      of the matrix A has not changed since the previous call the
      SymbolicFactorization() (if any).

      At this point, LHS and RHS of the Epetra_LinearProblem are not required.
     \return Integer error code, set to 0 if successful.
  */
  int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! Solve the linear system for all vectors contained in the
      Epetra_MultiVector RHS(), and store the solution in LHS().

      By default, Solve() will solve the problem with A, and not with
      A<SUP>T</SUP>.  Users can solve the problem with A<SUP>T</SUP> by
      creating a Teuchos::ParameterList (say, AmesosList), set
      AmesosList.set("UseTranspose",true), and calling
      SetParameters(AmesosList).
      
     \return Integer error code, set to 0 if successful.
  */
  int Solve();

  //! Destroys all data associated with \sl this object.
  void Destroy();
  
  //  char * Label() const {return(Epetra_Object::Label());};

  //! If set true, solve the problem with A<SUP>T</SUP>
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Sets parameters for \sl this object from input list
  /*! Sets all the parameters for \sl this object, retriving them form
    the input Teuchos::ParameterList. This call can modify the input
    list; default values (not found in the list) are added. For a
    detailed overview of the available parameters, please refer to the
    Amesos Reference Guide (in \c
    Trilinos/packages/amesos/doc/AmesosReferenceGuide/AmesosReferenceGuide.pdf).
  */
  int SetParameters(Teuchos::ParameterList &ParameterList );
  
  //@}

  //@{ \name Print functions
  
  //! Prints timing information.
  /*! In the destruction phase, prints out detailed information about
    the various phases: symbolic and numeric factorization, solution,
    gather/scatter for vectors and matrices.
   */
  void PrintTiming();
  
  //! Prints information about the factorization and solution phases.
  /*! In the destruction phase, prints out some information furnished by
    MUMPS, like the amount of required memory, the MFLOPS.
   */
  void PrintStatus();

  //@}

  //@{ \name MUMPS' specify functions

  
  //! Returns the Schur complement matrix as an Epetra_CrsMatrix.
  /*! Returns the (dense) Schur complement matrix as an Epetra_CrsMatrix. This
      matrix is defined on all the processes in the Epetra Communicator. However,
      it has rows on the host process only.
      If \In flag : if \c true, MUMPS will compute the Schur complement matrix,
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
  
  
  //! Set prescaling.
  /*! Use double precision vectors of size N (global dimension of the matrix) as
      scaling for columns and rows. \c ColSca and \c RowSca must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(8) = -1.

      Both input vectors are \c float with --enable-amesos-smumps, \c double otherwise.
      
  */
  int SetPrecscaling(AMESOS_TYPE * ColSca, AMESOS_TYPE * RowSca )
  {
    ColSca_ = ColSca;
    RowSca_ = RowSca;
    return 0;
  }

  //! Set row scaling
  /*! Use double precision vectors of size N (global dimension of the matrix) for row
      scaling.

      \param RowSca (In) -\c float pointer with --enable-amesos-smumps, \c double pointer otherwise.
  */
  int SetRowScaling(AMESOS_TYPE * RowSca )
  {
    RowSca_ = RowSca;
    return 0;
  }

  //! Set column scaling
  /*! Use double precision vectors of size N (global dimension of the matrix) for column
      scaling.

      \param ColSca (In) - \c float pointer with --enable-amesos-smumps, \c double pointer otherwise.
  */
  int SetColScaling(AMESOS_TYPE * ColSca )
  {
    ColSca_ = ColSca;
    return 0;
  }

  //! Sets ordering.
  /*! Use integer vectors of size N (global dimension of the matrix) as
      given ordering. \c PermIn must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(7) = 1.
  */
  int SetOrdering(int * PermIn)
  {
    PermIn_ = PermIn;
    return 0;
  }

  //! Sets the Maxis value (see MUMPS' manual)
  int SetMaxis(int Maxis)
  {
    Maxis_ = Maxis;
    return 0;
  }

  //! Sets the Maxs value (see MUMPS' manual)
  int SetMaxs( int Maxs) 
  {
    Maxs_ = Maxs;
    return 0;
  }

  //! Gets the pointer to the RINFO array (defined on all processes).
  /*! Gets the pointer to the internally stored RINFO array, of type \c
    float if option \c --enable-amesos-smumps is enabled, \c double
    otherwise.
   */
  AMESOS_TYPE * GetRINFO() 
  {
    return (MDS.rinfo);
  }

  //! Gets the pointer to the INFO array (defined on all processes).
  /*! Gets the pointer to the internally stored INFO array, of type \c int.
   */
  int * GetINFO() 
  {
    return (MDS.info);
  }

  //! Gets the pointer to the RINFOG array (defined on host only).
  /*! Gets the pointer to the internally stored RINFOG array (defined on
    the host process only), of type \c float if option \c
    --enable-amesos-smumps is enabled, \c double otherwise.
   */
  AMESOS_TYPE * GetRINFOG()
  {
    return (MDS.rinfog);
  }

  //! Get the pointer to the INFOG array (defined on host only).
  /*! Gets the pointer to the internally stored INFOG (defined on the
    host process only) array, of type \c int.
   */
  int * GetINFOG()
  {
    return (MDS.infog);
  }

  //! Copies the input array (of size 40) into the internally stored ICNTL array.
  int SetICNTL(int * ictnl);

  //! Set ICNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetICNTL(int pos, int value);

  //! Copy the input array (of size 5) into the internally stored CNTL array.
  int SetCNTL(double * ctnl);

  //! Set CNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetCNTL(int pos, double value);

  void SetUseMpiCommSelf() {
    UseMpiCommSelf_ = true;
  }

  //@}
  
protected:
  
  //! Converts to MUMPS format (COO format).
  int ConvertToTriplet();     

  //! Converts to MUMPS format (for the values only).
  int ConvertToTripletValues();

  //! Performs the symbolic factorization.      
  int PerformSymbolicFactorization(); 

  //! Performs the numeric factorization
  int PerformNumericFactorization(); 

  //! Checks for MUMPS error, prints them if any. See MUMPS' manual.
  void CheckError();

  //! Check input parameters.
  void CheckParameters();
  
  void SetICNTLandCNTL();

  //! Redistributed input matrix over the specified number of processes.
  void RedistributeMatrix(const int NumProcs);

  //! Redistributed matrix values over the specified number of processes.
  void RedistributeMatrixValues(const int NumProcs);
  
  //! \c true if SymbolicFactorization has been done
  bool SymbolicFactorizationOK_;
  //! \c true if NumericFactorization has been done
  bool NumericFactorizationOK_;
  //! \c true if matrix has already been converted to COO format
  bool IsConvertToTripletOK_;
  //! \c true if the Schur complement has been computed (need to free memory)
  bool IsComputeSchurComplementOK_;
  //! \c true if only local entries must be considered
  bool UseMpiCommSelf_;

#ifndef HAVE_AMESOS_SMUMPS  
  //! Mumps data structure for double-precision
  DMUMPS_STRUC_C MDS;
#else
  //! Mumps data structure for single-precision
  SMUMPS_STRUC_C MDS;
#endif
  
  //! row indices of nonzero elements
  Epetra_IntSerialDenseVector* Row;
  //! column indices of nonzero elements
  Epetra_IntSerialDenseVector* Col;
  //! values of nonzero elements
  Epetra_SerialDenseVector* Val;

#ifdef HAVE_AMESOS_SMUMPS
  //! single-precision values of nonzero elements
  float * SVal;
  //! single-precision solution vector (on host only)
  float * SVector;
#endif
  
  //!  Number of non-zero entries in Problem_->GetOperator()
  int numentries_;         
  //!  Number of rows and columns in the Problem_->GetOperator()
  int NumGlobalElements_;  

  bool  KeepMatrixDistributed_;          /*!< this governs the ICNTL(18) parameter.
                                            If false, then matrix is redistributed
                                            to proc 0 before converting it to
                                            triplet format. Then, MUMPS will take care
                                            of reditribution. If true, the input
                                            distributed matrix is passed to MUMPS. */

  //! Maximum number of processors in the MUMPS' communicator
  int MaxProcs_;
  //! Maximum number of processors that contains at least one element of A
  int MaxProcsInputMatrix_;
  
  //! Maps to redistribute from Matrix' Map to MaxProcs_ Map
  const Epetra_Map * Map_;

  //! actual number of global nonzeros in the matrix
  int NumMUMPSNonzeros_;
  //! actual number of local nonzeros in the matrix
  int NumMyMUMPSNonzeros_;
  
  //! If \c true, solve the problem with AT.
  bool UseTranspose_;
  //! If \c true, add a the value AddToDiag_ on the diagonal
  bool AddDiagElement_;
  
  //! Value to add to the diagonal if specified 
  double AddToDiag_;
  
  //! If \c true, print timing in the destruction phase
  bool PrintTiming_;
  //! If \c true, print status in the destruction phase
  bool PrintStatus_;
  //! If \c true, compute the norms of solution and RHS.
  bool ComputeVectorNorms_;
  //! If \c true, compute the true residual after solution
  bool ComputeTrueResidual_;

  //! Set the matrix property.
  /*! Matrix property can be 
    - 0 : general unsymmetric matrix;
    - 1 : SPD;
    - 2 : general symmetric matrix.
    */
  int MatrixProperty_;        
  //! Discard all elements whose absolute value is below this value
  double Threshold_;
  
  int icntl_[40];        
  AMESOS_TYPE cntl_[5];  
 
  //! Row and column scaling
  AMESOS_TYPE * RowSca_, * ColSca_; 

  //! PermIn for MUMPS.
  int * PermIn_;
  //! MAXIS for MUMPS.
  int Maxis_;
  // MAXS for MUMPS.
  int Maxs_;  

  //! Number of rows in the Schur complement (if required)
  int NumSchurComplementRows_;
  //! Rows for the Schur complement (if required)
  int * SchurComplementRows_;

  //! Pointer to the Schur complement, as CrsMatrix.
  Epetra_CrsMatrix * CrsSchurComplement_; 
  //! Pointer to the Schur complement,as DenseMatrix.
  Epetra_SerialDenseMatrix * DenseSchurComplement_;

  //! ID of calling process.
  int MyPID_;
  //! Number of processes in computation.
  int NumProcs_;
  //! Output level.
  int verbose_;
  
  EpetraExt_Redistor * Redistor_;
  
  Epetra_RowMatrix * OldMatrix_;

  Epetra_MultiVector * TargetVector_;

  //! time to convert to MUMPS format
  double ConTime_;
  //! time for symbolic factorization
  double SymTime_;
  //! time for numeric factorization
  double NumTime_;
  //! time for solution
  double SolTime_;
  //! time to redistribute vectors
  double VecTime_;
  //! time to redistribute matrix
  double MatTime_;
  
  //! Number of symbolic factorization phases
  int NumSymbolicFact_;
  //! Number of symbolic numeric phases
  int NumNumericFact_;
  //! Number of symbolic solution phases
  int NumSolve_;

  //! Used to track timing
  Epetra_Time * Time_;
  
#ifdef EPETRA_MPI
  //! MPI communicator used by MUMPS
  MPI_Comm MUMPSComm_;
#endif
  
};  // class Amesos_Mumps

#endif /* AMESOS_MUMPS_H */
