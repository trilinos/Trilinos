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

#ifndef AMESOS_MUMPS_H
#define AMESOS_MUMPS_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

class Epetra_Import;
class Epetra_RowMatrix;
class Epetra_MultiVector;
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;
class Amesos_EpetraInterface;

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Teuchos_RCP.hpp"
#include <map>
using namespace Teuchos;

//! Amesos_Mumps:  An object-oriented wrapper for the double precision version of MUMPS.
/*!  Amesos_Mumps is an interface to the the double precision version of 
  the sparse parallel direct
  solver MUMPS.  Given an Epetra_RowMatrix A, and two
  Epetra_MultiVectors X and B, the solution with Amesos_Mumps reads as
  follows:

  -# Epetra_LinearProblem Problem; Amesos_BaseSolver *
  Solver; Amesos Amesos_Factory;
  -# Solver = Amesos_Factory.Create("Amesos_Mumps", Problem);
  -# if( Solver == 0 ) std::cerr << "library not available" << std::endl;
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
  
  \warning This interface is compatible with MUMPS 4.5.4.

  \date Last modified 26-Jan-06

  \author Marzio Sala, ETHZ.
  
*/

extern "C" {
#include "dmumps_c.h"
}

class Amesos_Mumps: public Amesos_BaseSolver,
                    private Amesos_Time, 
                    private Amesos_NoCopiable, 
                    private Amesos_Utils,  
                    private Amesos_Control,  
                    private Amesos_Status { 

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

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //! Destroys all data associated with \sl this object.
  void Destroy();
  
  int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    if (UseTranspose_in)
      return (-1);
    return (0);
  };
  
  bool UseTranspose() const {return(UseTranspose_);};

  int SetParameters( Teuchos::ParameterList &ParameterList );
  
  //@}

  //@{ \name Solver status functions

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //@}

  //@{ \name Print functions
  
  //! Prints timing information.
  /*! In the destruction phase, prints out detailed information about
    the various phases: symbolic and numeric factorization, solution,
    gather/scatter for vectors and matrices.
   */
  void PrintTiming() const;
  
  //! Prints information about the factorization and solution phases.
  /*! In the destruction phase, prints out some information furnished by
    MUMPS, like the amount of required memory, the MFLOPS.
   */
  void PrintStatus() const;

  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }

  //@}

  //@{ \name MUMPS' specify functions

  
#if 0
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
#endif
  
  //! Set prescaling.
  /*! Use double precision vectors of size N (global dimension of the matrix) as
      scaling for columns and rows. \c ColSca and \c RowSca must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(8) = -1.

      Both input vectors are \c float with --enable-amesos-smumps, \c double otherwise.
      
  */
  int SetPrecscaling(double * ColSca, double * RowSca )
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
  int SetRowScaling(double * RowSca )
  {
    RowSca_ = RowSca;
    return 0;
  }

  //! Set column scaling
  /*! Use double precision vectors of size N (global dimension of the matrix) for column
      scaling.

      \param ColSca (In) - \c float pointer with --enable-amesos-smumps, \c double pointer otherwise.
  */
  int SetColScaling(double * ColSca )
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

  //! Gets the pointer to the RINFO array (defined on all processes).
  /*! Gets the pointer to the internally stored RINFO array, of type \c
    float if option \c --enable-amesos-smumps is enabled, \c double
    otherwise.
   */
  double * GetRINFO() ;

  //! Gets the pointer to the INFO array (defined on all processes).
  /*! Gets the pointer to the internally stored INFO array, of type \c int.
   */
  int * GetINFO() ;

  //! Gets the pointer to the RINFOG array (defined on host only).
  /*! Gets the pointer to the internally stored RINFOG array (defined on
    the host process only), of type \c float if option \c
    --enable-amesos-smumps is enabled, \c double otherwise.
   */
  double * GetRINFOG() ;

  //! Get the pointer to the INFOG array (defined on host only).
  /*! Gets the pointer to the internally stored INFOG (defined on the
    host process only) array, of type \c int.
   */
  int * GetINFOG() ;

  //! Set ICNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  void SetICNTL(int pos, int value);

  //! Set CNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  void SetCNTL(int pos, double value);

  //@}
  
  bool MatrixShapeOK() const
  {
  bool OK = true;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}

  
  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Gets a pointer to the Epetra_LinearProblem.
  const Epetra_LinearProblem * GetProblem() const { return(Problem_); };

protected:
  
  //! Returns a reference to the linear system matrix.
  Epetra_RowMatrix& Matrix();

  const Epetra_RowMatrix& Matrix() const;

  //! Returns a reference to the map for redistributed matrix.
  Epetra_Map& RedistrMap();

  //! Returns a reference for the redistributed importer.
  Epetra_Import& RedistrImporter();
  
  //! Returns a reference to the redistributed matrix, imports it is \c ImportMatrix is true.
  Epetra_RowMatrix& RedistrMatrix(const bool ImportMatrix = false);

  //! Returns a reference to the map with all elements on process 0.
  Epetra_Map& SerialMap();

  //! Returns a reference to the importer for SerialMap().
  Epetra_Import& SerialImporter();

  //! Converts to MUMPS format (COO format).
  int ConvertToTriplet(const bool OnlyValues);     

  //! Checks for MUMPS error, prints them if any. See MUMPS' manual.
  int CheckError();

  //! Check input parameters.
  void CheckParameters();
  
  void SetICNTLandCNTL();

  //! \c true if matrix has already been converted to COO format
  bool IsConvertToTripletOK_;
  //! \c true if the Schur complement has been computed (need to free memory)
  bool IsComputeSchurComplementOK_;


  bool NoDestroy_ ;  // Set true to prevent memory freeing
  
  //! row indices of nonzero elements
  std::vector <int> Row;
  //! column indices of nonzero elements
  std::vector<int> Col;
  //! values of nonzero elements
  std::vector<double> Val;

  //! Maximum number of processors in the MUMPS' communicator
  int MaxProcs_;
  
  //! If \c true, solve the problem with AT.
  bool UseTranspose_;

  //! Quick access pointers to the internal timers
  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_;  

  //! Row and column scaling
  double * RowSca_, * ColSca_; 

  //! PermIn for MUMPS.
  int * PermIn_;

  //! Number of rows in the Schur complement (if required)
  int NumSchurComplementRows_;
  //! Rows for the Schur complement (if required)
  int * SchurComplementRows_;

  //! Pointer to the Schur complement, as CrsMatrix.
  RCP<Epetra_CrsMatrix> CrsSchurComplement_; 
  //! Pointer to the Schur complement,as DenseMatrix.
  RCP<Epetra_SerialDenseMatrix> DenseSchurComplement_;

#ifdef HAVE_MPI
  //! MPI communicator used by MUMPS
  MPI_Comm MUMPSComm_;
#endif

  //! Pointer to the linear problem to be solved.
  const Epetra_LinearProblem* Problem_;

  //! Redistributed matrix.
  RCP<Epetra_Map> RedistrMap_;
  //! Redistributed importer (from Matrix().RowMatrixRowMap() to RedistrMatrix().RowMatrixRowMap()).
  RCP<Epetra_Import> RedistrImporter_;
  //! Redistributed matrix (only if MaxProcs_ > 1).
  RCP<Epetra_CrsMatrix> RedistrMatrix_;
  //! Map with all elements on process 0 (for solution and rhs).
  RCP<Epetra_Map> SerialMap_;
  //! Importer from Matrix.OperatorDomainMap() to SerialMap_.
  RCP<Epetra_Import> SerialImporter_;
  
  DMUMPS_STRUC_C MDS;

  std::map<int, int> ICNTL;
  std::map<int, double> CNTL;
};  // class Amesos_Mumps

#endif /* AMESOS_MUMPS_H */
