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

#ifndef HAVE_AMESOS_SMUMPS
#define AMESOS_TYPE double
#else
#define AMESOS_TYPE float
#endif


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

// Amesos_Mumps_Pimpl contains a pointer to the structures defined in 
// dmumps.h and smumps.h.  This prevents Amesos_Mumps.h 
// from having to include dmumps.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Mumps_Pimpl ; 
#endif

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
  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};
  
  bool UseTranspose() const {return(UseTranspose_);};

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
  AMESOS_TYPE * GetRINFO() ;

  //! Gets the pointer to the INFO array (defined on all processes).
  /*! Gets the pointer to the internally stored INFO array, of type \c int.
   */
  int * GetINFO() ;

  //! Gets the pointer to the RINFOG array (defined on host only).
  /*! Gets the pointer to the internally stored RINFOG array (defined on
    the host process only), of type \c float if option \c
    --enable-amesos-smumps is enabled, \c double otherwise.
   */
  AMESOS_TYPE * GetRINFOG() ;

  //! Get the pointer to the INFOG array (defined on host only).
  /*! Gets the pointer to the internally stored INFOG (defined on the
    host process only) array, of type \c int.
   */
  int * GetINFOG() ;

  //! Copies the input array (of size 40) into the internally stored ICNTL array.
  int SetICNTL(int * ictnl);

  //! Set ICNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetICNTL(int pos, int value);

  //! Copy the input array (of size 5) into the internally stored CNTL array.
  int SetCNTL(double * ctnl);

  //! Set CNTL[pos] to value. pos is expressed in FORTRAN style (starting from 1).
  int SetCNTL(int pos, double value);

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
  
  Teuchos::RefCountPtr<Amesos_Mumps_Pimpl> PrivateMumpsData_; 

  //! Returns a reference to the linear system matrix.
  Epetra_RowMatrix& Matrix();

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
  void CheckError();

  //! Check input parameters.
  void CheckParameters();
  
  void SetICNTLandCNTL();

  //! \c true if matrix has already been converted to COO format
  bool IsConvertToTripletOK_;
  //! \c true if the Schur complement has been computed (need to free memory)
  bool IsComputeSchurComplementOK_;


  bool NoDestroy_ ;  // Set true to prevent memory freeing
  
  //! row indices of nonzero elements
  vector <int> Row;
  //! column indices of nonzero elements
  vector<int> Col;
  //! values of nonzero elements
  vector<double> Val;

#ifdef HAVE_AMESOS_SMUMPS
  //! single-precision values of nonzero elements
  vector<float> SVal;
  //! single-precision solution vector (on host only)
  vector<float> SVector;
#endif
  
  //! Maximum number of processors in the MUMPS' communicator
  int MaxProcs_;
  
  //! If \c true, solve the problem with AT.
  bool UseTranspose_;
  
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

  //! Pointer to the linear problem to be solved.
  const Epetra_LinearProblem* Problem_;

  //! Redistributed matrix.
  Epetra_Map* RedistrMap_;
  //! Redistributed importer (from Matrix().RowMatrixRowMap() to RedistrMatrix().RowMatrixRowMap()).
  Epetra_Import* RedistrImporter_;
  //! Redistributed matrix (only if MaxProcs_ > 1).
  Epetra_CrsMatrix* RedistrMatrix_;
  //! Map with all elements on process 0 (for solution and rhs).
  Epetra_Map* SerialMap_;
  //! Importer from Matrix.OperatorDomainMap() to SerialMap_.
  Epetra_Import* SerialImporter_;

#ifdef EPETRA_MPI
  //! MPI communicator used by MUMPS
  MPI_Comm MUMPSComm_;
#endif
  
};  // class Amesos_Mumps

#endif /* AMESOS_MUMPS_H */
