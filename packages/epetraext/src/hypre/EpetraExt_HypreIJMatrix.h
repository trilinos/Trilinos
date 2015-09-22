//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_HYPREIJMATRIX_H_
#define EPETRAEXT_HYPREIJMATRIX_H_

// Trilinos source files
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"

//Hypre source files
#include "krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE.h"

class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Import;

//! EpetraExt_HypreIJMatrix: A class for constructing and using real-valued sparse compressed row matrices.

/*! The EpetraExt_HypreIJMatrix is a wrapper class for Hypre parallel IJ matrices.  It is
    derived from the Epetra_BasicRowMatrix class, and so provides Hypre users access to Trilinos solvers.
    This class is lightweight, i.e., there are no deep copies of matrix data.  Whenever possible, class
    methods utilize callbacks to native Hypre functions.  
*/    

#ifndef HYPRE_ENUMS
#define HYPRE_ENUMS
//! Enumerated type for Hypre solvers.
/*! This enumerated type is used to determine which functions are used to create, destroy, setup and solve Hypre solvers.
*/
  enum Hypre_Solver{ 
    BoomerAMG, /*! < A BoomerAMG solver and preconditioner*/
    ParaSails, /*! < A ParaSails preconditioner */
    Euclid,    /*! < A Euclid preconditioner */
    AMS,       /*! < An AMS solver and preconditioner */
    Hybrid,    /*! < A Hybrid solver */
    PCG,       /*! < A PCG solver */
    GMRES,     /*! < A GMRES solver */
    FlexGMRES, /*! < A FlexGMRES solver */
    LGMRES,    /*! < A LGMRES solver */
    BiCGSTAB   /*! < A BiCGSTAB solver */
  };

//! Enumerated type to choose to solve or precondition
  enum Hypre_Chooser{
    Solver,        /*! < Choose to solve the system */
    Preconditioner /*! < Choose to apply preconditioner */
  };
#endif //HYPRE_ENUMS

class EpetraExt_HypreIJMatrix: public Epetra_BasicRowMatrix  {
      
public:

  //! @name Constructors/Destructor
  //@{ 
  //! Epetra_HypreIJMatrix constructor.
  /*! Creates a Epetra_HypreIJMatrix object by encapsulating an existing Hypre matrix.
    
      \param matrix (In) - A completely constructed Hypre IJ matrix.
  */
  EpetraExt_HypreIJMatrix(HYPRE_IJMatrix matrix);

  //! EpetraExt_HypreIJMatrix Destructor
  virtual ~EpetraExt_HypreIJMatrix();
  //@}
  
  //! @name Extraction methods
  //@{ 

  //! Returns a copy of the specified local row in user-provided arrays.
  /*! 
   \param MyRow (In) - Local row to extract.
   \param Length (In) - Length of Values and Indices.
   \param NumEntries (Out) - Number of nonzero entries extracted.
   \param Values (Out) - Extracted values for this row.
   \param Indices (Out) - Extracted local column indices for the corresponding values.
	  
   \return Integer error code, set to 0 if successful.
  */
  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;
    
  //! Returns a reference to the ith entry in the matrix, along with its row and column index.
  /*! 
    \param CurEntry (In) - Index of local entry (from 0 to NumMyNonzeros()-1) to extract. 
    \param Value (Out) - Extracted reference to current values. 
    \param RowIndex (Out) - Row index for current entry. 
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful.
  */
   
  int ExtractMyEntryView(int CurEntry, double *&Value, int &RowIndex, int &ColIndex);
    
  //! Returns a const reference to the ith entry in the matrix, along with its row and column index.
  /*! 
  \param CurEntry (In) - Index of local entry (from 0 to NumMyNonzeros()-1) to extract. 
  \param Value (Out) - Extracted reference to current values.
  \param RowIndex (Out) - Row index for current entry. 
  \param ColIndex (Out) - Column index for current entry.
	  
  \return Integer error code, set to 0 if successful.
  */
   
  int ExtractMyEntryView(int CurEntry, const double *&Value, int &RowIndex, int &ColIndex) const;
  //@}

  //! @name Solver Setup Methods
  //@{

  //! Set a parameter that takes a single int.
   /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set maximum iterations would be &HYPRE_BoomerAMGSetMaxIter
   \param parameter (In) -The integer parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter);

  //! Set a parameter that takes a single double.
  /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set tolerance would be &HYPRE_BoomerAMGSetTol
   \param parameter (In) -The double parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter);

  //! Set a parameter that takes a double then an int.
  /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set relaxation weight for a given level would be &HYPRE_BoomerAMGSetLevelRelaxWt
   \param parameter1 (In) -The double parameter being set.
   \param parameter2 (In) - The integer parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2);

  //! Set a parameter that takes two int parameters.
  /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set relaxation type for a given level would be &HYPRE_BoomerAMGSetCycleRelaxType
   \param parameter1 (In) -The first integer parameter being set.
   \param parameter2 (In) - The second integer parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2);

  //! Set a parameter that takes a double*.
  /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set relaxation weight would be &HYPRE_BoomerAMGSetRelaxWeight
   \param parameter (In) -The double* parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter);

  //! Set a parameter that takes an int*.
  /*!
   \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
   \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
     An example is if the solver is BoomerAMG, the function to set grid relax type would be &HYPRE_BoomerAMGSetGridRelaxType
   \param parameter (In) -The int* parameter being set.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter);

  //! Sets the solver that is used by the Solve() and ApplyInverse() methods. Until this is called, the default solver is PCG.
  /*! 
  \param chooser (In) - A Hypre_Chooser enumerated type. If Solver, then we are selecting which solver, if Preconditioner, we are choosing which preconditioner to use.
  \param Solver (In) -A Hypre_Solver enumerated type to select the solver or preconditioner. Options for solver are:
   BoomerAMG, AMS, Hybrid, PCG, GMRES, FlexGMRES, LGMRES, and BiCGSTAB. See Hypre Ref Manual for more info on the solvers.
   Options for Preconditioner are: BoomerAMG, ParaSails, Euclid, and AMS.
  \param transpose (In) -Optional argument that selects to use a transpose solve. Currently BoomerAMG is the only solver available for transpose solve. It must be the argument for Solver if transpose is true.

  \return Integer error code, set to 0 if successful.
  */

  int SetParameter(Hypre_Chooser chooser, Hypre_Solver Solver, bool transpose=false);

  //! Sets the solver to use the selected preconditioner.
  /*! 
   \param UsePreconditioner (In) -A boolean, true use preconditioner, false do not use the supplied preconditioner with the solver.
   The solver and preconditioner must have been selected and the solver must be one of the following solvers:
     Hybrid, PCG, GMRES, FlexGMRES, LGMRES, BiCGSTAB.

   \return Integer error code, set to 0 if successful.
  */

  int SetParameter(bool UsePreconditioner);

  //! Choose to solve the problem or apply the preconditioner.
  /*! 
   \param chooser (In) -A Hypre_Chooser enumerated type, either Solver or Preconditioner.
   The chosen type must have been selected before this method is called.

   \return Integer error code, set to 0 if successful.
  */
  int SetParameter(Hypre_Chooser chooser);
  //@}
  //! @name Computational methods
  //@{ 

  //! Returns the result of a EpetraExt_HypreIJMatrix multiplied by a Epetra_MultiVector X in Y.
  /*! 
   \param TransA (In) -If true, multiply by the transpose of matrix, otherwise just use matrix.
   \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
   \param Y (Out) -A Epetra_MultiVector of dimension NumVectorscontaining result.

   \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a EpetraExt_HypreIJMatrix solving a Epetra_MultiVector X in Y.
  /*! 
   \param Upper (In) -If true, solve Ux = y, otherwise solve Lx = y.
   \param Trans (In) -If true, solve transpose problem.
   \param UnitDiagonal (In) -If true, assume diagonal is unit (whether it's stored or not).
   \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
   \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

   \return Integer error code, set to 0 if successful.
  */

  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Scales the EpetraExt_HypreIJMatrix on the left with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
      and j denotes the column number of A.
   \param X (In) -A Epetra_Vector to solve for.

   \return Integer error code, set to 0 if successful.
  */
  int LeftScale(const Epetra_Vector& X);


  //! Scales the EpetraExt_HypreIJMatrix on the right with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.
    \param X (In) -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
  int RightScale(const Epetra_Vector& X);
  //@}

  //! @name Additional methods required to support the Epetra_Operator interface
  //@{ 

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
   \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
   \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

   \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(EpetraExt_HypreIJMatrix::Multiply(EpetraExt_HypreIJMatrix::UseTranspose(), X, Y));};

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the EpetraExt_HypreIJMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

   \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
   \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

   \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(EpetraExt_HypreIJMatrix::Solve(EpetraExt_HypreIJMatrix::UpperTriangular(), EpetraExt_HypreIJMatrix::UseTranspose(), false, X, Y));};

    //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const {return(false);}

  //@}
  //! @name Additional methods required to implement RowMatrix interface
  //@{ 

  //! Return the current number of values stored for the specified local row.
  /*! Similar to NumMyEntries() except NumEntries is returned as an argument
      and error checking is done on the input value MyRow.
  \param MyRow (In) - Local row.
  \param NumEntries (Out) - Number of nonzero values.
	  
  \return Integer error code, set to 0 if successful.
  */
  int NumMyRowEntries(int MyRow, int & NumEntries) const;

  //! Return a reference to the Hypre matrix.
  HYPRE_IJMatrix& GetMatrix(){ return Matrix_;}

  //@}
protected:

  //! Set global variables to default values.
  int InitializeDefaults();
   
  // These methods are needed only because the create methods in Hypre sometimes take an MPI_Comm but not always. 
  // They simply call the create solver in the correct way.
  //! AMG Create passing function
  int Hypre_BoomerAMGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_BoomerAMGCreate(solver);}

  //! ParaSails Create passing function
  int Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParaSailsCreate(comm, solver);}

  //! Euclid Create passing function
  int Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_EuclidCreate(comm, solver);}

  //! AMS Create passing function
  int Hypre_AMSCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_AMSCreate(solver);}

  //! Hybrid Create passing function
  int Hypre_ParCSRHybridCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRHybridCreate(solver);}

  //! PCG Create passing function
  int Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRPCGCreate(comm, solver);}

  //! GMRES Create passing function
  int Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRGMRESCreate(comm, solver);}

  //! FlexGMRES Create passing function
  int Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRFlexGMRESCreate(comm, solver);}

  //! LGMRES Create passing function
  int Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRLGMRESCreate(comm, solver);}

  //! BiCGSTAB Create passing function
  int Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRBiCGSTABCreate(comm, solver);}

  //! Create the solver with selected type.
  int CreateSolver();
  //! Create the preconditioner with selected type.
  int CreatePrecond();
  // These two methods setup the solver or preconditioner by calling the pointer. 
  // They are const because they are called from the Solve() routine.
  // They really aren't const because they change the value of IsSolverSetup_. This is because it should only be called if it isn't setup.
  int SetupSolver() const;
  int SetupPrecond() const;

  // Hypre variables
  mutable HYPRE_IJMatrix Matrix_;
  mutable HYPRE_ParCSRMatrix ParMatrix_;
  mutable HYPRE_IJVector X_hypre;
  mutable HYPRE_IJVector Y_hypre;
  mutable HYPRE_ParVector par_x;
  mutable HYPRE_ParVector par_y;
  mutable hypre_ParVector *x_vec;
  mutable hypre_ParVector *y_vec;
  mutable hypre_Vector *x_local;
  mutable hypre_Vector *y_local;
  mutable HYPRE_Solver Solver_;
  mutable HYPRE_Solver Preconditioner_;
  // The following are pointers to functions to use the solver and preconditioner.
  int (EpetraExt_HypreIJMatrix::*SolverCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*SolverDestroyPtr_)(HYPRE_Solver);
  int (*SolverSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverPrecondPtr_)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  int (EpetraExt_HypreIJMatrix::*PrecondCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*PrecondDestroyPtr_)(HYPRE_Solver);
  int (*PrecondSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*PrecondSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
    
  //! Number of rows on local processor.
  int NumMyRows_;
  //! Number of rows across all processors.
  int NumGlobalRows_;
  //! Number of columns across all processors.
  int NumGlobalCols_;
  //! First row on local processor.
  int MyRowStart_;
  //! Last row on local processor.
  int MyRowEnd_;
  //! Hypre matrix type (parCSR).
  mutable int MatType_;
  //! Do a transpose solve, only BoomerAMG.
  bool TransposeSolve_;
  //! Choose to solve or apply preconditioner.
  Hypre_Chooser SolveOrPrec_; 
  //! Flag to know if solver needs to be destoyed.
  bool *IsSolverSetup_; 
  //! Flag to know if preconditioner needs to be destroyed.
  bool *IsPrecondSetup_;
};
#endif /* EPETRAEXT_HYPREIJMATRIX_H_ */
